"""Python side of the Section 5 Python/Julia hull bridge.

The transport is intentionally dependency-free on the Julia side: point and
hull coordinates use raw little-endian Float64 files, while a small TSV
manifest carries job metadata. Binary transport avoids CSV quoting bugs and
keeps the 100-million-point experiment practical.
"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
import tempfile
import time
from collections.abc import Sequence

import numpy as np


REQUIRED_JOB_KEYS = {"id", "points", "k", "mode"}


def _validate_text_field(value, *, name: str) -> str:
    value = str(value)
    if not value or any(character in value for character in "\t\r\n"):
        raise ValueError(f"{name} must be a non-empty string without tabs/newlines")
    return value


def _normalize_jobs(jobs: Sequence[dict]) -> list[dict]:
    normalized = []
    seen_ids = set()
    group_settings = {}
    for index, original in enumerate(jobs):
        missing = REQUIRED_JOB_KEYS - set(original)
        if missing:
            raise ValueError(f"job {index} is missing keys: {sorted(missing)}")
        job_id = _validate_text_field(original["id"], name=f"job {index} id")
        if job_id in seen_ids:
            raise ValueError(f"duplicate job id: {job_id!r}")
        seen_ids.add(job_id)
        group = _validate_text_field(original.get("group", job_id), name=f"job {job_id} group")
        points = np.asarray(original["points"], dtype=np.float64)
        if points.ndim != 2 or points.shape[1] != 2:
            raise ValueError(f"job {job_id}: points must have shape (n, 2), got {points.shape}")
        if not np.all(np.isfinite(points)):
            raise ValueError(f"job {job_id}: points contain NaN or infinity")
        k = int(original["k"])
        if k < 1:
            raise ValueError(f"job {job_id}: k must be >= 1")
        mode = str(original["mode"])
        if mode not in {"seq", "par"}:
            raise ValueError(f"job {job_id}: mode must be 'seq' or 'par'")
        time_it = bool(original.get("time_it", False))
        repeats = int(original.get("repeats", 1))
        if repeats < 1:
            raise ValueError(f"job {job_id}: repeats must be >= 1")
        if time_it:
            setting = (mode, repeats)
            previous = group_settings.setdefault(group, setting)
            if previous != setting:
                raise ValueError(
                    f"timed group {group!r} mixes mode/repeat settings: {previous} vs {setting}"
                )
        normalized.append(
            {
                "id": job_id,
                "group": group,
                "points": points,
                "points_source": original["points"],
                "k": k,
                "mode": mode,
                "time_it": time_it,
                "repeats": repeats,
            }
        )
    return normalized


def _canonical_vertices(vertices) -> np.ndarray:
    vertices = np.asarray(vertices, dtype=np.float64).reshape(-1, 2)
    if len(vertices) > 1 and np.array_equal(vertices[0], vertices[-1]):
        vertices = vertices[:-1]
    if len(vertices) == 0:
        return vertices
    vertices = np.unique(vertices, axis=0)
    order = np.lexsort((vertices[:, 1], vertices[:, 0]))
    return vertices[order]


def _reference_hull(points) -> np.ndarray:
    """Qhull reference, with a deterministic fallback for degenerate data."""
    from scipy.spatial import ConvexHull, QhullError

    points = np.unique(np.asarray(points, dtype=np.float64), axis=0)
    if len(points) <= 2:
        return points
    try:
        return points[ConvexHull(points).vertices]
    except QhullError:
        span = np.ptp(points, axis=0)
        axis = int(np.argmax(span))
        order = np.argsort(points[:, axis], kind="stable")
        return points[order[[0, -1]]]


def assert_hulls_match_scipy(jobs: Sequence[dict], results: dict) -> None:
    """Require the returned convex region to equal SciPy/Qhull's region.

    The original 4k implementation can retain a numerically collinear point on
    an edge after feature standardization, while Qhull reports only the two
    edge endpoints.  Re-hulling the returned vertices removes only those
    redundant boundary points; missing or incorrect extreme vertices still
    fail this check.
    """
    references = {}
    for job in jobs:
        source_key = id(job["points_source"])
        if source_key not in references:
            references[source_key] = _canonical_vertices(_reference_hull(job["points"]))
        expected = references[source_key]
        actual = _canonical_vertices(_reference_hull(results[job["id"]]["hull"]))
        if not np.array_equal(actual, expected):
            raise AssertionError(
                f"backend/SciPy convex-hull mismatch for job {job['id']!r}: "
                f"backend has {len(actual)} vertices, SciPy has {len(expected)}"
            )


class HullClientBase:
    def run_batch(self, jobs, workdir=None):
        raise NotImplementedError


class JuliaHullClient(HullClientBase):
    """Call the original Julia sequential or threaded hull implementation."""

    def __init__(
        self,
        bridge_dir=None,
        julia_bin="julia",
        threads=2,
        verbose=True,
        verify_with_scipy=False,
    ):
        self.bridge_dir = bridge_dir or os.path.dirname(os.path.abspath(__file__))
        self.julia_bin = julia_bin
        self.threads = str(threads)
        self.verbose = bool(verbose)
        self.verify_with_scipy = bool(verify_with_scipy)

    def run_batch(self, jobs, workdir=None):
        jobs = _normalize_jobs(jobs)
        if not jobs:
            return {}
        if workdir is None:
            with tempfile.TemporaryDirectory(prefix="hull_bridge_") as temporary_directory:
                return self._run_batch(jobs, temporary_directory)
        os.makedirs(workdir, exist_ok=True)
        return self._run_batch(jobs, os.path.abspath(workdir))

    def _run_batch(self, jobs, workdir):
        points_path = os.path.join(workdir, "points.f64")
        manifest_path = os.path.join(workdir, "manifest.tsv")
        hulls_path = os.path.join(workdir, "hulls.f64")
        meta_path = os.path.join(workdir, "results.tsv")

        # Reuse a point range when sequential and parallel jobs reference the
        # same ndarray object.  This halves bridge storage for the synthetic
        # benchmark without relying on an expensive content hash.
        source_ranges = {}
        next_row = 1
        with open(points_path, "wb") as points_file:
            for job in jobs:
                source_key = id(job["points_source"])
                cached = source_ranges.get(source_key)
                if cached is not None and cached[0] is job["points_source"]:
                    job["start_row"], job["end_row"] = cached[1], cached[2]
                    continue
                points = np.ascontiguousarray(job["points"], dtype="<f8")
                start_row = next_row
                points.tofile(points_file)
                next_row += len(points)
                end_row = next_row - 1
                job["start_row"], job["end_row"] = start_row, end_row
                source_ranges[source_key] = (job["points_source"], start_row, end_row)

        with open(manifest_path, "w", newline="", encoding="utf-8") as manifest_file:
            writer = csv.writer(manifest_file, delimiter="\t", lineterminator="\n")
            writer.writerow(
                ["id", "group", "start_row", "end_row", "k", "mode", "time_it", "repeats"]
            )
            for job in jobs:
                writer.writerow(
                    [
                        job["id"],
                        job["group"],
                        job["start_row"],
                        job["end_row"],
                        job["k"],
                        job["mode"],
                        int(job["time_it"]),
                        job["repeats"],
                    ]
                )

        bridge_script = os.path.join(self.bridge_dir, "hull_bridge.jl")
        command = [
            self.julia_bin,
            "--threads",
            self.threads,
            bridge_script,
            points_path,
            manifest_path,
            hulls_path,
            meta_path,
        ]
        if self.verbose:
            print("  [JuliaHullClient]", shlex.join(command), flush=True)
        try:
            subprocess.run(command, check=True)
        except FileNotFoundError as error:
            raise RuntimeError(
                f"Julia executable not found: {self.julia_bin!r}; pass --julia-bin with its full path"
            ) from error

        results = _read_binary_results(hulls_path, meta_path)
        _validate_results(jobs, results)
        if self.verify_with_scipy:
            assert_hulls_match_scipy(jobs, results)
            if self.verbose:
                print(f"  [JuliaHullClient] {len(jobs)} hull region(s) match SciPy/Qhull")
        return results


class MockHullClient(HullClientBase):
    """Pure-Python port for smoke tests only; its timings are not reportable."""

    def __init__(self, verbose=True, verify_with_scipy=False):
        from eightk_original_port import alg_8ksided4ch

        self._algorithm = alg_8ksided4ch
        self.verbose = bool(verbose)
        self.verify_with_scipy = bool(verify_with_scipy)

    def _run_hull(self, points, k):
        points = np.asarray(points, dtype=np.float64)
        if len(points) <= 2:
            return np.unique(points, axis=0)
        raw = self._algorithm(points, k)
        return np.asarray([np.asarray(point) for point in raw], dtype=np.float64).reshape(-1, 2)

    def run_batch(self, jobs, workdir=None):
        del workdir
        jobs = _normalize_jobs(jobs)
        results = {}
        timed_groups = {}
        for job in jobs:
            if job["time_it"]:
                timed_groups.setdefault(job["group"], []).append(job)
            else:
                hull = self._run_hull(job["points"], job["k"])
                results[job["id"]] = _result(hull, len(job["points"]), [])

        for group_jobs in timed_groups.values():
            runs = []
            final_hulls = {}
            for _ in range(group_jobs[0]["repeats"]):
                start = time.perf_counter()
                for job in group_jobs:
                    final_hulls[job["id"]] = self._run_hull(job["points"], job["k"])
                runs.append(time.perf_counter() - start)
            for job in group_jobs:
                results[job["id"]] = _result(
                    final_hulls[job["id"]], len(job["points"]), runs
                )

        _validate_results(jobs, results)
        if self.verify_with_scipy:
            assert_hulls_match_scipy(jobs, results)
        if self.verbose:
            print(
                f"  [MockHullClient] {len(jobs)} job(s); timings are smoke-test only, not reportable"
            )
        return results


def _result(hull, n_input, timing_runs):
    timing_runs = [float(value) for value in timing_runs]
    return {
        "hull": np.asarray(hull, dtype=np.float64).reshape(-1, 2),
        "time_s": float(np.median(timing_runs)) if timing_runs else None,
        "time_runs_s": timing_runs,
        "n_verts": len(hull),
        "n_input": int(n_input),
    }


def _read_binary_results(hulls_path, meta_path):
    coordinates = np.fromfile(hulls_path, dtype="<f8")
    if len(coordinates) % 2:
        raise RuntimeError("Julia hull output contains an odd number of Float64 values")
    vertices = coordinates.reshape(-1, 2)
    results = {}
    with open(meta_path, encoding="utf-8", newline="") as meta_file:
        for row in csv.DictReader(meta_file, delimiter="\t"):
            start = int(row["start_vertex"])
            end = int(row["end_vertex"])
            hull = vertices[start - 1 : end].copy() if end >= start else np.empty((0, 2))
            timing_runs = [
                float(value) for value in row["time_runs_s"].split(";") if value
            ]
            results[row["id"]] = _result(hull, int(row["n_input"]), timing_runs)
    return results


def _validate_results(jobs, results):
    expected_ids = {job["id"] for job in jobs}
    if set(results) != expected_ids:
        missing = sorted(expected_ids - set(results))
        extra = sorted(set(results) - expected_ids)
        raise RuntimeError(f"bridge result IDs differ; missing={missing}, extra={extra}")
    for job in jobs:
        result = results[job["id"]]
        hull = result["hull"]
        if hull.ndim != 2 or hull.shape[1] != 2 or not np.all(np.isfinite(hull)):
            raise RuntimeError(f"job {job['id']!r}: invalid hull array")
        if result["n_input"] != len(job["points"]):
            raise RuntimeError(f"job {job['id']!r}: bridge reported the wrong input size")
        if result["n_verts"] != len(hull):
            raise RuntimeError(f"job {job['id']!r}: bridge reported the wrong vertex count")
        expected_runs = job["repeats"] if job["time_it"] else 0
        if len(result["time_runs_s"]) != expected_runs:
            raise RuntimeError(
                f"job {job['id']!r}: expected {expected_runs} timing runs, "
                f"got {len(result['time_runs_s'])}"
            )


def make_client(backend="julia", **kwargs):
    if backend == "julia":
        return JuliaHullClient(**kwargs)
    if backend == "mock":
        return MockHullClient(**kwargs)
    raise ValueError(f"unknown backend: {backend!r}; expected 'julia' or 'mock'")


def _smoke_test():
    parser = argparse.ArgumentParser(description="Smoke-test the Julia hull bridge")
    parser.add_argument("--backend", choices=["julia", "mock"], default="julia")
    parser.add_argument("--julia-bin", default="julia")
    parser.add_argument("--threads", default="2")
    args = parser.parse_args()

    rng = np.random.default_rng(42)
    points = rng.normal(size=(1_000, 2))
    jobs = [
        {
            "id": f"smoke-{mode}",
            "group": f"smoke-{mode}",
            "points": points,
            "k": 1,
            "mode": mode,
            "time_it": True,
            "repeats": 3,
        }
        for mode in ("seq", "par")
    ]
    if args.backend == "julia":
        client = JuliaHullClient(
            julia_bin=args.julia_bin,
            threads=args.threads,
            verify_with_scipy=True,
        )
    else:
        client = MockHullClient(verify_with_scipy=True)
    results = client.run_batch(jobs)
    if not np.array_equal(
        _canonical_vertices(results["smoke-seq"]["hull"]),
        _canonical_vertices(results["smoke-par"]["hull"]),
    ):
        raise AssertionError("sequential and parallel Julia hulls differ")
    print("hull_client.py: bridge smoke test passed")


if __name__ == "__main__":
    _smoke_test()
