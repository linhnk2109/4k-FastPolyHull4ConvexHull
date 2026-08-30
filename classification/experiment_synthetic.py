"""Section 5.2 - large-scale three-class Gaussian-blob experiment.

``--profile paper`` uses the six nominal sizes in the updated Tables 17-18,
from 1,000 through 100,000,000, with equal class counts ``floor(n/3)``.
``--profile extended`` generates exactly the five configured totals through
10,000,000 points and is retained for compatibility with earlier reruns.
"""

from __future__ import annotations

import argparse
import csv
import os
import time

import numpy as np
from sklearn.datasets import make_blobs
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import LinearSVC, SVC

from classification_rule import HullClassificationRule
from hull_client import _canonical_vertices, make_client


HERE = os.path.dirname(os.path.abspath(__file__))
PAPER_SIZES = [1_000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000]
EXTENDED_SIZES = [1_000, 10_000, 100_000, 1_000_000, 10_000_000]
CENTERS = np.array([[0.0, 0.0], [6.0, 0.0], [3.0, 5.2]])
CLUSTER_STD = 2.0
RANDOM_STATE = 42
K_HULL = 1
DEFAULT_REPEATS = 3
DEFAULT_RBF_MAX_N = 20_000

METHOD_LABELS = {
    "knn": "k-NN (k=5)",
    "linear_svm": "Linear SVM",
    "rbf_svm": "RBF SVM",
    "hull_seq": "4k-FPH (sequential)",
    "hull_par": "P-4k-FPH (parallel)",
}

PAPER_ACCURACY_PCT = {
    2_000: {"knn": 85.00, "linear_svm": 85.83, "rbf_svm": 86.17, "hull_seq": 85.67, "hull_par": 85.67},
    20_000: {"knn": 87.45, "linear_svm": 88.78, "rbf_svm": 88.83, "hull_seq": 88.87, "hull_par": 88.87},
    200_000: {"knn": 86.79, "linear_svm": 88.43, "hull_seq": 88.42, "hull_par": 88.42},
    1_000_000: {"knn": 87.02, "linear_svm": 88.48, "hull_seq": 88.47, "hull_par": 88.47},
}

PAPER_TIME_S = {
    2_000: {"knn": 0.0028, "linear_svm": 0.0017, "rbf_svm": 0.0161, "hull_seq": 0.0061, "hull_par": 0.0277},
    20_000: {"knn": 0.0161, "linear_svm": 0.0108, "rbf_svm": 1.4875, "hull_seq": 0.0155, "hull_par": 0.0303},
    200_000: {"knn": 0.2311, "linear_svm": 0.1982, "hull_seq": 0.0741, "hull_par": 0.0883},
    1_000_000: {"knn": 1.4349, "linear_svm": 0.9473, "hull_seq": 0.3143, "hull_par": 0.2925},
}

COLORS = {
    "knn": "#eb6834",
    "linear_svm": "#1baf7a",
    "rbf_svm": "#eda100",
    "hull_seq": "#2a78d6",
    "hull_par": "#7149b8",
}


def generate_data(n, *, exact_total):
    if n < 3:
        raise ValueError("n must be at least 3")
    counts = [n // 3] * 3
    if exact_total:
        for index in range(n % 3):
            counts[index] += 1
    X, y = make_blobs(
        n_samples=counts,
        centers=CENTERS,
        cluster_std=CLUSTER_STD,
        random_state=RANDOM_STATE,
    )
    split = train_test_split(
        X,
        y,
        test_size=0.30,
        stratify=y,
        random_state=RANDOM_STATE,
    )
    return (*split, len(X))


def time_sklearn(model_factory, X_train, y_train, X_test, y_test, repeats):
    fit_runs, predict_runs, total_runs = [], [], []
    predictions = None
    for _ in range(repeats):
        model = model_factory()
        start = time.perf_counter()
        model.fit(X_train, y_train)
        after_fit = time.perf_counter()
        predictions = model.predict(X_test)
        after_predict = time.perf_counter()
        fit_runs.append(after_fit - start)
        predict_runs.append(after_predict - after_fit)
        total_runs.append(after_predict - start)
    return make_timing_result(fit_runs, predict_runs, total_runs, predictions, y_test)


def make_timing_result(fit_runs, predict_runs, total_runs, predictions, y_test):
    return {
        "fit_s": float(np.median(fit_runs)),
        "predict_s": float(np.median(predict_runs)),
        "total_s": float(np.median(total_runs)),
        "fit_runs_s": [float(value) for value in fit_runs],
        "predict_runs_s": [float(value) for value in predict_runs],
        "total_runs_s": [float(value) for value in total_runs],
        "accuracy": float(np.mean(predictions == y_test)),
    }


def build_hull_jobs(n, class_points, repeats):
    jobs = []
    for mode in ("seq", "par"):
        group = f"n-{n}__{mode}"
        for class_index, label in enumerate(sorted(class_points)):
            jobs.append(
                {
                    "id": f"{group}__class-{class_index}",
                    "group": group,
                    # The same ndarray object is shared by seq/par.  The binary
                    # bridge writes it once and lets both jobs reuse its range.
                    "points": class_points[label],
                    "k": K_HULL,
                    "mode": mode,
                    "time_it": True,
                    "repeats": repeats,
                }
            )
    return jobs


def time_centroids(class_points, repeats):
    runs = []
    centroids = None
    for _ in range(repeats):
        start = time.perf_counter()
        centroids = {label: points.mean(axis=0) for label, points in class_points.items()}
        runs.append(time.perf_counter() - start)
    return centroids, runs


def hull_result(
    n,
    mode,
    classes,
    class_points,
    backend_results,
    X_test,
    y_test,
    repeats,
    centroid_runs,
    centroids,
):
    group = f"n-{n}__{mode}"
    hulls = {
        label: backend_results[f"{group}__class-{class_index}"]["hull"]
        for class_index, label in enumerate(classes)
    }
    bridge_fit_runs = backend_results[f"{group}__class-0"]["time_runs_s"]
    if len(bridge_fit_runs) != repeats:
        raise RuntimeError(f"{group}: backend returned the wrong number of timing runs")
    fit_runs = [bridge_fit_runs[index] + centroid_runs[index] for index in range(repeats)]

    classifier = HullClassificationRule().fit_from_hulls(hulls, centroids)
    predict_runs = []
    predictions = None
    for _ in range(repeats):
        start = time.perf_counter()
        predictions = classifier.predict(X_test)
        predict_runs.append(time.perf_counter() - start)
    total_runs = [fit_runs[index] + predict_runs[index] for index in range(repeats)]
    return (
        make_timing_result(fit_runs, predict_runs, total_runs, predictions, y_test),
        hulls,
        predictions,
    )


def add_record(records, method, n, n_generated, result, backend_name):
    paper_accuracy = PAPER_ACCURACY_PCT.get(n, {}).get(method)
    paper_time = PAPER_TIME_S.get(n, {}).get(method)
    if result is None:
        records.append(
            {
                "method": method,
                "method_label": METHOD_LABELS[method],
                "n_total": n,
                "n_generated": n_generated,
                "fit_s": "",
                "predict_s": "",
                "total_s": "",
                "accuracy": "",
                "accuracy_pct": "",
                "paper_accuracy_pct": paper_accuracy if paper_accuracy is not None else "",
                "accuracy_delta_pp": "",
                "paper_total_s": paper_time if paper_time is not None else "",
                "timing_reportable": False,
                "fit_runs_s": "",
                "predict_runs_s": "",
                "total_runs_s": "",
            }
        )
        return

    accuracy_pct = 100.0 * result["accuracy"]
    records.append(
        {
            "method": method,
            "method_label": METHOD_LABELS[method],
            "n_total": n,
            "n_generated": n_generated,
            "fit_s": result["fit_s"],
            "predict_s": result["predict_s"],
            "total_s": result["total_s"],
            "accuracy": result["accuracy"],
            "accuracy_pct": accuracy_pct,
            "paper_accuracy_pct": paper_accuracy if paper_accuracy is not None else "",
            "accuracy_delta_pp": accuracy_pct - paper_accuracy if paper_accuracy is not None else "",
            "paper_total_s": paper_time if paper_time is not None else "",
            "timing_reportable": not (backend_name == "mock" and method.startswith("hull_")),
            "fit_runs_s": ";".join(f"{value:.9f}" for value in result["fit_runs_s"]),
            "predict_runs_s": ";".join(f"{value:.9f}" for value in result["predict_runs_s"]),
            "total_runs_s": ";".join(f"{value:.9f}" for value in result["total_runs_s"]),
        }
    )


def run_experiment(sizes, client, backend_name, repeats, rbf_max_n, profile):
    per_size = {}
    hull_jobs = []
    for n in sizes:
        # The current manuscript tables were generated with floor(n/3)
        # samples per class, so their actual total is 3*floor(n/3). Preserve
        # that behavior only in the explicit paper-reproduction profile; the
        # extended profile generates exactly the configured n.
        X_train, X_test, y_train, y_test, n_generated = generate_data(
            n, exact_total=profile != "paper"
        )
        classes = np.unique(y_train).tolist()
        class_points = {
            label: np.ascontiguousarray(X_train[y_train == label]) for label in classes
        }
        per_size[n] = {
            "X_train": X_train,
            "X_test": X_test,
            "y_train": y_train,
            "y_test": y_test,
            "classes": classes,
            "class_points": class_points,
            "n_generated": n_generated,
        }
        hull_jobs.extend(build_hull_jobs(n, class_points, repeats))

    print(
        f"Calling {backend_name} hull backend once: {len(hull_jobs)} jobs, "
        f"{repeats} timing runs/group"
    )
    wall_start = time.perf_counter()
    backend_results = client.run_batch(hull_jobs)
    print(f"Backend subprocess/batch wall time (not a table value): {time.perf_counter() - wall_start:.3f}s")

    records = []
    for n in sizes:
        data = per_size[n]
        X_train, X_test = data["X_train"], data["X_test"]
        y_train, y_test = data["y_train"], data["y_test"]
        print(
            f"\n--- n label={n:,}; generated={data['n_generated']:,}; "
            f"train={len(y_train):,}; test={len(y_test):,} ---"
        )

        knn = time_sklearn(
            lambda: KNeighborsClassifier(n_neighbors=5, n_jobs=1),
            X_train,
            y_train,
            X_test,
            y_test,
            repeats,
        )
        add_record(records, "knn", n, data["n_generated"], knn, backend_name)

        linear_svm = time_sklearn(
            lambda: LinearSVC(max_iter=20_000, random_state=RANDOM_STATE),
            X_train,
            y_train,
            X_test,
            y_test,
            repeats,
        )
        add_record(records, "linear_svm", n, data["n_generated"], linear_svm, backend_name)

        if n <= rbf_max_n:
            rbf_svm = time_sklearn(
                lambda: SVC(kernel="rbf"),
                X_train,
                y_train,
                X_test,
                y_test,
                repeats,
            )
            add_record(records, "rbf_svm", n, data["n_generated"], rbf_svm, backend_name)
        else:
            rbf_svm = None
            add_record(records, "rbf_svm", n, data["n_generated"], None, backend_name)

        centroids, centroid_runs = time_centroids(data["class_points"], repeats)
        sequential, sequential_hulls, sequential_predictions = hull_result(
            n,
            "seq",
            data["classes"],
            data["class_points"],
            backend_results,
            X_test,
            y_test,
            repeats,
            centroid_runs,
            centroids,
        )
        parallel, parallel_hulls, parallel_predictions = hull_result(
            n,
            "par",
            data["classes"],
            data["class_points"],
            backend_results,
            X_test,
            y_test,
            repeats,
            centroid_runs,
            centroids,
        )
        for label in data["classes"]:
            if not np.array_equal(
                _canonical_vertices(sequential_hulls[label]),
                _canonical_vertices(parallel_hulls[label]),
            ):
                raise AssertionError(f"n={n}, class={label}: sequential/parallel hull vertices differ")
        if not np.array_equal(sequential_predictions, parallel_predictions):
            raise AssertionError(f"n={n}: sequential/parallel hull predictions differ")
        add_record(records, "hull_seq", n, data["n_generated"], sequential, backend_name)
        add_record(records, "hull_par", n, data["n_generated"], parallel, backend_name)
        data["hulls"] = sequential_hulls

        for method, result in (
            ("knn", knn),
            ("linear_svm", linear_svm),
            ("rbf_svm", rbf_svm),
            ("hull_seq", sequential),
            ("hull_par", parallel),
        ):
            if result is None:
                print(f"  {METHOD_LABELS[method]:24s} skipped (n > {rbf_max_n:,})")
            else:
                print(
                    f"  {METHOD_LABELS[method]:24s} "
                    f"acc={100 * result['accuracy']:.2f}% total={result['total_s']:.6f}s"
                )
    return records, per_size


def write_long_csv(records, output_path):
    with open(output_path, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def write_wide_tables(records, results_dir, sizes):
    lookup = {(row["n_total"], row["method"]): row for row in records}
    methods = list(METHOD_LABELS)
    for value_key, filename in (
        ("accuracy_pct", "synthetic_accuracy_table.csv"),
        ("total_s", "synthetic_runtime_table.csv"),
    ):
        path = os.path.join(results_dir, filename)
        with open(path, "w", newline="", encoding="utf-8") as output:
            writer = csv.writer(output)
            writer.writerow(["n_total", *[METHOD_LABELS[method] for method in methods]])
            for n in sizes:
                writer.writerow([n, *[lookup[(n, method)][value_key] for method in methods]])

    comparison_path = os.path.join(results_dir, "synthetic_paper_comparison.csv")
    comparison_fields = [
        "n_total",
        "method",
        "measured_accuracy_pct",
        "paper_accuracy_pct",
        "accuracy_delta_pp",
        "rounded_accuracy_matches",
        "measured_total_s",
        "paper_total_s",
    ]
    with open(comparison_path, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=comparison_fields)
        writer.writeheader()
        for row in records:
            if row["paper_accuracy_pct"] == "":
                continue
            measured = row["accuracy_pct"]
            expected = row["paper_accuracy_pct"]
            writer.writerow(
                {
                    "n_total": row["n_total"],
                    "method": row["method_label"],
                    "measured_accuracy_pct": measured,
                    "paper_accuracy_pct": expected,
                    "accuracy_delta_pp": measured - expected if measured != "" else "",
                    "rounded_accuracy_matches": measured != "" and round(measured, 2) == round(expected, 2),
                    "measured_total_s": row["total_s"],
                    "paper_total_s": row["paper_total_s"],
                }
            )


def make_figures(records, per_size, sizes, figures_dir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    lookup = {(row["n_total"], row["method"]): row for row in records}
    markers = {"knn": "s", "linear_svm": "^", "rbf_svm": "D", "hull_seq": "o", "hull_par": "v"}

    def accuracy_panel(axis):
        for method in ("knn", "linear_svm", "rbf_svm", "hull_seq"):
            xs = [n for n in sizes if lookup[(n, method)]["accuracy_pct"] != ""]
            ys = [lookup[(n, method)]["accuracy_pct"] for n in xs]
            axis.plot(xs, ys, marker=markers[method], color=COLORS[method], label=METHOD_LABELS[method], linewidth=1.8)
        axis.set_xscale("log")
        axis.set_xlabel("Number of samples n (log scale)")
        axis.set_ylabel("Held-out accuracy (%)")
        axis.set_title("(c) Accuracy", loc="left", fontweight="bold")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=8)

    def runtime_panel(axis):
        for method in METHOD_LABELS:
            xs = [n for n in sizes if lookup[(n, method)]["total_s"] != ""]
            ys = [lookup[(n, method)]["total_s"] for n in xs]
            axis.plot(
                xs,
                ys,
                marker=markers[method],
                color=COLORS[method],
                label=METHOD_LABELS[method],
                linewidth=1.8,
                linestyle="--" if method == "hull_par" else "-",
            )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("Number of samples n (log scale)")
        axis.set_ylabel("Fit + predict time (s; median of runs)")
        axis.set_title("(d) Runtime", loc="left", fontweight="bold")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=8)

    def data_panel(axis, n):
        data = per_size[n]
        colors = ["#2a78d6", "#eb6834", "#1baf7a"]
        for class_index, label in enumerate(data["classes"]):
            points = data["X_train"][data["y_train"] == label]
            axis.scatter(points[:, 0], points[:, 1], s=5, color=colors[class_index], alpha=0.45, linewidths=0)
            hull = data["hulls"][label]
            axis.add_patch(Polygon(hull, closed=True, fill=False, edgecolor=colors[class_index], linewidth=1.8))
        axis.set_title(
            f"(a) Gaussian blobs and class convex hulls (n={n:,})",
            loc="left",
            fontweight="bold",
        )
        axis.set_xlabel("x1")
        axis.set_ylabel("x2")

    def decision_panel(axis, n):
        data = per_size[n]
        centroids = {label: data["class_points"][label].mean(axis=0) for label in data["classes"]}
        classifier = HullClassificationRule().fit_from_hulls(data["hulls"], centroids)
        all_points = np.vstack([data["X_train"], data["X_test"]])
        padding = 2.0
        xx, yy = np.meshgrid(
            np.linspace(all_points[:, 0].min() - padding, all_points[:, 0].max() + padding, 220),
            np.linspace(all_points[:, 1].min() - padding, all_points[:, 1].max() + padding, 220),
        )
        grid = np.column_stack([xx.ravel(), yy.ravel()])
        prediction = np.asarray(classifier.predict(grid), dtype=float).reshape(xx.shape)
        color_map = matplotlib.colors.ListedColormap(["#2a78d6", "#eb6834", "#1baf7a"])
        axis.pcolormesh(xx, yy, prediction, cmap=color_map, alpha=0.18, shading="auto")
        axis.scatter(data["X_test"][:, 0], data["X_test"][:, 1], s=7, c=data["y_test"], cmap=color_map, linewidths=0)
        axis.set_title("(b) Convex-hull classifier decision regions", loc="left", fontweight="bold")
        axis.set_xlabel("x1")
        axis.set_ylabel("x2")

    display_n = min(sizes)
    figure, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=180)
    data_panel(axes[0, 0], display_n)
    decision_panel(axes[0, 1], display_n)
    accuracy_panel(axes[1, 0])
    runtime_panel(axes[1, 1])
    figure.suptitle("Section 5.2 - synthetic classification benchmark", fontsize=13, fontweight="bold")
    figure.tight_layout(rect=[0, 0, 1, 0.96])
    figure.savefig(os.path.join(figures_dir, "synthetic_overview.png"), bbox_inches="tight")
    plt.close(figure)

    for panel_function, filename in (
        (accuracy_panel, "synthetic_accuracy.png"),
        (runtime_panel, "synthetic_runtime.png"),
    ):
        figure, axis = plt.subplots(figsize=(8, 5.5), dpi=180)
        panel_function(axis)
        figure.tight_layout()
        figure.savefig(os.path.join(figures_dir, filename), bbox_inches="tight")
        plt.close(figure)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=["paper", "extended"], default="paper")
    parser.add_argument("--sizes", type=int, nargs="+", help="override the selected profile")
    parser.add_argument("--backend", choices=["julia", "mock"], default="julia")
    parser.add_argument("--julia-bin", default="julia")
    parser.add_argument("--julia-threads", default="4")
    parser.add_argument("--verify-hulls", action="store_true")
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS)
    parser.add_argument("--rbf-max-n", type=int, default=DEFAULT_RBF_MAX_N)
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--results-dir", default=os.path.join(HERE, "results"))
    parser.add_argument("--figures-dir", default=os.path.join(HERE, "figures"))
    return parser.parse_args()


def main():
    args = parse_args()
    if args.repeats < 1:
        raise SystemExit("--repeats must be >= 1")
    sizes = args.sizes or (PAPER_SIZES if args.profile == "paper" else EXTENDED_SIZES)
    if len(set(sizes)) != len(sizes) or any(n < 3 for n in sizes):
        raise SystemExit("sizes must be distinct integers >= 3")

    results_dir = os.path.abspath(args.results_dir)
    figures_dir = os.path.abspath(args.figures_dir)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    client_kwargs = {"verify_with_scipy": args.verify_hulls}
    if args.backend == "julia":
        client_kwargs.update(julia_bin=args.julia_bin, threads=args.julia_threads)
    client = make_client(args.backend, **client_kwargs)
    records, per_size = run_experiment(
        sizes,
        client,
        args.backend,
        args.repeats,
        args.rbf_max_n,
        args.profile,
    )

    long_path = os.path.join(results_dir, "synthetic_benchmark_results.csv")
    write_long_csv(records, long_path)
    write_wide_tables(records, results_dir, sizes)
    if not args.no_plots:
        make_figures(records, per_size, sizes, figures_dir)
    print("\nwrote", long_path)
    if all(n in PAPER_SIZES for n in sizes):
        matches = [
            row
            for row in records
            if row["paper_accuracy_pct"] != ""
            and row["accuracy_pct"] != ""
            and round(row["accuracy_pct"], 2) == round(row["paper_accuracy_pct"], 2)
        ]
        comparable = [
            row for row in records if row["paper_accuracy_pct"] != "" and row["accuracy_pct"] != ""
        ]
        print(f"paper accuracy check (rounded to 2 decimals): {len(matches)}/{len(comparable)} match")


if __name__ == "__main__":
    main()
