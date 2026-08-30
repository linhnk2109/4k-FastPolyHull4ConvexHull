"""Regenerate Section 5 figures from saved benchmark results.

The large synthetic timings are read from CSV. Only the deterministic
``n=1,000`` example is recomputed so that the illustrative convex hulls and
decision regions can be redrawn without repeating the 100-million-point run.
"""

from __future__ import annotations

import csv
import os

from experiment_synthetic import DEFAULT_RBF_MAX_N, make_figures, run_experiment
from hull_client import make_client


HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results", "synthetic_paper")
FIGURES_DIR = os.path.join(HERE, "figures", "synthetic_paper")


def load_saved_records(path):
    records = []
    with open(path, encoding="utf-8", newline="") as source:
        for row in csv.DictReader(source):
            row["n_total"] = int(row["n_total"])
            for key in ("accuracy_pct", "total_s"):
                row[key] = float(row[key]) if row[key] else ""
            records.append(row)
    if not records:
        raise RuntimeError(f"No saved benchmark rows in {path}")
    return records


def main():
    records_path = os.path.join(RESULTS_DIR, "synthetic_benchmark_results.csv")
    records = load_saved_records(records_path)
    sizes = sorted({row["n_total"] for row in records})
    display_n = min(sizes)

    client = make_client("julia", threads="4")
    _, per_size = run_experiment(
        [display_n],
        client,
        "julia",
        repeats=1,
        rbf_max_n=DEFAULT_RBF_MAX_N,
        profile="paper",
    )
    os.makedirs(FIGURES_DIR, exist_ok=True)
    make_figures(records, per_size, sizes, FIGURES_DIR)
    print("Refreshed paper figures in", FIGURES_DIR)


if __name__ == "__main__":
    main()
