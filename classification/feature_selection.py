"""Section 5.3 - datasets and exhaustive two-feature selection.

The default ``strict`` protocol follows the leakage-free workflow:
perform the 70/30 split first, then run 5-fold CV only on the 70% training
partition.  ``--protocol paper`` reproduces the feature-selection heatmaps and
numbers printed in the current manuscript, which ran CV on the complete
dataset before making the final 70/30 split.  Keeping both modes makes the
methodological difference explicit instead of silently mixing their outputs.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import os

import numpy as np
from sklearn.datasets import load_breast_cancer, load_iris, load_wine
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler

from classification_rule import HullClassificationRule
from hull_client import make_client


HERE = os.path.dirname(os.path.abspath(__file__))
RANDOM_STATE = 42
K_HULL = 1
N_SPLITS = 5

DATASETS = {
    "iris": load_iris,
    "wine": load_wine,
    "breast_cancer": load_breast_cancer,
}

PAPER_FEATURE_SELECTION = {
    "iris": (2, 3, 0.9666666666666667),
    "wine": (0, 6, 0.9095238095238095),
    "breast_cancer": (22, 24, 0.9595870206489675),
}


def _split(X, y):
    return train_test_split(
        X,
        y,
        test_size=0.30,
        stratify=y,
        random_state=RANDOM_STATE,
    )


def build_cv_jobs(X, y, pairs, n_splits=N_SPLITS):
    """Build every (pair, fold, class) hull as one backend batch."""
    classes = np.unique(y).tolist()
    folds = list(
        StratifiedKFold(
            n_splits=n_splits,
            shuffle=True,
            random_state=RANDOM_STATE,
        ).split(X, y)
    )

    jobs = []
    job_key = {}
    fold_data = {}
    for feature_i, feature_j in pairs:
        pair = (feature_i, feature_j)
        for fold_index, (train_indices, validation_indices) in enumerate(folds):
            train_pair = X[train_indices][:, [feature_i, feature_j]]
            scaler = StandardScaler().fit(train_pair)
            scaled_train = scaler.transform(train_pair)
            fold_labels = y[train_indices]
            fold_data[(pair, fold_index)] = {
                "scaler": scaler,
                "scaled_train": scaled_train,
                "train_labels": fold_labels,
                "validation_indices": validation_indices,
            }
            for class_index, label in enumerate(classes):
                job_id = (
                    f"pair-{feature_i}-{feature_j}__fold-{fold_index}__class-{class_index}"
                )
                jobs.append(
                    {
                        "id": job_id,
                        "group": job_id,
                        "points": scaled_train[fold_labels == label],
                        "k": K_HULL,
                        "mode": "seq",
                        "time_it": False,
                    }
                )
                job_key[job_id] = (pair, fold_index, label)
    return jobs, job_key, folds, classes, fold_data


def cv_accuracy_all_pairs(X, y, feature_names, client):
    pairs = list(itertools.combinations(range(X.shape[1]), 2))
    jobs, job_key, folds, classes, fold_data = build_cv_jobs(X, y, pairs)
    print(
        f"    {len(pairs)} pairs x {N_SPLITS} folds -> "
        f"{len(jobs)} hull jobs in one backend call"
    )
    backend_results = client.run_batch(jobs)

    hulls_by_pair_fold = {}
    for job_id, (pair, fold_index, label) in job_key.items():
        hulls_by_pair_fold.setdefault((pair, fold_index), {})[label] = backend_results[
            job_id
        ]["hull"]

    rows = []
    for feature_i, feature_j in pairs:
        pair = (feature_i, feature_j)
        fold_accuracies = []
        for fold_index, (_, validation_indices) in enumerate(folds):
            data = fold_data[(pair, fold_index)]
            validation = data["scaler"].transform(
                X[validation_indices][:, [feature_i, feature_j]]
            )
            centroids = {
                label: data["scaled_train"][data["train_labels"] == label].mean(axis=0)
                for label in classes
            }
            classifier = HullClassificationRule().fit_from_hulls(
                hulls_by_pair_fold[(pair, fold_index)], centroids
            )
            predictions = classifier.predict(validation)
            fold_accuracies.append(float(np.mean(predictions == y[validation_indices])))

        rows.append(
            {
                "feat_i": feature_i,
                "feat_j": feature_j,
                "name_i": str(feature_names[feature_i]),
                "name_j": str(feature_names[feature_j]),
                "cv_mean_acc": float(np.mean(fold_accuracies)),
                "cv_std_acc": float(np.std(fold_accuracies)),
                "fold_accuracies": ";".join(f"{value:.12f}" for value in fold_accuracies),
            }
        )

    # Numerical ties are resolved by lower fold-to-fold variability and then
    # feature index.  Rounding prevents a ~1e-16 summation artifact from
    # deciding the winner.
    rows.sort(
        key=lambda row: (
            -round(row["cv_mean_acc"], 12),
            round(row["cv_std_acc"], 12),
            row["feat_i"],
            row["feat_j"],
        )
    )
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows


def make_heatmap(dataset_name, rows, best, figures_dir, protocol):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    dimension = max(max(row["feat_i"], row["feat_j"]) for row in rows) + 1
    names = [None] * dimension
    matrix = np.full((dimension, dimension), np.nan)
    for row in rows:
        feature_i, feature_j = row["feat_i"], row["feat_j"]
        names[feature_i], names[feature_j] = row["name_i"], row["name_j"]
        matrix[feature_i, feature_j] = matrix[feature_j, feature_i] = row["cv_mean_acc"]

    figure_width = max(6.0, 0.28 * dimension + 3.0)
    figure, axis = plt.subplots(figsize=(figure_width, figure_width * 0.92), dpi=180)
    image = axis.imshow(matrix, cmap="viridis", vmin=np.nanmin(matrix), vmax=np.nanmax(matrix))
    font_size = 5.5 if dimension > 15 else 8
    axis.set_xticks(range(dimension), labels=names, rotation=90, fontsize=font_size)
    axis.set_yticks(range(dimension), labels=names, fontsize=font_size)
    axis.set_xticks(np.arange(-0.5, dimension, 1), minor=True)
    axis.set_yticks(np.arange(-0.5, dimension, 1), minor=True)
    axis.grid(which="minor", color="white", linewidth=0.5)
    axis.tick_params(which="minor", bottom=False, left=False)
    for row_index, column_index in (
        (best["feat_i"], best["feat_j"]),
        (best["feat_j"], best["feat_i"]),
    ):
        axis.add_patch(
            Rectangle(
                (column_index - 0.5, row_index - 0.5),
                1,
                1,
                fill=False,
                edgecolor="red",
                linewidth=2.0,
            )
        )
    colorbar = figure.colorbar(image, ax=axis, fraction=0.035, pad=0.02)
    scope = "70% training partition" if protocol == "strict" else "complete dataset (paper protocol)"
    colorbar.set_label(f"Mean 5-fold CV accuracy - {scope}")
    axis.set_title(
        f"{dataset_name}: best = {best['name_i']} + {best['name_j']} "
        f"({best['cv_mean_acc']:.4f})",
        loc="left",
        fontsize=10,
        fontweight="bold",
    )
    figure.tight_layout()
    output_path = os.path.join(figures_dir, f"{dataset_name}_feature_heatmap.png")
    figure.savefig(output_path, bbox_inches="tight")
    plt.close(figure)
    print("  wrote", output_path)


def run_dataset(dataset_name, dataset, client, protocol, results_dir, figures_dir):
    X_train, X_test, y_train, y_test = _split(dataset.data, dataset.target)
    if protocol == "strict":
        selection_X, selection_y = X_train, y_train
    else:
        selection_X, selection_y = dataset.data, dataset.target

    print(f"\n=== {dataset_name} [{protocol}] ===")
    print(
        f"  fixed 70/30 split: train={len(y_train)}, test={len(y_test)}; "
        f"feature-selection rows={len(selection_y)}"
    )
    rows = cv_accuracy_all_pairs(
        selection_X,
        selection_y,
        list(dataset.feature_names),
        client,
    )
    best = rows[0]
    print(
        f"  selected: {best['name_i']} + {best['name_j']}; "
        f"CV={best['cv_mean_acc']:.6f} +/- {best['cv_std_acc']:.6f}"
    )

    pairs_path = os.path.join(results_dir, f"{dataset_name}_cv_all_pairs.csv")
    with open(pairs_path, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    make_heatmap(dataset_name, rows, best, figures_dir, protocol)

    if protocol == "paper":
        expected_i, expected_j, expected_accuracy = PAPER_FEATURE_SELECTION[dataset_name]
        if (best["feat_i"], best["feat_j"]) != (expected_i, expected_j):
            raise AssertionError(
                f"{dataset_name}: selected pair differs from the current paper: "
                f"got {(best['feat_i'], best['feat_j'])}, expected {(expected_i, expected_j)}"
            )
        if not np.isclose(best["cv_mean_acc"], expected_accuracy, atol=5e-7, rtol=0.0):
            raise AssertionError(
                f"{dataset_name}: CV accuracy {best['cv_mean_acc']} differs from paper {expected_accuracy}"
            )
        print("  paper feature pair/CV accuracy: MATCH")

    return {
        "dataset": dataset_name,
        "protocol": protocol,
        "feature_i": best["name_i"],
        "feature_j": best["name_j"],
        "feat_i": best["feat_i"],
        "feat_j": best["feat_j"],
        "cv_mean_acc": best["cv_mean_acc"],
        "cv_std_acc": best["cv_std_acc"],
        "selection_n": len(selection_y),
        "n_train": len(y_train),
        "n_test": len(y_test),
        "random_state": RANDOM_STATE,
    }


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", choices=["strict", "paper"], default="strict")
    parser.add_argument("--backend", choices=["julia", "mock"], default="julia")
    parser.add_argument("--julia-bin", default="julia")
    parser.add_argument(
        "--julia-threads",
        default="2",
        help="Julia thread count; 2 reproduces the evaluation setup stated in the paper",
    )
    parser.add_argument("--verify-hulls", action="store_true")
    parser.add_argument("--results-dir", default=os.path.join(HERE, "results"))
    parser.add_argument("--figures-dir", default=os.path.join(HERE, "figures"))
    return parser.parse_args()


def main():
    args = parse_args()
    results_dir = os.path.abspath(args.results_dir)
    figures_dir = os.path.abspath(args.figures_dir)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    if args.protocol == "paper":
        print(
            "WARNING: paper protocol performs feature-selection CV on the complete dataset. "
            "Use --protocol strict for the train-only workflow."
        )
    client_kwargs = {"verify_with_scipy": args.verify_hulls}
    if args.backend == "julia":
        client_kwargs.update(julia_bin=args.julia_bin, threads=args.julia_threads)
    client = make_client(args.backend, **client_kwargs)

    summary = []
    for dataset_name, loader in DATASETS.items():
        summary.append(
            run_dataset(
                dataset_name,
                loader(),
                client,
                args.protocol,
                results_dir,
                figures_dir,
            )
        )

    summary_path = os.path.join(results_dir, "feature_selection_summary.csv")
    with open(summary_path, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    print("\nSummary:")
    for row in summary:
        print(
            f"  {row['dataset']:15s} {row['feature_i']} + {row['feature_j']} "
            f"CV={row['cv_mean_acc']:.4f}"
        )
    print("wrote", summary_path)


if __name__ == "__main__":
    main()
