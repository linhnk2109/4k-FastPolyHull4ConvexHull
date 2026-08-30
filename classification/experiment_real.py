"""Section 5.4 - held-out classification on the three real datasets."""

from __future__ import annotations

import argparse
import csv
import os

import numpy as np
from sklearn.datasets import load_breast_cancer, load_iris, load_wine
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from classification_rule import HullClassificationRule
from hull_client import make_client


HERE = os.path.dirname(os.path.abspath(__file__))
RANDOM_STATE = 42
K_HULL = 1

DATASETS = {
    "iris": load_iris,
    "wine": load_wine,
    "breast_cancer": load_breast_cancer,
}

DISPLAY_NAMES = {
    "iris": "Iris",
    "wine": "Wine",
    "breast_cancer": "Breast Cancer Wisconsin",
}

PAPER_TEST_ACCURACY = {
    "iris": (0.9555555555555556, 0.9111111111111111, 0.9111111111111111),
    "wine": (0.8888888888888888, 0.8888888888888888, 0.8888888888888888),
    "breast_cancer": (0.9473684210526315, 0.9473684210526315, 0.9415204678362573),
}


def load_selected_pairs(results_dir, protocol):
    path = os.path.join(results_dir, "feature_selection_summary.csv")
    if not os.path.exists(path):
        raise SystemExit(
            f"Missing {path}. Run feature_selection.py --protocol {protocol} "
            f"--results-dir {shlex_path(results_dir)} first."
        )
    with open(path, encoding="utf-8", newline="") as source:
        selected = {row["dataset"]: row for row in csv.DictReader(source)}
    missing = set(DATASETS) - set(selected)
    if missing:
        raise SystemExit(f"Feature-selection summary is missing datasets: {sorted(missing)}")
    for dataset_name, row in selected.items():
        if row.get("protocol") != protocol:
            raise SystemExit(
                f"Summary protocol mismatch for {dataset_name}: file has {row.get('protocol')!r}, "
                f"command requested {protocol!r}"
            )
        if int(row.get("random_state", -1)) != RANDOM_STATE:
            raise SystemExit(f"Summary random_state mismatch for {dataset_name}")
    return selected


def shlex_path(path):
    import shlex

    return shlex.quote(path)


def run_dataset(dataset_name, dataset, selected, client, protocol):
    feature_i = int(selected["feat_i"])
    feature_j = int(selected["feat_j"])
    feature_names = (selected["feature_i"], selected["feature_j"])

    X_train_full, X_test_full, y_train, y_test = train_test_split(
        dataset.data,
        dataset.target,
        test_size=0.30,
        stratify=dataset.target,
        random_state=RANDOM_STATE,
    )
    raw_train = X_train_full[:, [feature_i, feature_j]]
    raw_test = X_test_full[:, [feature_i, feature_j]]
    scaler = StandardScaler().fit(raw_train)
    X_train = scaler.transform(raw_train)
    X_test = scaler.transform(raw_test)

    classes = np.unique(y_train).tolist()
    jobs = []
    class_points = {}
    for class_index, label in enumerate(classes):
        class_points[label] = X_train[y_train == label]
        jobs.append(
            {
                "id": f"{dataset_name}__class-{class_index}",
                "group": dataset_name,
                "points": class_points[label],
                "k": K_HULL,
                "mode": "seq",
                "time_it": False,
            }
        )
    backend_results = client.run_batch(jobs)
    hulls = {
        label: backend_results[f"{dataset_name}__class-{class_index}"]["hull"]
        for class_index, label in enumerate(classes)
    }
    centroids = {label: class_points[label].mean(axis=0) for label in classes}
    hull_classifier = HullClassificationRule().fit_from_hulls(hulls, centroids)
    hull_predictions = hull_classifier.predict(X_test)
    hull_accuracy = float(np.mean(hull_predictions == y_test))

    knn = KNeighborsClassifier(n_neighbors=5).fit(X_train, y_train)
    knn_accuracy = float(knn.score(X_test, y_test))
    rbf_svm = SVC(kernel="rbf").fit(X_train, y_train)
    rbf_svm_accuracy = float(rbf_svm.score(X_test, y_test))

    print(
        f"\n{DISPLAY_NAMES[dataset_name]} [{protocol}] - {feature_names[0]} + {feature_names[1]}\n"
        f"  train/test={len(y_train)}/{len(y_test)}; "
        f"Convex hull={hull_accuracy:.4f}, k-NN={knn_accuracy:.4f}, "
        f"RBF SVM={rbf_svm_accuracy:.4f}"
    )

    if protocol == "paper":
        actual = (hull_accuracy, knn_accuracy, rbf_svm_accuracy)
        expected = PAPER_TEST_ACCURACY[dataset_name]
        if not np.allclose(actual, expected, atol=5e-13, rtol=0.0):
            raise AssertionError(
                f"{dataset_name}: test accuracies differ from paper; actual={actual}, expected={expected}"
            )
        print("  paper held-out accuracies: MATCH")

    row = {
        "dataset": dataset_name,
        "protocol": protocol,
        "feature_i": feature_names[0],
        "feature_j": feature_names[1],
        "cv_mean_acc": float(selected["cv_mean_acc"]),
        "n_train": len(y_train),
        "n_test": len(y_test),
        "hull_test_acc": hull_accuracy,
        "knn_test_acc": knn_accuracy,
        "rbf_svm_test_acc": rbf_svm_accuracy,
    }
    plot_data = {
        "X_train": X_train,
        "y_train": y_train,
        "X_test": X_test,
        "y_test": y_test,
        "hull_predictions": hull_predictions,
        "hulls": hulls,
        "classes": classes,
        "class_names": list(dataset.target_names),
        "feature_names": feature_names,
        "accuracy": hull_accuracy,
    }
    return row, plot_data


def make_figure(all_plot_data, output_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch, Polygon

    colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
    dataset_names = list(all_plot_data)
    figure, axes = plt.subplots(1, len(dataset_names), figsize=(6 * len(dataset_names), 5.4), dpi=180)
    if len(dataset_names) == 1:
        axes = [axes]

    for axis, dataset_name in zip(axes, dataset_names):
        data = all_plot_data[dataset_name]
        legend_handles = []
        for class_index, label in enumerate(data["classes"]):
            color = colors[class_index % len(colors)]
            train_points = data["X_train"][data["y_train"] == label]
            class_name = data["class_names"][int(label)]
            axis.scatter(
                train_points[:, 0],
                train_points[:, 1],
                s=10,
                color=color,
                alpha=0.55,
                linewidths=0,
                zorder=2,
            )
            hull = data["hulls"][label]
            if len(hull) >= 3:
                axis.add_patch(
                    Polygon(
                        hull,
                        closed=True,
                        facecolor=color,
                        edgecolor=color,
                        alpha=0.16,
                        linewidth=1.8,
                        zorder=1,
                    )
                )
            elif len(hull) == 2:
                axis.plot(hull[:, 0], hull[:, 1], color=color, linewidth=1.8)
            legend_handles.append(Patch(facecolor=color, edgecolor=color, alpha=0.35, label=class_name))

        correct = data["hull_predictions"] == data["y_test"]
        axis.scatter(
            data["X_test"][correct, 0],
            data["X_test"][correct, 1],
            marker="x",
            s=25,
            color="black",
            linewidths=1.0,
            zorder=4,
        )
        axis.scatter(
            data["X_test"][~correct, 0],
            data["X_test"][~correct, 1],
            marker="x",
            s=34,
            color="red",
            linewidths=1.3,
            zorder=5,
        )
        legend_handles.extend(
            [
                Line2D([], [], marker="x", linestyle="None", color="black", label="correct test"),
                Line2D([], [], marker="x", linestyle="None", color="red", label="misclassified test"),
            ]
        )
        first_feature, second_feature = data["feature_names"]
        axis.set_xlabel(f"{first_feature} (standardized)", fontsize=8)
        axis.set_ylabel(f"{second_feature} (standardized)", fontsize=8)
        axis.set_title(
            f"{DISPLAY_NAMES[dataset_name]}\nConvex-hull test accuracy = {data['accuracy']:.3f}",
            loc="left",
            fontsize=10,
            fontweight="bold",
        )
        axis.legend(handles=legend_handles, fontsize=6, frameon=False, ncol=2)
        axis.grid(alpha=0.20)
        axis.spines[["top", "right"]].set_visible(False)

    figure.suptitle("Section 5.4 - convex-hull classification", fontsize=12, fontweight="bold")
    figure.tight_layout(rect=[0, 0, 1, 0.94])
    figure.savefig(output_path, bbox_inches="tight")
    plt.close(figure)
    print("wrote", output_path)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", choices=["strict", "paper"], default="strict")
    parser.add_argument("--backend", choices=["julia", "mock"], default="julia")
    parser.add_argument("--julia-bin", default="julia")
    parser.add_argument("--julia-threads", default="2")
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
    selected = load_selected_pairs(results_dir, args.protocol)

    client_kwargs = {"verify_with_scipy": args.verify_hulls}
    if args.backend == "julia":
        client_kwargs.update(julia_bin=args.julia_bin, threads=args.julia_threads)
    client = make_client(args.backend, **client_kwargs)

    rows = []
    plot_data = {}
    for dataset_name, loader in DATASETS.items():
        row, data = run_dataset(
            dataset_name,
            loader(),
            selected[dataset_name],
            client,
            args.protocol,
        )
        rows.append(row)
        plot_data[dataset_name] = data

    output_csv = os.path.join(results_dir, "real_data_results.csv")
    with open(output_csv, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    make_figure(plot_data, os.path.join(figures_dir, "real_data_hulls.png"))
    print("wrote", output_csv)


if __name__ == "__main__":
    main()
