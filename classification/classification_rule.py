"""Section 5.1 - the three-case convex-hull classification rule.

This module deliberately does not compute convex hulls.  It accepts the
vertices returned by Julia and applies exactly the rule stated in the paper:

1. inside exactly one class hull -> that class;
2. inside several hulls -> nearest class centroid among those hulls;
3. outside every hull -> nearest hull boundary.

Membership is inclusive of the boundary, as required by ``x in conv(X_c)``.
The implementation is vectorized over query points so it remains practical for
the million-point synthetic experiment.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np


DEFAULT_ATOL = 1.0e-12
DEFAULT_RTOL = 1.0e-12
DEFAULT_CHUNK_SIZE = 250_000


def _as_points(values, *, name: str) -> np.ndarray:
    points = np.asarray(values, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError(f"{name} must have shape (n, 2); got {points.shape}")
    if not np.all(np.isfinite(points)):
        raise ValueError(f"{name} contains NaN or infinity")
    return points


def _clean_hull(vertices) -> np.ndarray:
    """Remove repeated vertices without changing the cyclic order."""
    hull = _as_points(vertices, name="hull")
    if len(hull) == 0:
        return hull.copy()

    kept = [hull[0]]
    for point in hull[1:]:
        if not np.array_equal(point, kept[-1]):
            kept.append(point)
    hull = np.asarray(kept, dtype=np.float64)
    if len(hull) > 1 and np.array_equal(hull[0], hull[-1]):
        hull = hull[:-1]

    # A valid convex hull should not repeat non-consecutive vertices.  Removing
    # such repetitions here makes the rule robust to a bridge/backend that
    # happens to close or concatenate its chains redundantly.
    unique = []
    seen = set()
    for point in hull:
        key = (float(point[0]), float(point[1]))
        if key not in seen:
            seen.add(key)
            unique.append(point)
    return np.asarray(unique, dtype=np.float64).reshape(-1, 2)


def point_to_segment_distances(points, a, b) -> np.ndarray:
    """Euclidean distance from every row of *points* to segment ``[a,b]``."""
    points = _as_points(points, name="points")
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ab = b - a
    denom = float(np.dot(ab, ab))
    if denom == 0.0:
        return np.linalg.norm(points - a, axis=1)
    t = np.clip(((points - a) @ ab) / denom, 0.0, 1.0)
    projection = a + t[:, None] * ab
    return np.linalg.norm(points - projection, axis=1)


def point_to_segment_dist(point, a, b) -> float:
    point = np.asarray(point, dtype=np.float64).reshape(1, 2)
    return float(point_to_segment_distances(point, a, b)[0])


def point_to_polygon_distances(
    points, polygon, *, chunk_size: int = DEFAULT_CHUNK_SIZE
) -> np.ndarray:
    """Distance to a polygon boundary for all points.

    Hulls with one or two vertices are treated as a point or segment.  Work is
    chunked to keep memory bounded for multi-million-point test sets.
    """
    points = _as_points(points, name="points")
    polygon = _clean_hull(polygon)
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    if len(polygon) == 0:
        return np.full(len(points), np.inf, dtype=np.float64)
    if len(polygon) == 1:
        return np.linalg.norm(points - polygon[0], axis=1)

    distances = np.empty(len(points), dtype=np.float64)
    edge_count = 1 if len(polygon) == 2 else len(polygon)
    for start in range(0, len(points), chunk_size):
        stop = min(start + chunk_size, len(points))
        chunk = points[start:stop]
        best_sq = np.full(len(chunk), np.inf, dtype=np.float64)
        for edge in range(edge_count):
            a = polygon[edge]
            b = polygon[edge + 1] if len(polygon) == 2 else polygon[(edge + 1) % len(polygon)]
            ab = b - a
            denom = float(np.dot(ab, ab))
            if denom == 0.0:
                delta = chunk - a
            else:
                t = np.clip(((chunk - a) @ ab) / denom, 0.0, 1.0)
                delta = chunk - (a + t[:, None] * ab)
            best_sq = np.minimum(best_sq, np.einsum("ij,ij->i", delta, delta))
        distances[start:stop] = np.sqrt(best_sq)
    return distances


def point_to_polygon_dist(point, polygon) -> float:
    point = np.asarray(point, dtype=np.float64).reshape(1, 2)
    return float(point_to_polygon_distances(point, polygon)[0])


def convex_polygon_contains_points(
    points,
    polygon,
    *,
    atol: float = DEFAULT_ATOL,
    rtol: float = DEFAULT_RTOL,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> np.ndarray:
    """Inclusive point-in-convex-polygon test.

    Unlike ``matplotlib.path.Path.contains_points``, this function defines
    boundary membership explicitly.  Both clockwise and counter-clockwise hull
    order are accepted.
    """
    points = _as_points(points, name="points")
    polygon = _clean_hull(polygon)
    if atol < 0.0 or rtol < 0.0:
        raise ValueError("atol and rtol must be non-negative")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    if len(polygon) == 0:
        return np.zeros(len(points), dtype=bool)
    if len(polygon) <= 2:
        scale = max(1.0, float(np.max(np.abs(polygon))))
        tolerance = atol + rtol * scale
        return point_to_polygon_distances(points, polygon, chunk_size=chunk_size) <= tolerance

    x, y = polygon[:, 0], polygon[:, 1]
    signed_twice_area = float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    polygon_scale = max(1.0, float(np.max(np.abs(polygon))))
    if abs(signed_twice_area) <= atol + rtol * polygon_scale * polygon_scale:
        # Degenerate collinear "polygon": membership means lying on its
        # boundary segments.
        tolerance = atol + rtol * polygon_scale
        return point_to_polygon_distances(points, polygon, chunk_size=chunk_size) <= tolerance

    ccw = signed_twice_area > 0.0
    result = np.empty(len(points), dtype=bool)
    for start in range(0, len(points), chunk_size):
        stop = min(start + chunk_size, len(points))
        chunk = points[start:stop]
        inside = np.ones(len(chunk), dtype=bool)
        for edge in range(len(polygon)):
            a = polygon[edge]
            b = polygon[(edge + 1) % len(polygon)]
            ab = b - a
            rel = chunk - a
            cross = ab[0] * rel[:, 1] - ab[1] * rel[:, 0]
            # The cross product has squared-coordinate units.  Scale the
            # tolerance by edge length and query distance to remain stable
            # after standardization or a change of coordinate magnitude.
            tolerance = atol + rtol * max(1.0, float(np.linalg.norm(ab))) * np.maximum(
                1.0, np.linalg.norm(rel, axis=1)
            )
            inside &= cross >= -tolerance if ccw else cross <= tolerance
            if not np.any(inside):
                break
        result[start:stop] = inside
    return result


class HullClassificationRule:
    """Classifier implementing the three cases in Section 5.1."""

    def __init__(
        self,
        *,
        atol: float = DEFAULT_ATOL,
        rtol: float = DEFAULT_RTOL,
        chunk_size: int = DEFAULT_CHUNK_SIZE,
    ):
        self.atol = float(atol)
        self.rtol = float(rtol)
        self.chunk_size = int(chunk_size)

    def fit_from_hulls(self, hulls: Mapping, centroids: Mapping):
        if not hulls:
            raise ValueError("hulls must contain at least one class")
        if set(hulls) != set(centroids):
            raise ValueError("hulls and centroids must contain the same class labels")
        try:
            classes = sorted(hulls)
        except TypeError:
            classes = list(hulls)

        self.classes_ = classes
        self.class_array_ = np.asarray(classes)
        self.hulls_ = {label: _clean_hull(hulls[label]) for label in classes}
        self.centroids_ = {}
        for label in classes:
            centroid = np.asarray(centroids[label], dtype=np.float64)
            if centroid.shape != (2,) or not np.all(np.isfinite(centroid)):
                raise ValueError(f"centroid for class {label!r} must be a finite (2,) vector")
            if len(self.hulls_[label]) == 0:
                raise ValueError(f"hull for class {label!r} is empty")
            self.centroids_[label] = centroid
        return self

    def _require_fitted(self):
        if not hasattr(self, "classes_"):
            raise RuntimeError("call fit_from_hulls before predict")

    def predict_one(self, point):
        point = np.asarray(point, dtype=np.float64)
        if point.shape != (2,):
            raise ValueError(f"point must have shape (2,); got {point.shape}")
        return self.predict(point.reshape(1, 2))[0]

    def predict(self, X):
        self._require_fitted()
        X = _as_points(X, name="X")
        n_points = len(X)
        n_classes = len(self.classes_)
        if n_points == 0:
            return np.empty(0, dtype=self.class_array_.dtype)

        inside = np.empty((n_points, n_classes), dtype=bool)
        for class_index, label in enumerate(self.classes_):
            inside[:, class_index] = convex_polygon_contains_points(
                X,
                self.hulls_[label],
                atol=self.atol,
                rtol=self.rtol,
                chunk_size=self.chunk_size,
            )

        predictions = np.empty(n_points, dtype=self.class_array_.dtype)
        number_containing = inside.sum(axis=1)

        single_mask = number_containing == 1
        if np.any(single_mask):
            class_indices = np.argmax(inside[single_mask], axis=1)
            predictions[single_mask] = self.class_array_[class_indices]

        multiple_mask = number_containing > 1
        if np.any(multiple_mask):
            points = X[multiple_mask]
            centroids = np.asarray([self.centroids_[label] for label in self.classes_])
            squared_distances = np.sum((points[:, None, :] - centroids[None, :, :]) ** 2, axis=2)
            squared_distances[~inside[multiple_mask]] = np.inf
            predictions[multiple_mask] = self.class_array_[np.argmin(squared_distances, axis=1)]

        outside_mask = number_containing == 0
        if np.any(outside_mask):
            points = X[outside_mask]
            best_distance = np.full(len(points), np.inf, dtype=np.float64)
            best_class_index = np.zeros(len(points), dtype=np.int64)
            for class_index, label in enumerate(self.classes_):
                distance = point_to_polygon_distances(
                    points, self.hulls_[label], chunk_size=self.chunk_size
                )
                # Strict '<' intentionally keeps the first class in class order
                # when distances tie, making the otherwise unspecified argmin
                # deterministic.
                better = distance < best_distance
                best_distance[better] = distance[better]
                best_class_index[better] = class_index
            predictions[outside_mask] = self.class_array_[best_class_index]

        return predictions


def _self_test() -> None:
    hulls = {
        "A": np.array([[0, 0], [2, 0], [2, 2], [0, 2]], dtype=float),
        "B": np.array([[1.5, 1.5], [3.5, 1.5], [3.5, 3.5], [1.5, 3.5]], dtype=float),
    }
    centroids = {"A": np.array([1.0, 1.0]), "B": np.array([2.5, 2.5])}
    classifier = HullClassificationRule().fit_from_hulls(hulls, centroids)
    cases = [
        ([0.5, 0.5], "A", "inside exactly A"),
        ([3.0, 3.0], "B", "inside exactly B"),
        ([1.8, 1.8], "B", "overlap; nearest centroid"),
        ([5.0, 5.0], "B", "outside; nearest boundary"),
        ([-5.0, -5.0], "A", "outside; nearest boundary"),
        ([0.0, 1.0], "A", "boundary is inside"),
    ]
    points = np.asarray([case[0] for case in cases], dtype=float)
    vector_predictions = classifier.predict(points)
    for index, (point, expected, description) in enumerate(cases):
        scalar = classifier.predict_one(point)
        vector = vector_predictions[index]
        if scalar != expected or vector != expected:
            raise AssertionError(
                f"{description}: expected {expected!r}, got scalar={scalar!r}, vector={vector!r}"
            )

    # Degenerate hulls are useful for small smoke tests even though all paper
    # classes have many training samples.
    degenerate = HullClassificationRule().fit_from_hulls(
        {0: np.array([[0.0, 0.0]]), 1: np.array([[2.0, 0.0], [4.0, 0.0]])},
        {0: np.array([0.0, 0.0]), 1: np.array([3.0, 0.0])},
    )
    if not np.array_equal(degenerate.predict([[0, 0], [3, 0], [5, 0]]), [0, 1, 1]):
        raise AssertionError("degenerate point/segment hull handling failed")
    print(f"classification_rule.py: {len(cases) + 3} checks passed")


if __name__ == "__main__":
    _self_test()
