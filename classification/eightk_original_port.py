"""Pure-Python reference port of the bundled Julia implementations.

The algorithm uses 2k regular directions to find 4k extreme points (the
"8k-sided polygon" construction), partitions the remaining points into 4k
candidate regions with ``find_P``, and applies one of four direction-specific
reduction and monotone-chain rules to each region.

Variable names and control flow intentionally mirror the Julia source, with
indices translated from one-based to zero-based notation. This module is used
only for smoke tests and correctness checks; its timings are not reportable.
"""
import math
import numpy as np


def orient(p1, p2, p3):
    dx1 = p2[0] - p1[0]
    dy1 = p2[1] - p1[1]
    dx2 = p3[0] - p1[0]
    dy2 = p3[1] - p1[1]
    return dx1 * dy2 - dy1 * dx2


def starting_vertices(X, k):
    """Returns Is, a dict keyed 1..4k+1 (matching Julia's 1-based Is[j]),
    each value a pair [p_min_or_max, p_max_or_min] (2-element list of
    length-2 numpy arrays), exactly mirroring starting_vertices() /
    starting_vertices_parallel_ver1() in the Julia source (both are
    logically identical -- ver1 just runs the j-loop with @threads)."""
    n = X.shape[0]
    Is = {}
    for j in range(1, 2 * k + 1):
        d = math.pi * (j - 1) / (2 * k)
        e_x = math.cos(d)
        e_y = math.sin(d)
        max_j = e_x * X[0, 0] + e_y * X[0, 1]
        min_j = max_j

        p_max = np.array([X[0, 0], X[0, 1]], dtype=float)
        p_min = np.array([X[0, 0], X[0, 1]], dtype=float)
        p_max2 = np.array([X[0, 0], X[0, 1]], dtype=float)
        p_min2 = np.array([X[0, 0], X[0, 1]], dtype=float)

        if 1 <= j <= k:
            for i in range(1, n):
                x = X[i, 0]
                y = X[i, 1]
                t = e_x * x + e_y * y

                max_j_eq_t = (max_j == t)
                max_j_lt_t = (max_j < t)
                y_gt_y_max = y > p_max[1]
                y_lt_y_min = y < p_min[1]

                if max_j_lt_t:
                    max_j = t
                    p_max[0] = x
                    p_min[0] = x
                    p_max[1] = y
                    p_min[1] = y

                if max_j_eq_t and y_lt_y_min:
                    p_min[0] = x
                    p_min[1] = y
                if max_j_eq_t and y_gt_y_max:
                    p_max[0] = x
                    p_max[1] = y

                min_j_eq_t = (min_j == t)
                min_j_gt_t = (min_j > t)
                y_gt_y_max2 = y > p_max2[1]
                y_lt_y_min2 = y < p_min2[1]

                if min_j_gt_t:
                    min_j = t
                    p_min2[0] = x
                    p_max2[0] = x
                    p_min2[1] = y
                    p_max2[1] = y

                if min_j_eq_t and y_lt_y_min2:
                    p_min2[0] = x
                    p_min2[1] = y
                if min_j_eq_t and y_gt_y_max2:
                    p_max2[0] = x
                    p_max2[1] = y

            Is[j] = [p_min, p_max]
            Is[j + 2 * k] = [p_max2, p_min2]
        else:
            for i in range(1, n):
                x = X[i, 0]
                y = X[i, 1]
                t = e_x * x + e_y * y

                max_j_eq_t = (max_j == t)
                max_j_lt_t = (max_j < t)
                x_gt_x_max = x > p_max[0]
                x_lt_x_min = x < p_min[0]

                if max_j_lt_t:
                    max_j = t
                    p_max[0] = x
                    p_min[0] = x
                    p_max[1] = y
                    p_min[1] = y

                if max_j_eq_t and x_lt_x_min:
                    p_min[0] = x
                    p_min[1] = y
                if max_j_eq_t and x_gt_x_max:
                    p_max[0] = x
                    p_max[1] = y

                min_j_eq_t = (min_j == t)
                min_j_gt_t = (min_j > t)
                x_gt_x_max2 = x > p_max2[0]
                x_lt_x_min2 = x < p_min2[0]

                if min_j_gt_t:
                    min_j = t
                    p_min2[0] = x
                    p_max2[0] = x
                    p_min2[1] = y
                    p_max2[1] = y

                if min_j_eq_t and x_lt_x_min2:
                    p_min2[0] = x
                    p_min2[1] = y
                if min_j_eq_t and x_gt_x_max2:
                    p_max2[0] = x
                    p_max2[1] = y

            Is[j] = [p_max, p_min]
            Is[j + 2 * k] = [p_min2, p_max2]

    Is[4 * k + 1] = [Is[1][0]]
    return Is


def find_P(X, Is, k):
    Pp = {}
    for j in range(1, 4 * k + 1):
        a = Is[j][1]
        b = Is[j + 1][0]
        alpha = b[1] - a[1]
        beta = b[0] - a[0]
        gamma = alpha * a[0] - beta * a[1]
        mask = X[:, 0] * alpha - X[:, 1] * beta > gamma
        extra = np.array([a, b])
        Pp[j] = np.vstack([X[mask], extra])
    return Pp


def reduce_1(P):
    idx = np.argsort(-P[:, 0], kind="stable")
    P = P[idx]
    Pp = [P[0].copy()]
    x_min = P[0, 0]
    y_max = P[0, 1]
    for i in range(1, len(P)):
        x, y = P[i, 0], P[i, 1]
        if y > y_max:
            if x != x_min:
                Pp.append(P[i].copy())
                x_min, y_max = x, y
            else:
                Pp[-1] = P[i].copy()
                y_max = y
    return Pp


def reduce_2(P):
    idx = np.argsort(P[:, 0], kind="stable")
    P = P[idx]
    Pp = [P[0].copy()]
    x_max = P[0, 0]
    y_max = P[0, 1]
    for i in range(1, len(P)):
        x, y = P[i, 0], P[i, 1]
        if y > y_max:
            if x != x_max:
                Pp.append(P[i].copy())
                x_max, y_max = x, y
            else:
                Pp[-1] = P[i].copy()
                y_max = y
    return Pp


def reduce_3(P):
    idx = np.argsort(P[:, 0], kind="stable")
    P = P[idx]
    Pp = [P[0].copy()]
    x_max = P[0, 0]
    y_min = P[0, 1]
    for i in range(1, len(P)):
        x, y = P[i, 0], P[i, 1]
        if y < y_min:
            if x != x_max:
                Pp.append(P[i].copy())
                x_max, y_min = x, y
            else:
                Pp[-1] = P[i].copy()
                y_min = y
    return Pp


def reduce_4(P):
    idx = np.argsort(-P[:, 0], kind="stable")
    P = P[idx]
    Pp = [P[0].copy()]
    x_min = P[0, 0]
    y_min = P[0, 1]
    for i in range(1, len(P)):
        x, y = P[i, 0], P[i, 1]
        if y < y_min:
            if x != x_min:
                Pp.append(P[i].copy())
                x_min, y_min = x, y
            else:
                Pp[-1] = P[i].copy()
                y_min = y
    return Pp


def _monotone_scan(P, pop_le):
    """algorithm_1/3 use `orient(...) <= 0` as the pop condition;
    algorithm_2/4 use `orient(...) >= 0`. pop_le=True selects the
    former, False the latter -- exact translation of the Julia stack
    loop (push/pop via a preallocated, resizing array in Julia; a
    plain Python list here, semantically identical)."""
    l = len(P)
    if l < 3:
        return list(P)
    V = [P[0]]
    for i in range(1, l):
        while len(V) >= 2:
            o = orient(V[-2], V[-1], P[i])
            if (pop_le and o <= 0) or ((not pop_le) and o >= 0):
                V.pop()
            else:
                break
        V.append(P[i])
    return V


def algorithm_1(P):
    return _monotone_scan(P, pop_le=True)


def algorithm_2(P):
    return _monotone_scan(P, pop_le=False)


def algorithm_3(P):
    return _monotone_scan(P, pop_le=True)


def algorithm_4(P):
    return _monotone_scan(P, pop_le=False)


def alg_8ksided4ch(X, k):
    """Sequential port of alg_8ksided4CH(X, k) from 8ksidedPolygon4CH.jl."""
    X = np.asarray(X, dtype=float)
    if X.shape[0] < 3:
        return [X[i] for i in range(X.shape[0])]

    Is = starting_vertices(X, k)
    P = find_P(X, Is, k)
    V = {}
    for i in range(1, 4 * k + 1):
        if 1 <= i <= k:
            Pp = reduce_1(P[i])
            V[i] = algorithm_1(Pp)
        elif k < i <= 2 * k:
            Pp = reduce_2(P[i])
            V[i] = list(reversed(algorithm_2(Pp)))
        elif 2 * k < i <= 3 * k:
            Pp = reduce_3(P[i])
            V[i] = algorithm_3(Pp)
        else:
            Pp = reduce_4(P[i])
            V[i] = list(reversed(algorithm_4(Pp)))

    hull = list(V[1])
    for i in range(2, 4 * k + 1):
        if np.array_equal(hull[-1], V[i][0]):
            hull.extend(V[i][1:])
        else:
            hull.extend(V[i])
    if len(hull) > 1 and np.array_equal(hull[0], hull[-1]):
        hull.pop()
    return hull


def alg_8ksided4ch_parallel_ver1(X, k):
    """The 'ver1' parallel Julia file is logically identical to the
    sequential one (same starting_vertices/find_P math, same reduce_i/
    algorithm_i per-region rules) -- only the two outer loops
    (over 1:2k in starting_vertices, and 1:4k in the main function) are
    wrapped in @threads. Since each iteration writes to a distinct
    array slot and reads only shared input X, the result is
    mathematically identical to the sequential version; this Python
    port therefore reuses the same functions and exists to make that
    equivalence explicit/testable rather than to reimplement anything
    different."""
    return alg_8ksided4ch(X, k)
