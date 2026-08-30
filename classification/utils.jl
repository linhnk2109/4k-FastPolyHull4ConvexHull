## Minimal utility module shared by the sequential and threaded algorithms.
## Only orient() is required by the classification pipeline. Data generators,
## reporting helpers, and their CSV/DataFrames/StatsBase/JLD2 dependencies are
## intentionally omitted so the Julia bridge uses only the standard library.
@inline function orient(p1, p2, p3)
    dx1 = p2[1] - p1[1]
    dy1 = p2[2] - p1[2]
    dx2 = p3[1] - p1[1]
    dy2 = p3[2] - p1[2]
    return dx1 * dy2 - dy1 * dx2
end
