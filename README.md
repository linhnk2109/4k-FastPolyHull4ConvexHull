# 4k-FastPolyHull4ConvexHull
4k-FastPolyHull and P-4k-FastPolyHull for finding convex hull of a finite set of points in plane

Below are some of the codes we created to support the paper “Efficient and Parallel Algorithms for Determining Convex Hulls of Large-Scale 2D Data.”

1. 8ksidedPolygon4CH.jl is the code for the \textsc{$4k$-FastPolyHull} algorithm, which uses the 8k-sided polygon as a preprocessing step.
2. 8ksidedPolygon4CH_parallel_ver1.jl is a parallel version of the \textsc{$4k$-FastPolyHull} algorithm, used in the presentation of the parallelization results of the \textsc{$4k$-FastPolyHull} algorithm in the paper “Efficient and Parallel Algorithms for Determining Convex Hulls of Large-Scale 2D Data.”
3. 8ksidedPolygon4CH_parallel_ver2.jl is another parallel version of the \textsc{$4k$-FastPolyHull} algorithm.
4. quickhull.jl is the code to implement the Quickhull algorithm from the Qhull library ( http://www.qhull.org/).

### Install libraries
- Open a terminal and go to the directory containing the code
>        julia install.jl

**Run the programs**
- Open a terminal and go to the directory containing the code.
>        julia main.jl

- Run Algorithms in parallel mode.
>        julia -t numberOfThreads main.jl

**Note**
Creat a file "result" in the directory containing the codes.

**Setting**
Benchmarking mode
Set benchmarking = true in the main functions.

Export the convex hull to file
Set benchmarking = false and exportResult = true in the main functions.
   
