# 4k-FastPolyHull for Planar Convex Hulls

This repository provides the Julia implementations of 4k-FastPolyHull and its
parallel variants for computing the convex hull of a finite set of points in
the plane. It also includes the Python/Julia experiments for the classification
application presented in Section 5 of the accompanying paper.

## Repository layout

### Convex-hull algorithms

| File | Purpose |
|---|---|
| `8ksidedPolygon4CH.jl` | Sequential 4k-FastPolyHull implementation. |
| `8ksidedPolygon4CH_parallel_ver1.jl` | P-4k-FastPolyHull implementation used for the reported parallel results. |
| `8ksidedPolygon4CH_parallel_ver2.jl` | Alternative parallel implementation. |
| `quickHull.jl` | Wrapper around QHull.jl for comparison. |
| `main.jl` | Original synthetic-data benchmark driver. |
| `compareOutput.jl` | Utility for comparing exported hull results. |
| `utils.jl` | Data generators, geometry helpers, and reporting utilities. |

### Classification application

The [`classification`](classification/) directory reproduces Sections 5.1-5.4:

- the three-case convex-hull classification rule;
- Gaussian-blob experiments from 1,000 to 100,000,000 samples;
- train-only, five-fold feature selection for Iris, Wine, and Breast Cancer
  Wisconsin; and
- held-out comparison with k-NN and RBF SVM.

It includes the Python/Julia bridge, pinned Python dependencies, canonical CSV
results, and scripts that regenerate the figures. See the
[`classification/README.md`](classification/README.md) for complete setup,
validation, and reproduction instructions.

## Quick start: classification experiments

The classification bridge uses library-safe copies of the sequential and
parallel algorithms and requires no third-party Julia package.

```bash
git clone https://github.com/linhnk2109/4k-FastPolyHull4ConvexHull.git
cd 4k-FastPolyHull4ConvexHull/classification

python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python classification_rule.py
python hull_client.py --backend julia --threads 4
```

## Running the original Julia benchmark driver

The root-level driver preserves the original benchmark configuration. Install
the Julia packages imported by `main.jl` and `utils.jl`, create its output
directory, review the sizes and switches near the top of `main.jl`, and then
run:

```bash
mkdir -p result
julia --threads 4 main.jl
```

The default configuration includes sample sizes up to 100 million points and
is therefore resource intensive. Reduce `sizes` in `main.jl` for a validation
run before launching the complete benchmark.

The main switches are:

- `benchmarking = true`: measure execution time;
- `benchmarking = false`: run the export path;
- `exportResult = true`: export hull vertices when benchmarking is disabled;
- `setK`: values of the 4k-FastPolyHull parameter `k`; and
- `dataType`: the synthetic point distribution.

## Reproducibility notes

- Synthetic random generation uses a fixed seed where specified by the
  experiment driver.
- The Section 5 workflow stores the nominal and actual generated sample counts
  separately because equal class sizes use `3 * floor(n/3)` observations.
- P-4k-FastPolyHull must be started with multiple Julia threads to execute the
  threaded path.
- Timing results depend on hardware, Julia/Python versions, and system load.

## Code availability

The implementations and experiment scripts supporting the paper are publicly
available in this repository:

<https://github.com/linhnk2109/4k-FastPolyHull4ConvexHull>

## License

No license file is currently included. Add an appropriate open-source license
before allowing reuse or redistribution beyond the terms stated by the
authors.
