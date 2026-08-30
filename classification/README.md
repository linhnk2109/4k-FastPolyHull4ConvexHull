# Classification Experiments (Section 5)

This directory contains the Python and Julia code used for the classification
application in Section 5 of the 4k-FastPolyHull paper. Python handles data
generation, splitting, standardization, feature selection, classification,
evaluation, and plotting. Julia is used only to compute convex hull vertices
with 4k-FastPolyHull or P-4k-FastPolyHull.

The default workflow implements the leakage-free protocol reported in the
revised manuscript:

- stratified 70/30 train/test splitting before feature selection;
- stratified five-fold cross-validation on the training partition only;
- a separate `StandardScaler` fitted inside every cross-validation fold;
- six synthetic sizes from 1,000 through 100,000,000 samples;
- three timing repetitions summarized by the median; and
- four Julia threads for P-4k-FastPolyHull in the synthetic experiment.

## Directory contents

| File | Section | Purpose |
|---|---:|---|
| `classification_rule.py` | 5.1 | Implements the three-case convex-hull classification rule. |
| `experiment_synthetic.py` | 5.2 | Runs the Gaussian-blob benchmark and exports tables and figures. |
| `feature_selection.py` | 5.3 | Selects a feature pair using train-only five-fold cross-validation. |
| `experiment_real.py` | 5.4 | Evaluates the selected pair on the held-out test set. |
| `hull_client.py` | 5.1-5.4 | Implements the Python side of the binary Python/Julia bridge. |
| `hull_bridge.jl` | 5.1-5.4 | Calls the sequential or threaded Julia algorithm and returns hull vertices. |
| `refresh_paper_figures.py` | 5.2 | Redraws the synthetic figures from saved CSV results. |
| `8ksidedPolygon4CH.jl` | bridge | Library-safe copy of 4k-FastPolyHull. |
| `8ksidedPolygon4CH_parallel_ver1.jl` | bridge | Library-safe copy of P-4k-FastPolyHull. |
| `utils.jl` | bridge | Minimal orientation helper required by the Julia implementations. |
| `eightk_original_port.py` | testing | Python-only fallback for smoke tests; its timings are not reportable. |
| `requirements.txt` | setup | Tested Python dependency versions. |

The bridge transports coordinates as little-endian Float64 values and uses a
small TSV manifest for job metadata. This avoids CSV precision and file-size
problems at the largest sample sizes. The Julia bridge uses only the Julia
standard library.

## Requirements

- Python 3.12 or later
- Julia 1.x with threading support
- the Python packages listed in `requirements.txt`

The reported rerun used Python 3.12.5, scikit-learn 1.7.0, and Julia 1.7.3.

From the repository root:

```bash
cd classification
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
julia --version
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

If Julia is not on `PATH`, pass its absolute executable path to an experiment:

```bash
--julia-bin /absolute/path/to/julia
```

## Quick validation

Test the classification rule:

```bash
python classification_rule.py
```

Expected final line:

```text
classification_rule.py: 9 checks passed
```

Test both Julia code paths with four threads and compare their convex regions
with SciPy/Qhull:

```bash
python hull_client.py --backend julia --threads 4
```

Expected final line:

```text
hull_client.py: bridge smoke test passed
```

If Julia is temporarily unavailable, the bridge wiring can be checked with the
Python reference port:

```bash
python hull_client.py --backend mock --threads 2
```

Mock-backend timings are for diagnostics only and must not be reported.

## Sections 5.3 and 5.4: real datasets

Run feature selection before the final evaluation because
`experiment_real.py` reads `feature_selection_summary.csv`:

```bash
python feature_selection.py \
  --protocol strict \
  --backend julia \
  --julia-threads 2 \
  --verify-hulls \
  --results-dir results/strict \
  --figures-dir figures/strict

python experiment_real.py \
  --protocol strict \
  --backend julia \
  --julia-threads 2 \
  --verify-hulls \
  --results-dir results/strict \
  --figures-dir figures/strict
```

The `strict` protocol performs the split first. Every candidate feature pair
is then evaluated by stratified five-fold cross-validation on the training
partition. Each fold fits its own scaler. The held-out test partition is used
only after the winning feature pair has been selected.

Expected deterministic results with the pinned dependencies are:

| Dataset | Training-selected pair | Mean CV | Convex hull | k-NN | RBF SVM |
|---|---|---:|---:|---:|---:|
| Iris | sepal length + petal length | 97.14% | 91.11% | 91.11% | 95.56% |
| Wine | alcohol + flavanoids | 94.37% | 88.89% | 88.89% | 88.89% |
| Breast Cancer Wisconsin | worst area + worst concave points | 95.98% | 94.74% | 92.98% | 93.57% |

## Section 5.2: synthetic experiment

Run a small validation before starting the large benchmark:

```bash
python experiment_synthetic.py \
  --profile paper \
  --sizes 1000 10000 100000 \
  --backend julia \
  --julia-threads 4 \
  --verify-hulls \
  --results-dir results/synthetic_smoke \
  --figures-dir figures/synthetic_smoke
```

Run all six manuscript sizes:

```bash
python experiment_synthetic.py \
  --profile paper \
  --backend julia \
  --julia-threads 4 \
  --results-dir results/synthetic_paper \
  --figures-dir figures/synthetic_paper
```

The nominal sizes are 1,000, 10,000, 100,000, 1,000,000, 10,000,000, and
100,000,000. To keep all three classes equal, the paper profile generates
`3 * floor(n/3)` observations. The CSV output stores both the nominal and
actual totals.

RBF SVM is evaluated only up to the default limit of 20,000 samples. Each
method is timed three times after Julia warm-up, and the median fit-plus-predict
time is reported. Do not enable `--verify-hulls` for the 100-million-point run
unless sufficient additional memory is available, because verification also
invokes SciPy/Qhull.

The 100-million-point experiment is resource intensive. Run it on a machine
with adequate RAM, free disk space for the temporary binary bridge files, and
no competing heavy workloads. Wall-clock timings are machine dependent;
accuracies are deterministic under the pinned software versions.

## Redrawing the synthetic figures

After `results/synthetic_paper/synthetic_benchmark_results.csv` is available,
redraw all synthetic panels without repeating the largest runs:

```bash
python refresh_paper_figures.py
```

The script reads all six saved rows and recomputes only the deterministic
1,000-point illustration used for the data and decision-region panels.

## Outputs

`feature_selection.py` writes:

- `{dataset}_cv_all_pairs.csv`;
- `feature_selection_summary.csv`; and
- `{dataset}_feature_heatmap.png`.

`experiment_real.py` writes:

- `real_data_results.csv`; and
- `real_data_hulls.png`.

`experiment_synthetic.py` writes:

- `synthetic_benchmark_results.csv`;
- `synthetic_accuracy_table.csv`;
- `synthetic_runtime_table.csv`;
- `synthetic_paper_comparison.csv`;
- `synthetic_overview.png`;
- `synthetic_accuracy.png`; and
- `synthetic_runtime.png`.

Canonical CSV results are stored in `results/strict` and
`results/synthetic_paper`. Generated figures are stored in the corresponding
subdirectories under `figures`.

## Useful options

- `--backend mock` checks the Python workflow without Julia; its hull timings
  are not reportable.
- `--verify-hulls` compares returned convex regions with SciPy/Qhull.
- `--repeats N` changes the number of synthetic timing repetitions.
- `--rbf-max-n N` changes the largest size evaluated by RBF SVM.
- `--sizes ...` overrides the selected synthetic size profile.
- `--julia-threads N` sets the Julia thread count for the parallel algorithm.
- `--no-plots` skips plotting during a trial timing run.

The optional `--protocol paper` mode in the real-data scripts is retained only
to audit an earlier full-dataset feature-selection workflow. The revised
manuscript and the default commands use the leakage-free `strict` protocol.

If Matplotlib cannot write to its default cache directory, select a writable
cache location before running an experiment:

```bash
export MPLCONFIGDIR=/tmp/matplotlib-cache
mkdir -p "$MPLCONFIGDIR"
```
