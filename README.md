# SNAP-HiC-Diff

**SNP-Aware Pipeline for Hi-C Differential Analysis**

A Snakemake-based bioinformatics workflow for detecting differential chromatin interactions in Hi-C sequencing data. Processes multi-resolution `.mcool` contact map files through normalization, sliding-window decomposition, gap detection, and statistical comparison.

---

## Overview

The pipeline takes ICE-balanced Hi-C contact maps (`.mcool`) as input and produces:

- O/E (Observed over Expected) normalized matrices
- Sliding-window decomposition with gap-filling and Gaussian smoothing
- Gap quality control metrics
- Differential analysis with MSE, SSIM, Spearman correlation, NCC, and Tracy-Widom statistics
- Breakpoint detection via NaN density patterns

Supported samples: **HC04_A01-A13** (A compartments) and **HC04_D01-D13** (D compartments), 26 chromosomes total.

---

## Quick Start

### 1. Setup environment

On an ARM login node with Conda channel access, build every workflow environment. Repeat this command after changing any `env/*.yaml` file:

```bash
bash scripts/conda_envs.sh bootstrap
```

Compute nodes are offline. Validate the prepared cache without channel access:

```bash
CONDA_OFFLINE=true bash scripts/conda_envs.sh verify
```

### 2. Configure

Edit `config/config.yaml` to set paths and parameters:

```yaml
# Sample manifest (TSV)
manifest: data/manifest.tsv

# Matrix resolution (bp)
resolution: 4000

# Resolution-derived genomic windows
window_analysis:
  window_bins: [25, 50, 75, 100, 125]
  step_fraction: 0.20

# Smoothing sigma
smoothing:
  sigma: 3
```

For the CWorld/Meltron insulation workflow, CWorld IS is calculated from the balanced `.cool` resolution configured at
`meltron.cworld_is.source_resolution_bp`. Set this value to the source cooler
resolution: an integer of at least 1,000 bp (for example `4000` or `10000`).
Meltron rejects lower resolutions. The converter writes each
bedGraph bin start as Meltron `viewpoint` and verifies that every complete
bedGraph bin has exactly that width, so Meltron infers the configured
resolution rather than a midpoint-derived half-resolution.

Edit `data/manifest.tsv` with your sample paths.

### 3. Run locally (dry-run)

```bash
CONDA_OFFLINE=true snakemake -n --use-conda --conda-prefix "$PWD/.snakemake/conda" --configfile config/config.yaml
```

### 4. Run on HPC

```bash
# Submit only after the login-node bootstrap and offline verification succeed.
dsub -s run_snakemake.sh
```

See [docs/dsub_hpc_standard.md](docs/dsub_hpc_standard.md) for HPC setup details.

## Multi-window analysis and portability

Each run analyses one configured matrix resolution through five independent
genomic window scales: 25, 50, 75, 100, and 125 bins. The scale is always
expressed in bins, so a valid integer resolution changes physical sizes without
changing workflow, Python, or R code.

### Bin-derived arithmetic

For resolution `R`, window bins `W`, and `window_analysis.step_fraction` `F`:

```text
window_size_bp = R * W
step_bins      = ceil(W * F)
step_size_bp   = R * step_bins
```

The ceiling is mandatory: a 3.4-bin step becomes 4 bins before coordinates are
created. With the default `F = 0.20`, the 4 kb scales are:

| Window bins | Window size | Step bins | Step size |
|---:|---:|---:|---:|
| 25 | 100 kb | 5 | 20 kb |
| 50 | 200 kb | 10 | 40 kb |
| 75 | 300 kb | 15 | 60 kb |
| 100 | 400 kb | 20 | 80 kb |
| 125 | 500 kb | 25 | 100 kb |

At **10 kb**, the same tags remain `w25`, `w50`, `w75`, `w100`, and `w125`,
but physical windows become 250, 500, 750, 1,000, and 1,250 kb; step bins
remain 5, 10, 15, 20, and 25, with 50, 100, 150, 200, and 250 kb steps.

### KNN background presets

`window_analysis` is the sole source of genomic-window arithmetic:

```yaml
window_analysis:
  window_bins: [25, 50, 75, 100, 125]
  step_fraction: 0.20
```

`background_pool` is the producer: it selects cohorts and creates large KNN
background tables. SNAP-HiC-Diff is a strict consumer: it neither invokes nor
copies `background_pool`; it selects, validates, and consumes declared files.

`knn_empirical_pvalue.background_presets` requires these five exact mappings
for the default 4 kb configuration:

```yaml
knn_empirical_pvalue:
  background_presets:
    - {resolution_bp: 4000, window_bins: 25, step_bins: 5, file: /path/to/r4000_w25.tsv}
    - {resolution_bp: 4000, window_bins: 50, step_bins: 10, file: /path/to/r4000_w50.tsv}
    - {resolution_bp: 4000, window_bins: 75, step_bins: 15, file: /path/to/r4000_w75.tsv}
    - {resolution_bp: 4000, window_bins: 100, step_bins: 20, file: /path/to/r4000_w100.tsv}
    - {resolution_bp: 4000, window_bins: 125, step_bins: 25, file: /path/to/r4000_w125.tsv}
```

Each active tuple has exactly one preset. A referenced TSV must be non-empty
and include `comparison`, `chrom`, `window_start`, `window_end`, `EC`, `IB`,
`SP`, `MSE`, `SSIM`, `Spearman`, `tw_pvalue`, and `NCC`. Validation also checks
finite KNN rows, physical window width, step grid, and duplicate coordinates.
Missing presets fail during preflight with the expected tuple; missing or
malformed files fail before KNN allocation and identify the tuple and file.

### Per-scale selection and cross-scale union

Each `comparison x window_bins` branch is evaluated independently:

1. KNN calculates raw empirical P values for `MSE`, `d_SSIM`, `d_Spearman`, `tw_score`, and `d_NCC`.
2. A window is KNN-passing when **any** raw P value is `< 0.05`.
3. It is Meltron-passing when `meltron_best_score >= 5.0`.
4. A candidate must pass both tests in its own scale.

KNN selection uses the five raw values; no cross-window multiple-testing
correction is applied. A candidate requires no minimum number of genomic
window-scale supports: a valid one-scale candidate remains in the union.

Genomic `support_scale_count` is descriptive only. It is distinct from CWorld
IS scale-band support: the existing IS-change caller still requires exactly two
or more CWorld IS scale bands for a primary biological call.

### Output contract

All window-dependent paths are qualified by resolution and bin window scale:

```text
results/window_scales/r<resolution>/w<window_bins>/
results/multiscale/r<resolution>/<comparison>/all_scale_candidates.tsv
results/multiscale/r<resolution>/<comparison>/candidate_union.tsv
results/multiscale/r<resolution>/<comparison>/is_change/<comparison>_is_change_loci.tsv
results/multiscale/r<resolution>/<comparison>/is_change/figures/manifest.tsv
```

`all_scale_candidates.tsv` is the source-window evidence ledger; it preserves
the scale identity, Meltron score, and five raw KNN P values. `candidate_union.tsv`
is the merged region-level handoff to IS-change. Every configured scale must
complete before aggregation publishes those outputs; valid empty scales retain
header-only schemas.

## CWorld insulation scores for Meltron

CWorld IS is produced only for Meltron from the configured balanced `.cool`
inputs. Its configured source resolution must be an integer of at least 1,000
bp and must match the input cooler. The configured windows are 32,000, 60,000,
92,000, 120,000, 152,000, 180,000, 212,000, 240,000, 272,000, and 300,000 bp;
each uses `--ids = --is / 2`. A CWorld filename may show `is32001`, but the
manifest and final table retain requested distance `32000`.

Each sample has one `compute_IS_chrom` task per available chromosome, each
processing all ten windows sequentially, followed by one `merge_IS` task.
The tasks use `meltron.cworld_is.chrom_job_threads` and
`meltron.cworld_is.merge_job_threads`, respectively. Donau memory is
CPU-driven; adjust `meltron.cworld_is.chrom_job_threads`, not a memory field.
Matrices are retained in `tmp/cworld_matrix/<sample>/`; bedGraphs and manifests
in `results/cworld_insulation/<sample>/`; final Meltron input is
`results/IS/<sample>/<sample>_IS.tsv` with `chrom`, `viewpoint`, `IS_distance`,
and `IS`. One failed chromosome blocks only that sample merge and comparison.

Preview selected-sample cleanup first; add `--execute` only after reviewing its output, then submit the workflow:

```bash
python scripts/cleanup_cworld_is.py --outdir results --samples H1394 H1395 H1396 H1397 H1398 H1406 H1409 H1410 H1436 L021
python scripts/cleanup_cworld_is.py --outdir results --samples H1394 H1395 H1396 H1397 H1398 H1406 H1409 H1410 H1436 L021 --execute
dsub -s run_snakemake.sh
```

---

## Directory Structure

```
.
├── Snakefile                    # Snakemake workflow orchestrator
├── run_snakemake.sh             # Entry point script
├── config/
│   ├── config.yaml              # Main pipeline parameters
│   ├── chroms.yaml             # Chromosome list (26 chromosomes)
│   └── dsub_submit.sh         # HPC job submission wrapper
├── data/
│   ├── manifest.tsv            # Sample metadata (paths, conditions)
│   └── window_compare_pairs.yaml # Comparison pairs
├── rules/                       # Snakemake rule definitions
│   ├── validate_inputs.smk
│   ├── extract_chrom.smk
│   ├── extract_expected.smk
│   ├── compute_oe_matrix.smk
│   ├── compute_gap_qc.smk
│   ├── sliding_window.smk
│   ├── build_windows.smk
│   └── window_compare.smk
├── scripts/                     # Python execution scripts
│   ├── validate_inputs.py
│   ├── extract_expected.py
│   ├── compute_oe_matrix.py
│   ├── compute_gap_qc.py
│   ├── sliding_window.py
│   ├── build_windows.py
│   ├── compare_gaussian_smoothed_windows.py
│   ├── compare_iced_windows_tw.py
│   ├── compute_window_pair_covariates.py
│   ├── merge_window_pair_metrics.py
│   └── detect_breakpoint_by_nan_density.py
├── jupyter/                     # Exploratory notebooks
│   ├── gap_example.ipynb
│   ├── h1398_a11_breakpoint_explorer.ipynb
│   ├── window_comparison_visualization.ipynb
│   └── windows_normalize.ipynb
├── docs/                        # Documentation
│   ├── project_design.md       # Architecture & pipeline design
│   ├── rule_index.md            # All 8 rules reference
│   ├── snakemake_v8_30_0_standard.md  # Snakemake coding standards
│   ├── dsub_hpc_standard.md    # HPC/dSub usage guide
│   └── snakemake/wildcards.md  # Wildcard conventions
└── env/
    └── SNAP-HiC-Diff.yaml      # Conda environment
```

---

## Pipeline Stages

| Stage | Rule | Key Script | Output |
|-------|------|-----------|--------|
| 1. Validate inputs | `validate_inputs.smk` | `validate_inputs.py` | Validated `.mcool` files |
| 2. Extract chromosome | `extract_chrom.smk` | — | Per-chromosome cool files |
| 3. Extract expected | `extract_expected.smk` | `extract_expected.py` | Expected contact vectors (NPZ) |
| 4. O/E matrix | `compute_oe_matrix.smk` | `compute_oe_matrix.py` | O/E matrices (NPZ) |
| 5. Gap QC | `compute_gap_qc.smk` | `compute_gap_qc.py` | Gap region annotations |
| 6. Sliding window | `sliding_window.smk` | `sliding_window.py` | HDF5 (ln_oe, gap_filled, zscore, gaussian) |
| 7. Build windows | `build_windows.smk` | `build_windows.py` | Window coordinates (parquet) |
| 8. Compare | `window_compare.smk` | `compare_*.py` | Differential metrics (CSV) |

---

## Key Algorithms

### Gap-Filling (`sliding_window.py`)
Two-stage nearest-neighbor 2D interpolation across multiple mask scales. Handles chromosomal gaps in Hi-C data where contact information is absent.

### Z-Score Normalization (`sliding_window.py`)
Global robust z-score using **median + MAD** (median absolute deviation), resilient to outliers in sparse regions.

### Breakpoint Detection (`compute_gap_qc.py`)
Gradient ratio thresholding on per-window coverage profiles to identify abrupt transitions indicating structural gaps.

### NaN Density Analysis (`detect_breakpoint_by_nan_density.py`)
Scans upper-triangle diagonals for 5x5 zero-density patterns, identifying breakpoints at chromosomal boundaries.

### Tracy-Widom Statistics (`compare_iced_windows_tw.py`)
Extracts the largest eigenvalue from ICE-normalized matrices, computes TW statistic and p-value for significance testing.

---

## Configuration Reference

### `config/config.yaml`

| Section | Key Parameters |
|---------|---------------|
| `resolution` | Bin size in bp (e.g., `10000`) |
| `window_analysis.window_bins` | Genomic window widths in bins (`[25, 50, 75, 100, 125]` by default) |
| `window_analysis.step_fraction` | Fraction of a window used for the ceiling-rounded bin step (`0.20` by default) |
| `knn_empirical_pvalue.background_presets` | One external KNN background preset for every active resolution/window/step tuple |
| `gap_detection.gradient_threshold` | Coverage drop ratio to flag a gap |
| `gap_detection.max_consecutive_nan` | Max NaN bins before gap |
| `smoothing.sigma` | Gaussian kernel sigma |
| `parallel.threads_per_job` | CPU cores per job |
| `dsub.queue` | HPC queue (e.g., `arm`) |

### `data/window_compare_pairs.yaml`

Defines which sample windows to compare:

```yaml
pairs:
  - name: "rep1_vs_rep2"
    sample_a: "H1398_rep1"
    sample_b: "L128_rep1"
```

---

## HPC Deployment

The pipeline uses **dsub** for HPC job submission. See [`docs/dsub_hpc_standard.md`](docs/dsub_hpc_standard.md) for the complete guide.

Key points:

- `config/dsub_submit.sh` sets `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1` to prevent nested parallelism
- Jobs run on the `arm` queue with 10–80 cores per job
- Logs go to `%J.out` / `%J.err`
- Pipeline uses Snakemake `--dsub` profile for job coordination

---

## Development

### Snakemake Standards

All rules follow Snakemake v8.30.0 conventions documented in [`docs/snakemake_v8_30_0_standard.md`](docs/snakemake_v8_30_0_standard.md):

- Field order: `input` → `output` → `log` → `benchmark` → `threads` → `params` → `conda` → `shell`
- No f-strings; use `.format()` with `{{wildcard}}` escaping
- Wildcards extracted in shell blocks: `SAMPLE="{wildcards.sample}"`
- Bash vars after `_` or alphanumeric: `${{VAR}}`
- Benchmark required on every rule

### Adding a New Rule

1. Create `rules/new_rule.smk` following the standard template
2. Add corresponding Python script in `scripts/`
3. Register in `Snakefile`
4. Update `config/config.yaml` if new parameters needed
5. Document in `docs/rule_index.md`

---

## Documentation

| Document | Purpose |
|----------|---------|
| [`docs/project_design.md`](docs/project_design.md) | Pipeline architecture, 11 stages, data flow |
| [`docs/rule_index.md`](docs/rule_index.md) | All 8 Snakemake rules reference |
| [`docs/snakemake_v8_30_0_standard.md`](docs/snakemake_v8_30_0_standard.md) | Snakemake coding conventions |
| [`docs/dsub_hpc_standard.md`](docs/dsub_hpc_standard.md) | HPC job submission guide |
| [`docs/snakemake/wildcards.md`](docs/snakemake/wildcards.md) | Wildcard and shell placeholder conventions |

---

## Knowledge Graph

Interactive visualization of the codebase architecture is available. Run the dashboard:

```bash
cd /path/to/snap-hic-diff
node .understand-anything-plugin/packages/dashboard/serve.mjs dist
```

Then open the reported URL to explore layers, dependencies, and guided tours.

---

## License

Internal research pipeline. Contact authors for access.
