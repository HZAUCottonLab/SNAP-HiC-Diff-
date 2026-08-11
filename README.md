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

# Resolution (bp)
resolution: 4000

# Window size & step (bp)
window_size: 300000
step_size: 30000

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
| `window_size` | Genomic window width (bp) |
| `step_size` | Sliding step (bp) |
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
