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

# Multiscale window settings
window_analysis:
  window_bins: [25, 50, 75, 100, 125]
  step_fraction: 0.2

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

### Multiscale window-analysis contract

At the 4 kb production resolution, the configured window scales are:

| Label | Window bins | Window bp | Step bins | Step bp |
| --- | ---: | ---: | ---: | ---: |
| `w25` | 25 | 100,000 | 5 | 20,000 |
| `w50` | 50 | 200,000 | 10 | 40,000 |
| `w75` | 75 | 300,000 | 15 | 60,000 |
| `w100` | 100 | 400,000 | 20 | 80,000 |
| `w125` | 125 | 500,000 | 25 | 100,000 |

The five production absolute background paths live only in
`config/config.yaml` and are read in place; the workflow does not copy, link,
regenerate, or edit these backgrounds. The scale root
`results/window_scales/r4000/w*/` holds the scale-dependent windows,
comparisons, KNN results, diagnostics, and Meltron outputs.

For each comparison, `results/multiscale/r4000/<comparison>/all_scale_candidates.tsv`
is the evidence ledger. `candidate_union.tsv` and `candidate_union.bed` are the
any-scale union and retain support columns that show the contributing scales
and windows. Primary IS-change calls still require support from at least two
CWorld IS scale bands. A valid header-only empty result is a successful result,
not a failed workflow.

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

One manual `dsub -s run_snakemake.sh` starts one full DAG; every scale branch
whose dependencies are ready is eligible to run concurrently.

See [docs/dsub_hpc_standard.md](docs/dsub_hpc_standard.md) for HPC setup details.

## 3DNetMod TADs and boundary regions

The default workflow runs raw and balanced 3DNetMod D calls at 4 kb, 8 kb,
and 20 kb for every manifest sample. Raw calls use positive pixel counts;
balanced calls use `count * weight1 * weight2` rescaled by one constant to raw-count depth. No chromosome is excluded by
default.

Build or refresh all Conda prefixes on the login node before submission:

```bash
bash scripts/conda_envs.sh bootstrap
CONDA_OFFLINE=true bash scripts/conda_envs.sh verify
```

Compute nodes do not create environments or access Conda channels.

Stable outputs are:

```text
results/tads/{sample}.{4kb|8kb|20kb}.{raw|balanced}.3dnetmod.bed
results/tads/{sample}.{4kb|8kb|20kb}.{raw|balanced}.3dnetmod.bedpe
results/tad_boundaries/{sample}.{4kb|8kb|20kb}.{raw|balanced}.boundaries.bed
results/tad_boundaries/{sample}.multires.{raw|balanced}.boundaries.bed
```

Run one sample by targeting its two merged boundary files; Snakemake will pull
the six required TAD calls:

```bash
snakemake --use-conda --conda-prefix "$PWD/.snakemake/conda" \
  --configfile config/config.yaml -- \
  results/tad_boundaries/H1394.multires.raw.boundaries.bed \
  results/tad_boundaries/H1394.multires.balanced.boundaries.bed
```

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

## IS-change and TAD-boundary annotation

This optional downstream layer annotates IS-change loci with a single,
provenance-aware collection of TAD-boundary events. It is a boundary-level
analysis: it does not infer complete-TAD split, fusion, or nested-domain
topology. For `H1398_vs_L134`, direction is always A-to-B
(`H1398` to `L134`): `gained_in_B` is present in L134 only and `lost_in_B` is
present in H1398 only.

Position and strength are deliberately separate axes. `position_state` reports
whether a boundary is `conserved`, `gained_in_B`, `lost_in_B`, `shifted`,
`complex`, or `ambiguous`. `strength_state` is evaluated only where useful
from the multiscale IS curves: a conserved boundary can therefore be
`weakened_in_B` or `strengthened_in_B` while its position remains conserved.
IS and 3DNetMod both use the same Hi-C experiments, so their agreement is
cross-method consistency evidence, not an independent biological replicate.

Raw and balanced 3DNetMod boundary outputs remain intact. The annotation stage
uses a downstream union only: overlapping calls retain `both` provenance and
unmatched calls retain `raw_only` or `balanced_only` provenance. The raw and
balanced calls are not pooled as replicates. Localization jitter is the 95th
percentile of eligible technical offsets (at least 100 offsets); otherwise it
falls back to 20 kb. The primary shifted-boundary search radius is 100 kb, and
60/80/100 kb calls are retained as a sensitivity check. A shifted event is
associated with an IS locus only through either endpoint support interval,
never through the corridor between its old and new positions.

The primary and `single_band_preliminary` IS-locus tables are annotated by the
same algorithm but are published separately. Use the primary table for main
summaries; preliminary calls are a sensitivity output and must not be pooled
with it.

Stable annotation outputs are:

```text
results/tad_boundaries/{sample}.multires.unified.boundaries.tsv
results/tad_boundaries/{sample}.boundary_jitter_offsets.tsv
results/tad_boundaries/{sample}.annotation_views.tsv
results/tad_boundaries/{sample}.tad_unevaluable.tsv
results/tad_boundary_events/{comparison}.boundary_events.tsv
results/tad_boundary_events/{comparison}.boundary_event_stability.tsv
results/tad_boundary_events/{comparison}.matching_qc.tsv
results/tad_boundary_events/{comparison}.boundary_strength_by_scale.tsv
results/tad_boundary_events/{comparison}.boundary_strength_summary.tsv
results/tad_boundary_events/{comparison}.annotation_qc.tsv
results/tad_boundary_events/{comparison}.annotation_summary.tsv
results/multiscale/r4000/{comparison}/is_change/{comparison}_is_change_loci.boundary_annotated.tsv
results/multiscale/r4000/{comparison}/is_change/{comparison}_is_change_loci_preliminary.boundary_annotated.tsv
results/multiscale/r4000/{comparison}/is_change/figures/manifest.tsv
results/tad_boundary_events/{comparison}.primary.permutation.tsv  # optional
results/tad_boundary_events/{comparison}.preliminary.permutation.tsv  # optional
```

Target the H1398-to-L134 annotated primary table and its diagnostic figure
manifest as follows. This schedules its required TAD, IS, unification,
matching, annotation, and plotting dependencies; it does not claim that data
are available or that a production run has completed.

```bash
CONDA_OFFLINE=true snakemake --use-conda --conda-prefix "$PWD/.snakemake/conda" \
  --configfile config/config.yaml -- \
  results/multiscale/r4000/H1398_vs_L134/is_change/H1398_vs_L134_is_change_loci.boundary_annotated.tsv \
  results/multiscale/r4000/H1398_vs_L134/is_change/figures/manifest.tsv
```

The set-level permutation analysis is disabled by default and is not a
locus-level P value. Enable it only for a deliberate sensitivity analysis:

```yaml
is_boundary_annotation:
  permutation:
    enabled: true
```

For a deterministic visual audit, use the annotation-mode figure manifest's
`locus_id` (not the non-unique `region_id`) together with its
`boundary_event_id`; an empty event ID explicitly records no endpoint event.
The plot selector retains one deterministic representative for every available
`position_state` before it fills the global top-N list. Filter to simple
classes and inspect up to 20 loci from every sufficiently populated class.
Record `boundary_event_id`, reviewer, `pass/fail`, and a note in
`results/tad_boundary_events/{comparison}.visual_review.tsv`. This manual QC
artifact is neither a Snakemake completion dependency nor a statistical test.

The rationale and implementation contract are in the
[design](docs/superpowers/specs/2026-09-03-is-change-tad-boundary-annotation-design.md),
[implementation plan](docs/superpowers/plans/2026-09-03-is-change-tad-boundary-annotation.md),
[ADR-0001](docs/adr/0001-3dnetmod-raw-and-balanced-signals.md), and
[ADR-0002](docs/adr/0002-downstream-unified-tad-boundary-view.md).

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
