# Multi-window KNN, Meltron, and IS-change analysis

## Goal

Extend the workflow from one fixed genomic window to five resolution-derived
window scales in a single run. Each scale is analysed independently by the
existing matrix-comparison, KNN empirical-P-value, and Meltron stages. Passing
windows from all five scales are then combined and sent to the existing
multi-scale CWorld insulation-score (IS) change caller.

The workflow must remain resolution-portable: changing the configured matrix
resolution must change physical window and step sizes through bin arithmetic,
without changing Python, Snakemake, or R code.

This design deliberately keeps KNN background production outside this
repository. The separate `background_pool` workflow produces background
tables; SNAP-HiC-Diff consumes five explicitly declared, compatible presets and
fails clearly when a preset is absent or incompatible.

## Scope and invariants

- One configured matrix resolution is analysed per workflow run.
- Five window scales are expressed as bin counts, not hard-coded base-pair
  lengths.
- The default scale set is `[25, 50, 75, 100, 125]` bins.
- The step is 20% of each window and is always rounded upward to a whole bin.
- KNN and Meltron run independently for every `comparison x window_bins`.
- The ten configured CWorld IS distances are generated once per sample and
  shared by all five Meltron/window-scale branches.
- Existing KNN, Meltron, and IS-change thresholds are not relaxed or replaced.
- A final candidate does not need support from two or more window scales. The
  union of all five scales is intentional.
- Existing unqualified single-window results are neither deleted nor used by
  the new multi-window workflow.

Generating KNN backgrounds, running several matrix resolutions in one
invocation, introducing a multiple-window-testing correction, and changing the
biological IS-change caller are out of scope.

## Window arithmetic

Let:

- `R` be the configured matrix resolution in base pairs;
- `W` be one configured window size in bins; and
- `F` be the configured step fraction, default `0.20`.

Then:

```text
window_size_bp = R * W
step_bins      = ceil(W * F)
step_size_bp   = R * step_bins
```

`ceil` is mandatory. Fractional bins must never reach coordinate construction.
For example, a calculated step of 3.4 bins is 4 bins.

At 4 kb resolution, the defaults produce:

| Window bins | Window size | Step bins | Step size |
|---:|---:|---:|---:|
| 25 | 100 kb | 5 | 20 kb |
| 50 | 200 kb | 10 | 40 kb |
| 75 | 300 kb | 15 | 60 kb |
| 100 | 400 kb | 20 | 80 kb |
| 125 | 500 kb | 25 | 100 kb |

At 10 kb resolution, the same configuration produces 250, 500, 750, 1,000,
and 1,250 kb windows while retaining steps of 5, 10, 15, 20, and 25 bins. No
physical size in either example may be encoded in workflow logic.

## Configuration contract

The single-window `window_size` and `step_size` settings are replaced by one
authoritative multi-window section:

```yaml
resolution: 4000

window_analysis:
  window_bins: [25, 50, 75, 100, 125]
  step_fraction: 0.20

knn_empirical_pvalue:
  background_presets:
    - resolution_bp: 4000
      window_bins: 25
      step_bins: 5
      file: /path/to/r4000_w25_background.tsv
    - resolution_bp: 4000
      window_bins: 50
      step_bins: 10
      file: /path/to/r4000_w50_background.tsv
    - resolution_bp: 4000
      window_bins: 75
      step_bins: 15
      file: /path/to/r4000_w75_background.tsv
    - resolution_bp: 4000
      window_bins: 100
      step_bins: 20
      file: /path/to/r4000_w100_background.tsv
    - resolution_bp: 4000
      window_bins: 125
      step_bins: 25
      file: /path/to/r4000_w125_background.tsv
```

The workflow derives the expected `(resolution_bp, window_bins, step_bins)`
tuples and requires exactly one preset for each tuple. Duplicate presets and
unused presets for the active resolution are configuration errors. A preset's
metadata is the producer/consumer contract; the table coordinates provide an
additional physical-window and step-grid check.

The default list has five entries, but the implementation should iterate the
configured list rather than contain five copied rule bodies. Validation
requires positive integer resolution and window-bin values, a finite positive
step fraction, unique window-bin values, and a derived step of at least one
bin.

Legacy `window_size` and `step_size` are removed from the active configuration
instead of being maintained as a second source of truth. A helpful preflight
error explains the new settings when the multi-window section is missing.

## Responsibility boundary for KNN backgrounds

The two repositories remain independent:

- `background_pool` owns background cohort selection and background-table
  generation for each resolution/window tuple.
- SNAP-HiC-Diff owns preset selection, validation, KNN calibration of target
  comparisons, and downstream result aggregation.

SNAP-HiC-Diff does not invoke `background_pool`, copy its rules, or attempt an
automatic fallback calculation. Large background inputs therefore remain
outside this portable analysis repository, while the required background
interface is explicit and testable.

Before a KNN job can use a preset, one validation job per window scale checks:

- the declared file exists, is readable, and is non-empty;
- provenance/coordinate columns `comparison`, `chrom`, `window_start`, and
  `window_end` exist;
- KNN feature/score columns `EC`, `IB`, `SP`, `MSE`, `SSIM`, `Spearman`,
  `tw_pvalue`, and `NCC` exist;
- numeric cleaning leaves at least `knn_k` valid rows;
- every complete window has `window_end - window_start == resolution_bp *
  window_bins`;
- positive start differences within each `comparison x chrom` group are
  integer multiples of the expected `step_size_bp`, allowing windows removed
  by upstream QC; and
- the table has no duplicate `(comparison, chrom, window_start, window_end)`
  rows. Reuse of the same genomic coordinates by different background
  comparisons is valid and expected.

Validation writes a small scale-qualified audit artifact containing the
preset path, expected tuple, row counts, and observed coordinate properties.
The corresponding KNN job depends on that artifact. Missing preset mappings
fail during workflow preflight; malformed tables fail before expensive KNN
work and name the offending tuple and file.

A background table with a 300 kb physical span is not automatically valid for
every way of producing 300 kb. Its declared resolution and window-bin tuple
must match the active preset. The existing single 300 kb table can therefore
serve at most one declared scale and cannot substitute for the other four.

## Workflow topology

Resolution-level products remain shared:

```text
manifest and comparisons
          |
          +--> balanced/expected/OE products (once per sample and resolution)
          |
          +--> ten-distance CWorld IS products (once per sample)
          |
          +--> for each configured window_bins
                    |
                    +--> windows and scale-specific QC
                    +--> sliding-window matrices and five comparison metrics
                    +--> matching KNN background validation
                    +--> KNN empirical P values
                    +--> Meltron intervals and scores using shared IS
                    +--> scale-specific jointly passing windows
          |
          +--> concatenate five scale-specific candidate tables
                    |
                    +--> coordinate union with scale provenance
                              |
                              +--> existing IS-change caller
```

Every window-dependent Snakemake wildcard and script invocation receives
`window_bins` and its derived window/step values explicitly. Module-level
single-window constants must not silently leak into scale-dependent jobs.

All five branches are required for the final aggregate. A failure in one scale
blocks aggregation rather than publishing a partial result. A valid scale with
zero passing windows produces a header-only candidate table; if all scales
have zero candidates, the final union and downstream IS-change outputs are
valid empty artifacts rather than workflow errors.

## Output and identity contract

All window-dependent artifacts, logs, benchmarks, and validation reports are
qualified by active resolution and window bins:

```text
results/window_scales/r4000/w25/
results/window_scales/r4000/w50/
results/window_scales/r4000/w75/
results/window_scales/r4000/w100/
results/window_scales/r4000/w125/
```

Within each directory, the workflow retains the existing logical stages such
as windows, window comparison, KNN empirical P value, and Meltron. Shared
balanced matrices and IS data remain in their existing single-copy locations.

The stable cross-stage window identity is:

```text
(resolution_bp, window_bins, chrom, window_start, window_end)
```

`window_id` may remain as a readable per-scale field, but it is not globally
unique and must not be used alone when joining or aggregating scales. Meltron
preparation and results carry the source `window_start` and `window_end`
separately from Meltron's algorithm-facing `start` and `end`. KNN and Meltron
join on the full scale-qualified source identity. Their algorithm-facing
intervals retain the current start-consistency and end-delta audit rather than
silently redefining the source window.

Old unqualified output paths are not compatibility aliases. They are left on
disk for provenance but are not dependencies of the new DAG, preventing stale
300 kb files from satisfying new multi-window targets.

## Per-scale statistical selection

Each `comparison x window_bins` branch preserves the current statistics:

1. Derive the five KNN score directions for `MSE`, `d_SSIM`, `d_Spearman`,
   `tw_score`, and `d_NCC` from the existing raw metrics and configuration.
2. Compute the five raw KNN empirical P values with the matching background.
3. Mark a window as KNN-passing when any one of the five values is `< 0.05`.
4. Mark it as Meltron-passing when `meltron_best_score >= 5`.
5. Select the window only when both KNN and Meltron conditions pass.

There is no additional correction across the five window scales. Documentation
must describe the resulting union as candidate screening, not as a new
family-wise or globally corrected significance statement.

## Cross-scale union

Scale-specific selected tables are first concatenated without information
loss. The audit table records one row per source window and adds
`resolution_bp`, `window_bins`, `window_size_bp`, `step_bins`, and
`step_size_bp`.

Selected intervals are then sorted by chromosome and coordinates and merged
when they overlap or are BED-adjacent. No minimum number of window scales is
required. Each merged region records:

- `support_window_bins`: sorted unique supporting bin counts;
- `support_scale_count`: number of unique supporting window scales;
- `source_window_count`: number of contributing source windows;
- scale-qualified member window and Meltron identifiers;
- `best_meltron_score`: maximum contributing score; and
- the minimum empirical P value for each of the five KNN metrics.

Recommended aggregate outputs are:

```text
results/multiscale/r4000/<comparison>/all_scale_candidates.tsv
results/multiscale/r4000/<comparison>/candidate_union.bed
results/multiscale/r4000/<comparison>/candidate_union.tsv
results/multiscale/r4000/<comparison>/is_change/
```

The source-window table is the evidence ledger; the union table is the
region-level input to IS-change analysis. This separation prevents merged
regions from hiding which window sizes generated the evidence.

## Integration with the existing IS-change caller

Only candidate construction changes. The current IS-change scoring and
classification remain intact:

- ten configured CWorld IS distances are reused;
- existing small/medium/large IS scale bands are retained;
- the current requirement for support from at least two IS scale bands remains
  the primary-call rule;
- existing transition, focal-pulse, mixed, and broad-shift classifications are
  retained; and
- existing local-extremum refinement/recovery behavior remains unchanged.

`support_scale_count` refers to genomic window scales and is descriptive. It
must not be confused with the number of CWorld IS scale bands, which remains a
biological-calling criterion downstream.

## Failure handling

- Invalid window arithmetic or duplicate window-bin entries fail during
  configuration preflight.
- A missing exact background preset identifies the expected resolution,
  window bins, step bins, and configuration key.
- A missing or malformed background file blocks only after producing a focused
  validation error and before KNN allocation.
- KNN and Meltron inputs with duplicate scale-qualified identities or
  coordinate disagreement fail instead of being many-to-many joined.
- Scale-dependent outputs are written atomically so a failed rerun cannot
  masquerade as a completed scale.
- Aggregation refuses an incomplete set of configured scale outputs.
- Empty-but-valid candidate results preserve schemas and continue cleanly.

## Documentation

README updates must describe:

- the five bin-derived window defaults and the ceiling step formula;
- 4 kb and 10 kb examples;
- the separate responsibilities of `background_pool` and SNAP-HiC-Diff;
- the background preset schema, required table columns, and validation errors;
- the KNN, Meltron, cross-scale union, and IS-change decision sequence;
- the distinction between window-scale support and CWorld IS-scale support;
- scale-qualified and aggregate output locations; and
- the fact that the union is candidate screening with unchanged raw thresholds.

## Verification

Focused automated tests cover:

1. exact window and step arithmetic at 4 kb and 10 kb;
2. ceiling behavior for non-integer products, including 3.4 bins becoming 4;
3. arbitrary valid integer resolutions without code changes;
4. rejection of invalid, duplicate, or legacy-only window configuration;
5. exact background preset selection for every configured tuple;
6. missing, duplicate-logical-row, undersized, malformed, wrong-span, and
   off-grid background fixtures, including valid coordinate reuse by separate
   background comparisons;
7. collision-free scale-qualified paths and stable identities;
8. preservation of the current five KNN P-value and Meltron score thresholds;
9. scale-specific KNN/Meltron joins and zero-candidate behavior;
10. cross-scale overlap/adjacency merging, provenance annotations, and P-value
    and Meltron summaries;
11. a single-scale-supported region remaining eligible; and
12. regression of the existing IS-change configuration, calling, extremum,
    and plotting suites.

A small synthetic integration fixture exercises the Snakemake DAG through all
five scale branches and the aggregate target. Dry runs are performed for 4 kb
and 10 kb configurations to prove that target paths and physical sizes derive
from configuration rather than hard-coded base-pair constants. Heavy real-data
KNN, CWorld, and Meltron computations are not required for the unit-level
acceptance suite.

## Acceptance criteria

- One normal Snakemake invocation expands all five configured window scales for
  every comparison.
- Every scale consumes exactly one matching, successfully validated KNN
  background preset.
- Changing from 4 kb to 10 kb or another valid integer resolution requires
  configuration changes only.
- No scale overwrites or reuses another scale's artifacts, logs, or identities.
- Current KNN raw-P-value, Meltron score, and IS-change thresholds are
  unchanged.
- Passing windows from any one or more scales enter the annotated coordinate
  union.
- The union feeds the existing ten-distance, multi-band IS-change analysis and
  retains complete source-window provenance.
- Valid empty results complete with stable schemas; incomplete scale execution
  never publishes a final union.
- Updated README examples, focused tests, full relevant regression tests, and
  4 kb/10 kb Snakemake dry runs pass.

## Isolation

Design and implementation occur on a `codex/` feature branch. Existing
uncommitted runtime outputs and unrelated primary-worktree changes are not
staged, deleted, or incorporated into feature commits.
