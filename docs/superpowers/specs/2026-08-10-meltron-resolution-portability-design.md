# Resolution-portable CWorld-to-Meltron coordinates

## Goal

Make the CWorld insulation-score export compatible with Meltron at any configured source resolution.  A 4 kb CWorld bedGraph must export Meltron viewpoints on 4 kb starts, a 10 kb bedGraph on 10 kb starts, and no code path may assume a particular resolution.

The repair also regenerates the derived IS and Meltron outputs for `H1398_vs_L134` so its final significant-window count is statistically interpretable.

## Problem

The converter currently writes each bedGraph bin midpoint as `viewpoint`.  For a 4 kb bin this produces `2000, 6000, ...`.  Meltron determines matrix resolution with the greatest common divisor of `viewpoint`; it therefore infers 2 kb and doubles each interval's expected number of IS values.  The resulting apparent missingness is about 50%, so all Meltron scores become `NA`.

## Coordinate contract

`scripts/convert_cworld_IS_to_meltron_IS.py` is the single adapter between CWorld bedGraph and Meltron TSV.

- Input data rows are `chrom`, `start`, `end`, `insulation`.
- The modal complete bedGraph width is the observed resolution `W`.
- Rows whose width differs from `W` are accepted only when they are terminal incomplete bins; all other mismatches fail with a contextual error.
- Every retained row must have `start % W == 0`.
- Meltron `viewpoint` is the bedGraph `start`, never the bin midpoint.
- The converter receives the workflow's `source_resolution_bp` and fails if `W` differs from it.
- `(chrom, viewpoint, IS_distance)` must be unique in the produced TSV.
- Existing atomic output semantics remain: failures do not leave a partial `results/IS/<sample>/<sample>_IS.tsv`.

Using starts gives Meltron `GCD(viewpoint) == W` for all positive integer source resolutions, including 4 kb, 10 kb, and non-round test values.

## Workflow integration

`compute_IS_chrom` remains responsible for CWorld calculation and per-chromosome manifests.  It does not change computational topology or thread allocation.

`merge_IS` passes `CWORLD_IS.source_resolution_bp` to the converter.  Its output remains the per-sample Meltron input TSV.

For the real-data validation, only derived outputs are regenerated:

1. Force `merge_IS` for `H1398` and `L134` to rewrite their IS TSVs from existing manifests and bedGraphs.
2. Force `meltron_per_comparison` for `H1398_vs_L134` to rebuild its interval, directional, and merged outputs.
3. Count final calls using `meltron_is_significant == 1`; report candidate count, callable-score count, and final significant count separately.

CWorld matrix generation and `compute_IS_chrom` are not rerun solely for this coordinate repair.

## Error handling

The workflow fails before Meltron when a manifest is incomplete, a bedGraph has incompatible full-bin width, coordinates are not resolution-aligned, or output keys are duplicated.  These errors name the sample, chromosome, source file, observed width, and configured resolution.

The post-run report treats an all-`NA` Meltron score set as an invalid/uncallable analysis, not as evidence of zero biological rearrangements.

## Verification and acceptance

Automated tests cover:

- 4 kb, 10 kb, and a non-predeclared resolution fixture; each exports bin starts and has the correct inferred Meltron resolution.
- A width/configuration mismatch and an off-grid start, both of which fail.
- Preservation of the incomplete-terminal-bin behavior.
- Rejection of duplicate Meltron coordinates and existing manifest integrity failures.

The actual `H1398_vs_L134` rerun passes when both IS tables contain ten configured distances with viewpoints divisible by 4,000; Meltron no longer has the systematic approximately 0.50 missingness artefact; and each `NA` score, if any, is attributable to an interval's real missingness rather than coordinate resolution inference.  The biological number of significant windows is not predetermined.

## Isolation

Implementation occurs only on branch `codex/meltron-resolution-portability` in its isolated worktree.  Existing uncommitted changes in the primary checkout are not staged, moved, or incorporated.
