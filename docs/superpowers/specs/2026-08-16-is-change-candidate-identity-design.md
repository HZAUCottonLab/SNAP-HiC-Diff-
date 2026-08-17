# IS-change candidate identity design

## Purpose

Ensure that Meltron-filtered windows enter multiscale IS-change analysis even
when Meltron reports an interval endpoint using a different bin convention from
the original KNN sliding window.

This is a narrow repair to candidate-region construction.  It does not alter
the KNN or Meltron thresholds, IS-change scoring, missing-bin policy, region
context, or downstream interpretation.

## Observed interface mismatch

For `H1398_vs_L134`:

- 28,981 KNN windows have at least one raw empirical P value below 0.05.
- The Meltron merged table contains 13,370 prefiltered windows.
- 7,329 windows have `meltron_best_score >= 5.0`.
- Every Meltron row has a unique matching KNN row by `chrom + window_id`.
- Every paired row has the same start coordinate, but Meltron reports an end
  coordinate exactly 4 kb beyond the KNN end coordinate.

The existing selector uses `chrom + start + end` as its join key.  Consequently
it selects zero windows: the endpoint representation is being used as window
identity.

## Decision

### Window identity and coordinates have distinct roles

- **Identity:** `chrom + window_id` is the required, unique KNN--Meltron join
  key.
- **Candidate core interval:** KNN `window_start` and `window_end` are the
  canonical BED-style `[start, end)` coordinates.  They define the region
  merged from selected windows and the interval in which an IS-change summit
  may be reported.
- **Meltron interval coordinates:** they are diagnostic metadata only.  They
  must not gate the join and must not replace the canonical KNN core interval.

The selector retains the declared final gate: at least one KNN raw empirical P
value below the configured threshold and `meltron_best_score` at or above the
configured threshold.  In the observed comparison, all 7,329 Meltron-score
passing rows already satisfy the KNN condition, so all 7,329 enter IS-change.

### Region geometry

Selected KNN core intervals are merged per chromosome when they overlap or are
BED-adjacent.  For the observed comparison, the 7,329 selected windows produce
1,364 merged candidate regions.

The caller continues to read 300 kb of IS context on both sides of each merged
core region.  Context supports flank calculations only; it never expands the
set of reportable IS-change loci.  Only events whose summits or intervals
overlap the merged KNN core are reported.

## Validation and failure behaviour

Before filtering, candidate construction must verify that both inputs contain
`chrom` and `window_id`, and that each `chrom + window_id` key is unique within
each input.  It must reject duplicated or missing identity keys rather than
silently create many-to-many matches.

The selector must retain KNN start/end coordinates.  It should audit Meltron
coordinate differences in its log: matching starts are expected; endpoint
differences are reported but allowed because Meltron may use an inclusive final
IS bin.  A start-coordinate mismatch is a true mapping error and must fail.

## Tests

- A Meltron row with the same `chrom + window_id`, the same start, and an end
  one 4 kb bin beyond KNN must be selected and retain KNN coordinates.
- Exact-coordinate rows must remain supported.
- Duplicate identity keys, unmatched IDs, and shifted starts must fail with a
  clear error.
- Existing overlap, adjacency, chromosome-boundary, and positive-gap merging
  tests remain unchanged.
- The real `H1398_vs_L134` candidate build must yield 7,329 selected windows
  and 1,364 merged regions before IS-change calling.

## Out of scope

- Changing the KNN raw-P or Meltron-score thresholds.
- Subtracting 4 kb from Meltron files or changing Meltron itself.
- Treating the 4 kb endpoint convention as IS-change evidence.
- Altering IS curve scoring, masking, scale integration, or TAD-boundary
  annotation.
