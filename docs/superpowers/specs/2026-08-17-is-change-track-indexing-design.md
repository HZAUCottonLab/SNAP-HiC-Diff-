# IS-change Track Indexing Design

## Purpose

Make full IS-change calling practical after repaired candidate selection admits
the real `H1398_vs_L134` candidate set.  This is a performance-only change:
it removes repeated scans of the complete IS tracks while preserving every
scientific decision in the existing caller.

## Evidence

The repaired candidate interface produces 7,329 selected windows and 1,364
merged KNN-core regions for `H1398_vs_L134`.  Each sample's IS track has
5,742,821 rows.  The current caller reloads each track once, but then, for
every candidate region and every configured scale, filters all 5.7 million
rows to obtain one chromosome/scale subset.  A 2-core Donau validation showed
stable memory use below 1 GB but made slow, single-core progress.  The login
node's `137` termination was therefore a local-resource limit, not evidence
that the IS algorithm or the repaired candidate coordinates are invalid.

## Decision

After the existing `read_is_track` validation succeeds, construct a lookup
whose key is `(chrom, IS_distance)` and whose value contains exactly that
chromosome/scale pair's `start` and `IS` columns.  Build it once per sample.

`_event_rows_for_region` will retrieve the relevant two small lookup tables
for each scale rather than evaluating a boolean filter over the whole track.
It will perform the same left merge by bin start into the existing context
grid.  Candidate regions, balanced-bin availability, context extent, masking,
score calculation, percentile thresholds, event classification, clustering,
ranking, and atomic output writing remain unchanged.

The lookup is created only after validation, so malformed IS input retains its
current failure behaviour.  The full raw track object is discarded once its
lookup has been built, so the lookup does not intentionally retain two copies
of the full input table.

## Alternatives considered

1. **Chromosome-by-chromosome streaming:** lowers memory further, but requires
   rereading each 236 MB IS file for every chromosome and adds I/O-heavy
   control flow.
2. **One dsub job per chromosome:** may shorten elapsed time, but requires
   independent output merging and complicates atomic all-or-nothing results.
3. **Chosen: one-time in-memory lookup:** keeps a single input pass and the
   caller's current serial output contract while eliminating the repeated
   global scans.

## Validation

- Add a focused unit test for the lookup's chromosome/scale partitioning.
- Retain the caller's existing integration tests to prove output schemas,
  masks, ranking, and score behaviour remain intact.
- Run all `test/is_change` tests only through a Donau `dsub` compute-node job.
- Run a real `H1398_vs_L134` candidate build and IS-change replay only through
  `dsub`, initially with `cpu=2` (10 GB on `arm`).  Record the final primary,
  preliminary, and masked row counts and scheduler `MEM_MAX_SUM`.
- Use observed CPU and memory use to choose the eventual Snakemake resource
  request; do not infer a biological result from runtime or memory metrics.

## Non-goals

- Changing KNN/Meltron candidate gates or their identity repair.
- Changing any IS value, log transformation, missing-bin rule, scale band,
  threshold, event label, or TAD-boundary interpretation.
- Altering Meltron coordinate files.
- Introducing chromosome-parallel execution in this change.
