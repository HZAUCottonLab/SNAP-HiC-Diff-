# Multiscale IS-change locus screening design

## Purpose and scope

Given the final KNN-plus-Meltron candidate windows for one A-versus-B comparison,
identify locations or intervals where the two samples' multiscale insulation-score
(IS) curves show large, structured changes.  The primary output is named an
**IS-change locus** or **IS-change interval**.

This stage does not call TAD boundaries, does not call structural variants, and
does not claim population-level statistical significance.  There is one biological
replicate per material.  A later workflow may annotate IS-change loci with TAD
boundary calls, SV calls, and contact-map evidence, but those inputs must not alter
the initial IS-change call.

## Inputs

- Final candidate windows: the 7,329 windows that pass both the KNN criterion
  (at least one of five empirical P values below 0.05) and
  `meltron_best_score >= 5.0`.
- A and B cworld IS tracks at 4 kb bins and ten IS windows:
  32, 60, 92, 120, 152, 180, 212, 240, 272, and 300 kb.
- Per-sample bin-validity information from the balanced `.cool` files and the IS
  tracks themselves.

The cworld IS value is already a normalized log2 insulation score.  The workflow
must not take another log2 transform.  It may robustly standardize *differences*
within a chromosome and scale for cross-scale comparability.

## Candidate-region construction

1. Filter to the final 7,329 candidate windows before merging.
2. Within each chromosome, sort by `window_start`, then merge any intervals for
   which `next_start <= current_end`.  This merges overlapping and BED-adjacent
   intervals only; `merge_gap_bp` is fixed to zero by default.
3. Never merge windows across chromosomes or across a positive gap.
4. Retain a region-level provenance record: `region_start`, `region_end`, number
   of constituent windows, constituent IDs, per-metric minimum empirical P value,
   and maximum Meltron score.  P values and Meltron scores are summaries only;
   they are not recombined with IS evidence.
5. Read up to 300 kb of flanking IS context on each side for calculation.  Only
   positions within the merged candidate region may be reported as IS-change
   results.  A region may contain zero, one, or multiple IS-change events.

## Missing-bin policy

Keep genomic coordinates intact.  Do not delete rows, impute missing IS, replace
missing values with zero, or smooth across a missing interval.

- `IS == 0` is a valid normalized value and must not be masked.
- A bin is unavailable at a scale if either sample has non-finite IS there, or if
  the underlying balanced bin is invalid in either sample.
- The central bin of a called event must be available in both samples.
- A local three-bin amplitude calculation requires at least two valid bins.
- Each flank used for transition or pulse calculations requires at least 80%
  observed bins.  A run of two or more missing bins in a required flank makes
  that scale-position `not_evaluable`.
- A/B-specific empty bins are recorded as a separate coverage-state annotation;
  they never contribute numeric IS-change evidence.

This policy avoids both false signal from missingness and excessive loss of
large-scale calls caused by demanding a completely gap-free 300 kb neighborhood.

## Per-scale IS evidence

For an available bin `i` and IS scale `s`, compute:

`delta_s(i) = IS_B,s(i) - IS_A,s(i)`.

Within each chromosome and scale, robustly center and scale the valid delta values
using their median and MAD:

`Z_s(i) = (delta_s(i) - median(delta_s)) / (1.4826 * MAD(delta_s))`.

This is a descriptive background normalization, not a biological P value.

For each valid position, calculate three complementary properties.

1. **Local amplitude, A**: the median of `abs(Z)` across the central bin and its
   immediate genomic neighbors.  It asks whether A and B differ strongly here.
2. **Left-right transition, C**: the absolute difference between left- and
   right-flank medians of `Z`.  The flank size is
   `h_s = max(2 bins, round(s / (4 * 4 kb)))`.  It asks whether the A/B difference
   changes state at this position.
3. **Local pulse, P**: the absolute difference between the central `Z` value and
   the median of the two flanks.  It asks whether the position is a local peak or
   trough relative to its surroundings.  At least one adjacent valid bin must show
   the same change direction, so a singleton spike is not called.

Convert A, C, and P independently to chromosome-and-scale empirical percentiles
over valid bins.  Initial high-evidence thresholds are the 95th percentile for A,
C, and P; a 90th-percentile amplitude threshold is used to identify broad changes.
These ranks are configurable descriptive cutoffs, not significance tests.

## Event types within one IS scale

Evidence is classified by curve geometry before it is ranked.  Do not add A, C,
and P into an opaque single score.

| Event type | Per-scale evidence |
| --- | --- |
| `transition` | C is high and A is at least moderately high; the curve changes state from left to right. |
| `focal_pulse` | A and P are high while C is not dominant; the curve has a localized peak or trough. |
| `mixed` | A, C, and P are all high at the same local event. |
| `broad_shift` | A exceeds the 90th percentile over a contiguous run of at least `max(3 bins, s / (2 * 4 kb))` without a localized C or P maximum.  Report an interval, not an artificial point. |
| `not_evaluable` | Availability or flank-quality requirements fail. |

## Cross-scale integration

Treat the ten IS windows as correlated views of the same data, not independent
replicates.  Group them into three bands:

- small: 32, 60, 92 kb;
- medium: 120, 152, 180 kb;
- large: 212, 240, 272, 300 kb.

Within a merged candidate region, cluster per-scale events of the same or
compatible geometry whose coordinates are within 12 kb (three 4 kb bins).  For
`broad_shift`, cluster overlapping intervals rather than summits.

The primary list requires support from at least two scale bands.  Events with
strong evidence in one band only remain in a separately labelled
`single_band_preliminary` result.  A scale band counts once regardless of how many
of its closely related windows support the event.

For point-like events (`transition`, `focal_pulse`, `mixed`), report a `summit`
without adding the three properties into an opaque score: first choose the
coordinate with the most supporting scale bands; then break ties by the median
relevant percentile (C for `transition`, P for `focal_pulse`, and the minimum of
A/C/P for `mixed`); then use the central genomic coordinate.  Report the cluster
extent as the positional uncertainty interval.  For a broad event, report its
interval and an optional summit only as a visual aid.  If a candidate region
contains multiple spatially separate event clusters, retain each locus separately
and label the region `complex` at the region-summary level.

## Outputs

`candidate_regions.tsv` contains one row per merged candidate region with the
upstream KNN/Meltron provenance.

`is_change_loci.tsv` contains one row per high-confidence event with:

```text
comparison, region_id, chrom, event_start, event_end, summit,
is_change_type, confidence_tier, support_scale_bands, support_scales,
max_amplitude_percentile, max_transition_percentile, max_pulse_percentile,
broad_shift_length_bp, mask_status, coverage_state_summary, rank_within_region
```

`single_band_preliminary_is_change_loci.tsv` contains one-band events with the
same schema and a lower-confidence label.

`masked_or_not_evaluable_intervals.bed` records candidate positions that cannot
be safely assessed from IS because of missing-bin rules.

Diagnostic plots show A and B IS tracks and delta tracks for all ten scales,
highlight the candidate-region core, mask unavailable positions, and overlay the
reported event coordinates or broad intervals.

## Validation criteria

- Interval merging tests cover overlap, adjacency, positive gaps, and chromosome
  boundaries.
- Masking tests prove that `0` remains valid, `NA` is unavailable, coordinates
  are never compressed, and short missingness does not invalidate an otherwise
  well-covered flank.
- Synthetic IS tracks verify that known step, local pulse, broad-shift, mixed,
  and missing-data patterns receive the intended classifications.
- A manual review set overlays the highest-ranked calls on multiscale IS plots;
  this is descriptive validation, not replicate-based inference.

## Out of scope

- Calling or labelling TAD boundaries in this stage.
- Declaring biological P values or population-level significance.
- Using empty-bin state as numeric IS-change evidence.
- Inferring the molecular class of a structural rearrangement from IS alone.
