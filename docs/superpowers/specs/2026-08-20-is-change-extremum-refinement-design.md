# IS-change local-extremum refinement and recovery

## Goal

Improve the genomic coordinate of existing IS-change calls and recover strong,
previously uncalled IS-change extrema inside already selected candidate regions.
The KNN/Meltron candidate selection and the existing primary/preliminary
threshold logic remain the high-precision first stage.

This design addresses two separate failures:

1. An existing point event can fall beside the biologically informative local
   peak or valley of the multi-scale Delta-IS curve.
2. A strong local peak or valley inside a candidate core can fail the current
   amplitude/transition/pulse combination and therefore never become an event.

No new KNN/Meltron candidate windows are created. The 300 kb plotting context
is display-only and is never searched for new events.

## Terminology

- **Delta-IS**: `IS_B - IS_A` at one 4 kb bin and one IS scale.
- **Scale z-score**: the existing robustly standardized signed Delta-IS from
  `score_scale`.
- **Consensus curve**: the median of all finite scale z-scores at each bin.
- **Candidate core**: `[region_start, region_end)`; the only allowed search
  domain for refinement and recovery.
- **Continuous valid segment**: adjacent valid 4 kb bins within one candidate
  core. Empty bins and invalid IS values terminate a segment and cannot be
  crossed.
- **Plateau**: consecutive equal-valued bins constituting one local peak or
  valley.

## Stage 1: retain the existing caller

The pipeline continues to:

1. select KNN/Meltron candidate regions;
2. compute the current per-scale amplitude, transition, and pulse scores;
3. classify and cluster point/broad events; and
4. assign the current primary or single-band preliminary confidence tier.

This stage is not relaxed or re-thresholded.

## Stage 2: refine existing point-event coordinates

Only non-broad events are eligible. For each such event:

1. Build the consensus curve over its candidate core from the existing scale
   z-scores.
2. Determine target polarity from the signed consensus at the original summit:
   positive seeks a peak; negative seeks a valley. A non-finite or zero
   consensus is ambiguous and remains unchanged.
3. Scan left and right in the same continuous valid segment. For each side,
   identify the first local extremum of the target polarity. There is no
   arbitrary genomic search radius.
4. Choose the closer of the two extrema. On an equal distance, choose the one
   with larger absolute consensus Delta-IS. For a plateau, choose the plateau
   bin closest to the original summit.
5. Accept the move only if at least two configured scale bands have a signed
   z-score at the proposed bin consistent with the consensus polarity.

On acceptance, point-event `summit`, `event_start`, and `event_end` are set to
the refined 4 kb bin. Broad shifts retain their original intervals. If two
point events converge on one refined plateau, only the stronger multi-scale
evidence is retained in the final event table.

## Stage 3: recover residual candidate-core extrema

Run after Stage 2. It is a recall-oriented supplement, not a replacement for
the existing caller.

1. Find plateau-aware local peaks and valleys of the consensus curve in every
   continuous valid segment of every candidate core.
2. Exclude extrema already occupied by a refined or unchanged original event.
3. Calculate local prominence without a fixed genomic radius. A peak is
   compared with its nearest flanking valleys and a valley with its nearest
   flanking peaks; the limiting flank defines its prominence.
4. Pool residual extrema across the comparison's candidate cores. A residual
   extremum is eligible only when its prominence is at or above the existing
   `high_percentile` (currently P95) of that pooled distribution.
5. Require at least two configured scale bands whose signed z-scores agree
   with the consensus polarity and whose existing amplitude percentile is at
   least `high_percentile`.

Each eligible extremum becomes a new `primary` event with type
`recovered_peak` or `recovered_valley`. Existing calls are not demoted or
reclassified. Multiple residual rows mapped to one plateau are collapsed to
the strongest multi-scale evidence row.

## Output and provenance

Primary and preliminary TSVs retain their current fields and add compatible
audit columns:

- `discovery_mode`: `original`, `refined`, or `recovered`;
- `original_summit`: the caller's summit before refinement, missing for a
  recovered event;
- `refinement_status`: `refined`, `unchanged`, `ambiguous`, or `recovered`;
- `refinement_direction`: `peak`, `valley`, or missing where not applicable;
- `refinement_support_scale_bands`: bands supporting the final coordinate.

The diagnostic plotter continues to draw all non-broad events as red point
markers, so recovered events appear at their final peak/valley coordinate.
No change is made to broad-shift rendering, candidate-region files, or
KNN/Meltron output.

## Failure handling

- Never cross an invalid/empty bin or candidate-core boundary.
- Preserve the original point event when polarity is ambiguous, no qualifying
  extremum exists, or final multi-band directional support fails.
- Do not emit recovered events for invalid bins, candidate context outside the
  core, or extrema already occupied by an original event.
- Keep all final coordinates aligned to the 4 kb resolution.

## Validation

Add focused synthetic tests for:

1. peak and valley refinement to the nearest valid extremum;
2. choosing the nearest extremum rather than a farther, stronger extremum;
3. plateau selection, invalid-bin stopping, ambiguous polarity, and failed
   multi-band support preserving the original event;
4. point-event collision collapse and broad-shift immutability;
5. residual peak/valley recovery with two-band directional and amplitude
   support;
6. rejection of weak, single-band, occupied, invalid, and context-only
   extrema; and
7. output-schema provenance and final plot-marker coordinates.

Run the full IS-change suite and the normal workflow through Donau. The normal
result figures are the final visual review artifact; no region-specific manual
replay is required by this specification.

## Acceptance criteria

- Existing primary calls retain their original candidate provenance and
  confidence tier, but point coordinates move to the nearest supported
  consensus peak or valley when appropriate.
- Strong, unoccupied, multi-band local extrema inside candidate cores can be
  reported as recovered primary events without lowering the original caller's
  thresholds.
- No event is created outside a candidate core or through empty-bin gaps.
- Each output row exposes whether it was original, refined, or recovered, and
  every displayed point marker matches its final 4 kb coordinate.
