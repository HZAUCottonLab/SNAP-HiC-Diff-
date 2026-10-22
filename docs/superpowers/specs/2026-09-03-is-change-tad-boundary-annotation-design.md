# IS-change Locus and Multi-resolution TAD Boundary Annotation Design

## Status

Approved on 2026-09-03.

## Goal

Add a downstream annotation layer that determines whether each IS-change locus
is spatially consistent with a TAD-boundary event between two samples. The
annotation operates only at the boundary level. It does not compare complete
TAD intervals or infer split, fusion, or nested-domain topology.

For the initial comparison, sample A is H1398 and sample B is L134. Both are
mapped to the same reference genome and use the same 4 kb IS coordinate grid.
All directional labels are defined from A to B:

```text
gained_in_B       = boundary present in L134 but absent in H1398
lost_in_B         = boundary present in H1398 but absent in L134
strengthened_in_B = conserved boundary is stronger in L134
weakened_in_B     = conserved boundary is weaker in L134
```

The design has two deliberately separate axes:

1. `position_state` describes whether the boundary position is conserved,
   gained, lost, shifted, complex, ambiguous, or absent near the IS locus.
2. `strength_state` describes whether a positionally conserved boundary is
   strengthened, weakened, stable, or not consistently evaluable from the IS
   curves.

Consequently, a locus may be `position_state=conserved` and
`strength_state=weakened_in_B`. Positional conservation must not be interpreted
as stable insulation strength.

## Scientific Interpretation Boundary

The TAD calls and IS curves are derived from the same Hi-C experiments but use
different representations and algorithms. Agreement between them is
cross-method consistency evidence, not an independent experimental
replication. IS-derived strength must therefore be described as a mechanistic
subtype or supporting direction, not as independent validation of an
IS-selected locus.

There is one biological sample per condition. The workflow must not emit
population-level differential P values or claim that a locus is statistically
significant across biological replicates.

## Relationship to Existing Designs

This design consumes, but does not redefine, the outputs from:

- `2026-08-15-is-change-locus-design.md`, which defines IS-change event
  detection, scale bands, and the `primary` versus
  `single_band_preliminary` confidence tiers;
- `2026-08-21-3dnetmod-multires-tad-boundaries-design.md`, which defines
  4/8/20 kb 3DNetMod calls and per-signal multi-resolution boundary regions.

The main analysis annotates all `primary` IS-change loci. Preliminary loci are
annotated through the same algorithm and published separately as a sensitivity
output. They must not be pooled into the primary summary.

## Scope

In scope:

- combine raw and balanced multi-resolution boundary calls into one
  provenance-aware boundary collection per sample;
- construct one genome-wide boundary-event map per ordered sample comparison;
- distinguish conserved, gained, lost, shifted, complex, and ambiguous
  boundary-position states;
- calculate multi-scale IS valley strength at externally supplied 3DNetMod
  boundary positions;
- classify strength changes at conserved boundaries;
- assess whether IS strength changes have the expected direction for gained,
  lost, and shifted events;
- annotate primary and preliminary IS-change tables;
- publish threshold-sensitivity and stability diagnostics;
- optionally perform a set-level constrained permutation analysis.

Out of scope:

- full-TAD overlap, split, merge, fusion, or nested-domain reconstruction;
- replacing 3DNetMod with TADCompare, TADreg, DiffDomain, or another caller;
- re-calling IS-change loci from TAD boundaries;
- using IS strength to override a boundary-presence classification;
- locus-level inferential P values;
- population-level biological inference from the current one-sample-per-state
  comparison.

## Existing Input Contracts

### TAD boundaries

For each sample and signal, the existing workflow publishes:

```text
results/tad_boundaries/{sample}.multires.raw.boundaries.bed
results/tad_boundaries/{sample}.multires.balanced.boundaries.bed
```

Each row has seven fields:

```text
chrom  start  end  boundary_id  support_resolutions
n_supporting_resolutions  n_source_boundaries
```

The regions are BED half-open intervals. Within one signal, they are connected
components of actually overlapping 4, 8, and 20 kb endpoint intervals.
Adjacency alone is not merged. Transitive overlap can create a broad boundary
region and is retained as input evidence rather than silently split by this
stage.

The per-resolution boundary files remain available and are used only for
localization-jitter estimation and leave-one-resolution-out validation.

### IS curves

The current IS workflow reads the balanced 4 kb contact map and publishes ten
CWorld-normalized IS tracks at:

```text
32, 60, 92, 120, 152, 180, 212, 240, 272, and 300 kb
```

CWorld reports:

```text
IS = log2(smoothed_insulation / mean_smoothed_insulation)
```

Lower local IS minima indicate stronger local insulation. The existing scale
bands remain authoritative:

```text
small  = 32, 60, 92 kb
medium = 120, 152, 180 kb
large  = 212, 240, 272, 300 kb
```

### IS-change loci

The primary and preliminary tables already provide `chrom`, `event_start`,
`event_end`, `summit`, `is_change_type`, `confidence_tier`, scale support, and
coverage fields. Point-like events use the 4 kb summit bin as their annotation
anchor. A `broad_shift` uses its event interval.

## High-level Data Flow

```text
per-signal, multi-resolution TAD boundaries
  -> within-sample raw/balanced unification
  -> localization-jitter calibration
  -> ordered A/B boundary matching
  -> genome-wide boundary-event map
  -> IS valley-strength extraction at boundary coordinates
  -> primary and preliminary IS-locus annotation
  -> QC, sensitivity summaries, and optional set-level permutation
```

The boundary-event map is constructed before it is intersected with IS-change
loci. This prevents the local IS candidates from influencing which TAD
boundary events are called.

## Within-sample Raw/Balanced Unification

### Rationale

Raw and balanced boundary calls are two views of the same sample, not
biological replicates. The user-selected policy is a union: unmatched
boundaries are retained rather than discarded. Provenance must be preserved so
that `both`, `raw_only`, and `balanced_only` evidence can be distinguished.

### Matching rule

For each chromosome, sort the raw and balanced multi-resolution boundary
regions by coordinate. A raw/balanced pair is eligible only when their
half-open intervals actually overlap:

```text
max(raw_start, balanced_start) < min(raw_end, balanced_end)
```

Eligible one-to-one matches are selected with an order-preserving assignment
that maximizes the number of overlaps and then minimizes total center distance.
The union interval of a matched pair becomes one unified boundary with
`signal_support=both`. Unmatched inputs become `raw_only` or `balanced_only`
unified boundaries.

Support metadata are combined as follows:

```text
support_resolutions     = union of input resolution sets
n_supporting_resolutions = number of distinct resolutions
n_source_boundaries     = sum of original endpoint counts
```

### Complex components

If an overlap component contains more than one raw and more than one balanced
input, or one input has two equally close counterparts within one 4 kb source
bin, the component is not forced into arbitrary one-to-one simple boundaries.
It is emitted as one diagnostic union region with
`component_state=complex_within_sample` and the complete member-ID lists.

Complex component regions do not enter simple A/B matching. Any cross-sample
event or IS locus intersecting them is classified as `complex` unless later
evidence uniquely resolves the component through an explicitly reviewed
algorithm extension.

## Localization-jitter Calibration

### Purpose

A small separation between boundary centers may arise from binning,
normalization, or resolution rather than biological movement. The workflow
must estimate this technical localization scale before assigning `shifted`.

### Calibration observations

Pool high-confidence offsets from both samples using:

1. one-to-one, overlapping raw/balanced boundary pairs;
2. reciprocal one-to-one matches between the full multi-resolution set and
   each leave-one-resolution-out set;
3. reciprocal matches among stable 4, 8, and 20 kb boundary regions.

Exclude complex components, chromosome-end-truncated intervals, missing-data
regions, and pairs without a unique reciprocal closest match.

The primary jitter threshold is:

```text
J = max(4000 bp, empirical 95th percentile of eligible center offsets)
```

At least 100 pooled eligible pairs are required. If fewer are available, use a
20 kb fallback, record `jitter_source=fallback_max_tad_resolution`, and
downgrade shifted-event evidence by one tier.

If the estimated `J` is at least the primary maximum shift distance, the data
cannot distinguish technical jitter from a simple shifted event. In that case,
the workflow publishes the QC failure and reports otherwise eligible shifts as
`ambiguous`, rather than silently increasing the search radius.

## Cross-sample Boundary Matching

### Direction

All comparisons are ordered. For H1398 versus L134:

```text
sample_a = H1398
sample_b = L134
```

Swapping the comparison direction must swap gain/loss and
strengthened/weakened labels while leaving conserved and shifted geometry
unchanged.

### Stage 1: fixed-position matching

For each chromosome, simple unified boundaries are sorted by center. Eligible
overlapping A/B pairs are assigned one-to-one while preserving order.

- Overlapping matched uncertainty intervals are `conserved_exact`.
- Non-overlapping, mutually closest boundaries separated by no more than `J`
  are `conserved_within_jitter`.

Both subtypes publish `position_state=conserved`; the subtype is retained in
`position_substate`.

### Stage 2: local shifted matching

Only boundaries left unmatched by Stage 1 enter Stage 2. Use a banded,
order-preserving dynamic program with match, insertion, and deletion moves.
For eligible pair `(i, j)`:

```text
d_ij       = abs(center_a_i - center_b_j)
match_cost = d_ij / Dmax
gap_cost   = 0.5 per unmatched boundary
```

Matches with `d_ij > Dmax` are prohibited. Thus a pair at the maximum distance
costs the same as one A deletion plus one B insertion; closer pairs are
preferred, while the order constraint prevents crossing correspondences.

The primary maximum distance is:

```text
Dmax = 100000 bp
```

This is equivalent to five 20 kb bins, following the five-bin shifted-boundary
convention used by TADCompare while adapting it to the coarsest contributing
TAD resolution. The workflow must also rerun the matcher at 60 and 80 kb for
sensitivity reporting.

Stage-2 outcomes are:

```text
J < distance <= Dmax   -> shifted
unmatched A boundary   -> lost_in_B
unmatched B boundary   -> gained_in_B
```

### Match uniqueness and complex neighborhoods

For every Stage-2 pair, record its closest and second-closest eligible
counterpart distances. A match is geometrically ambiguous when:

- it is not a reciprocal closest match; or
- the closest and second-closest distances differ by less than one 4 kb bin.

A local eligible bipartite component with multiple A and multiple B
boundaries is `complex` when alternative non-crossing assignments satisfy the
same one-bin margin. The workflow must not convert such a region into a series
of apparently precise shifts.

### Position-state vocabulary

The final mutually exclusive states are:

```text
conserved
gained_in_B
lost_in_B
shifted
complex
ambiguous
```

`no_boundary_association` is added only at the later IS-locus annotation
layer; it is not itself a genome-wide boundary event.

## Boundary Event Representation

Each simple event stores separate A and B support intervals. There is no single
continuous support corridor for a shifted event.

- `conserved`: use each sample's matched uncertainty interval;
- `gained_in_B`: use the B boundary interval;
- `lost_in_B`: use the A boundary interval;
- `shifted`: store both the old A interval and the new B interval;
- `complex`: store the diagnostic component region;
- `ambiguous`: store every available interval plus the ambiguity reason.

For shifted events, `event_span_start` and `event_span_end` may describe the
full old-to-new span for visualization only. They must never be used as the
default interval-overlap target because doing so would label every position
between the old and new boundary as a boundary.

## IS-locus to Boundary-event Association

### Anchor intervals

```text
point-like IS event = [summit, summit + 4000)
broad_shift         = [event_start, event_end)
```

The TAD boundary inputs have already been expanded by their active TAD
resolution. No additional association padding is applied.

An IS anchor is associated when it overlaps at least one stored A or B support
interval of a boundary event. Merely falling between the old and new endpoints
of a shifted event is insufficient.

### One final event per IS locus

- zero associated events: `position_state=no_boundary_association`;
- one associated simple event: inherit that event's position state;
- multiple simple events, or one complex component: `position_state=complex`;
- one event with missing or conflicting evaluability: `position_state=ambiguous`.

The full candidate count and candidate IDs are retained. The workflow does not
silently choose a biologically convenient event type when multiple events are
present.

## IS-derived Boundary Strength

### Why local prominence is required

A single IS value mixes valley depth with the local baseline. Boundary strength
is therefore measured as the local prominence of a valley, consistent with the
cooltools practice of applying peak prominence to `-IS`.

For each sample, boundary interval, and IS scale `w`:

1. read the segment extending `w` bp on both sides of the interval;
2. apply the existing missing-value and flank-validity rules;
3. find plateau-aware local minima of IS within the boundary interval;
4. calculate their peak prominence on `-IS` within the local segment;
5. retain the largest eligible prominence as `S(sample, boundary, w)`.

If the interval and both flanks are sufficiently covered but contain no local
minimum, set strength to zero and record `strength_detection=no_valley`. If
coverage is insufficient, report `NA`; missing data must not be interpreted as
zero strength.

### Coordinate roles by event type

- conserved: evaluate A at the A interval and B at the B interval;
- gained in B: evaluate both IS tracks at the B boundary interval;
- lost in B: evaluate both tracks at the A boundary interval;
- shifted: evaluate both tracks at both the old A and new B intervals;
- complex/ambiguous: retain per-scale diagnostic values but do not issue a
  simple strength class.

### Robust effect standardization

At scale `w`, define:

```text
delta_strength_w = strength_b_w - strength_a_w
```

Use all genome-wide simple conserved boundary pairs, not only IS-change loci,
as the effect-size reference distribution. For each scale:

```text
robust_z_w =
  (delta_strength_w - median(delta_strength_w)) /
  (1.4826 * MAD(delta_strength_w))
```

This is an effect standardization, not a hypothesis test. At least 50
evaluable conserved pairs and a nonzero MAD are required at a scale. Otherwise
the scale is unevaluable for categorical strength calls.

For each small, medium, and large scale band, the band score is the median of
its valid `robust_z_w` values. A band requires at least two valid scales.

### Conserved-boundary strength classification

```text
at least 2 band scores >=  2 -> strengthened_in_B
at least 2 band scores <= -2 -> weakened_in_B
at least 2 evaluable bands and all abs(score) < 1 -> stable_strength
otherwise -> discordant_or_unevaluable
```

The output retains the number and names of supporting scale bands. Thresholds
1.5, 2.0, and 2.5 are evaluated in sensitivity analysis; 2.0 is primary.

### IS consistency for positional rearrangements

IS direction is recorded separately as `is_consistency_state`; it never
changes the boundary-derived `position_state`.

Expected patterns are:

```text
gained_in_B:
  strength_b - strength_a > 0 at the new B coordinate

lost_in_B:
  strength_b - strength_a < 0 at the old A coordinate

shifted:
  old A coordinate weakens in B, and
  new B coordinate strengthens in B
```

For gained and lost events, the expected direction must cross the primary
absolute robust-Z threshold in at least two scale bands to be
`direction_concordant`. For shifted events, the old and new endpoints are
judged with the same two-band threshold:

```text
both old-weakening and new-strengthening supported -> bidirectional_shift_support
only one endpoint supported                      -> partial_shift_support
neither supported or opposite direction          -> discordant
insufficient valid scales                         -> unevaluable
```

## Evidence Tier

The pipeline publishes all underlying evidence columns and uses the following
transparent summary tier. It does not calculate an opaque weighted score.

### High

- IS locus is `primary`;
- event is simple and uniquely matched;
- every present-side boundary is supported by both signals or at least two
  TAD resolutions;
- the position label is unchanged at 60, 80, and 100 kb matching radii;
- a rearranged event has direction-concordant IS support, or a conserved
  strength change is supported by at least two IS scale bands.

### Moderate

- IS locus is `primary` and the event is unique;
- exactly one of the high-tier support conditions is absent, for example a
  single-signal boundary or a shift that appears only at 80/100 kb.

### Low

- locus is `single_band_preliminary`; or
- a simple label depends on one signal, one resolution, or one IS scale band;
- jitter calibration used the fallback threshold.

Complex and ambiguous states retain `evidence_tier=not_applicable`. A clean,
fully evaluable `no_boundary_association` may receive a negative-evidence tier,
but it must not be presented as proof that no physical boundary exists.

## Evaluable Regions and Missing Data

An absence-based call (`gained_in_B`, `lost_in_B`, or clean
`no_boundary_association`) requires valid input coverage in the sample where a
boundary is declared absent. At minimum:

- the relevant chromosome and bins exist in both samples;
- the local 4 kb IS track passes the configured missing-value rules;
- every configured TAD resolution and signal completed successfully for the
  sample;
- the coordinate is not inside a declared bad or excluded region.

Failure of these conditions yields `ambiguous` or `unevaluable`, never a gain
or loss.

## Optional Set-level Spatial Association Test

The annotation itself requires no background distribution or per-locus P
value. An optional, explicitly exploratory permutation addresses a separate
set-level question:

> Are IS-change loci more frequently associated with TAD-boundary
> rearrangements than matched positions drawn from the same IS candidate
> search space?

The randomization universe is the union of original IS candidate regions after
applying the same evaluability mask. Each permutation preserves:

- the number of loci per chromosome;
- point versus broad-event geometry and interval length;
- primary or preliminary stratum;
- missingness/evaluability class;
- the original candidate-region restriction.

Run 10,000 permutations. Pre-specified statistics are:

```text
fraction of IS loci associated with any boundary event
fraction of boundary-associated loci that are rearranged rather than conserved
counts or fractions of gained, lost, and shifted events
median distance to the nearest eligible boundary
```

For an upper-tail statistic, report:

```text
empirical_p = (1 + number(null >= observed)) / (B + 1)
```

Also report observed/null fold enrichment and an empirical interval. Apply BH
correction only across the pre-specified collection of type-level statistics.
The interpretation is limited to spatial association in this specific H1398
versus L134 comparison.

## Configuration Contract

A proposed configuration section is:

```yaml
is_boundary_annotation:
  enabled: true
  association_padding_bp: 0
  raw_balanced_overlap_only: true
  base_bin_bp: 4000
  jitter_quantile: 0.95
  jitter_min_pairs: 100
  jitter_fallback_bp: 20000
  max_shift_bp: 100000
  shift_sensitivity_bp: [60000, 80000, 100000]
  match_tie_tolerance_bp: 4000
  strength:
    local_radius_scale_factor: 1.0
    min_reference_boundaries: 50
    min_valid_scales_per_band: 2
    robust_z_threshold: 2.0
    stable_abs_z_threshold: 1.0
    sensitivity_thresholds: [1.5, 2.0, 2.5]
  permutation:
    enabled: false
    n_permutations: 10000
    seed: 20260903
```

All distances are integer base pairs and exact multiples of the 4 kb source
grid. Sensitivity lists must be sorted, unique, and contain the primary value.
The maximum shift must exceed the jitter threshold before simple shift calls
are allowed.

## Output Contract

### Unified sample boundaries

```text
results/tad_boundaries/{sample}.multires.unified.boundaries.tsv
```

Required columns:

```text
chrom
start
end
unified_boundary_id
center_bp
component_state
signal_support
raw_boundary_ids
balanced_boundary_ids
support_resolutions
n_supporting_resolutions
n_source_boundaries
```

### Genome-wide comparison event map

```text
results/tad_boundary_events/{comparison}.boundary_events.tsv
```

Required columns include:

```text
comparison
sample_a
sample_b
chrom
boundary_event_id
position_state
position_substate
boundary_a_id
boundary_b_id
a_start
a_end
b_start
b_end
event_span_start
event_span_end
shift_bp
match_stage
match_is_reciprocal
second_best_margin_bp
n_local_candidates
jitter_threshold_bp
max_shift_bp
ambiguity_reason
```

### Annotated IS loci

```text
{MULTISCALE_ROOT}/{comparison}/is_change/
  {comparison}_is_change_loci.boundary_annotated.tsv

{MULTISCALE_ROOT}/{comparison}/is_change/
  {comparison}_is_change_loci_preliminary.boundary_annotated.tsv
```

Both retain every original IS-change column and add:

```text
boundary_event_id
position_state
position_substate
strength_state
is_consistency_state
evidence_tier
boundary_a_id
boundary_b_id
shift_bp
signal_support_a
signal_support_b
resolution_support_a
resolution_support_b
n_boundary_candidates
jitter_threshold_bp
max_shift_bp
annotation_reason
```

### Per-scale strength and QC

```text
results/tad_boundary_events/{comparison}.boundary_strength_by_scale.tsv
results/tad_boundary_events/{comparison}.annotation_qc.tsv
results/tad_boundary_events/{comparison}.annotation_summary.tsv
```

The long-format strength table contains one row per event, coordinate role,
and IS scale. The QC table records jitter calibration, counts at each filter
stage, missingness, complex components, label changes under sensitivity
settings, and any fallback behavior.

## Validation Strategy

### Deterministic unit cases

Tests must cover:

- exact overlap and BED adjacency semantics;
- raw/balanced one-to-one merge and unmatched retention;
- many-to-many components becoming complex;
- conserved-within-jitter versus shifted separation;
- one-to-one and non-crossing dynamic-programming matches;
- correct A-to-B gain/loss direction;
- failure to call absence when the opposite sample is unevaluable;
- no shifted-event corridor used for IS association;
- correct strength sign when an IS valley deepens or becomes shallow;
- zero strength for a valid interval without a valley versus `NA` for missing
  data.

### Direction-swap invariant

Running B versus A on the same inputs must satisfy exactly:

```text
gained_in_B <-> lost_in_B
strengthened_in_B <-> weakened_in_B
conserved -> conserved
shifted -> shifted
abs(shift_bp) unchanged
```

The signed shift changes sign. This invariant must pass for every simple event.

### Synthetic truth benchmark

Generate ordered boundary sets with known conserved, gained, lost, shifted,
and complex events. Sample technical offsets from the empirical jitter
distribution and test both sparse and crowded boundary neighborhoods.

Noiseless canonical cases must be classified exactly. Under empirical jitter,
simple events separated by more than `2 * Dmax` and injected shifts exceeding
`J + 4 kb` should achieve at least 0.90 macro-F1. Metrics must still be
reported when this target is missed; the implementation must not tune on the
evaluation set without recording the changed threshold.

### Internal stability analyses

Publish label concordance for:

- raw-only, balanced-only, and unified boundary inputs;
- leave-4-kb-out, leave-8-kb-out, and leave-20-kb-out boundary sets;
- 60, 80, and 100 kb maximum shift radii;
- strength robust-Z thresholds 1.5, 2.0, and 2.5.

Instability downgrades evidence rather than deleting records.

### Visual review

For each sufficiently populated simple class, review at least 20 randomly
selected primary loci; review all loci when a class has fewer than 20. Each
diagnostic view should align:

- the two Hi-C contact maps;
- raw and balanced boundary tracks;
- 4, 8, and 20 kb support;
- all ten IS curves;
- the assigned old/new boundary coordinates and final labels.

Complex and ambiguous examples must be deliberately sampled rather than
excluded from review.

## Error and Warning Semantics

Implementation errors should use stable leading codes:

```text
IS_BOUNDARY_CONFIG_INVALID
IS_BOUNDARY_INPUT_SCHEMA_INVALID
IS_BOUNDARY_REFERENCE_MISMATCH
IS_BOUNDARY_JITTER_CALIBRATION_FAILED
IS_BOUNDARY_MATCH_FAILED
IS_BOUNDARY_STRENGTH_FAILED
IS_BOUNDARY_ANNOTATION_FAILED
```

Recoverable scientific limitations, such as insufficient jitter pairs or an
unevaluable IS scale, are written to QC and reflected in the annotation state.
They should not abort unrelated chromosomes. Coordinate-system mismatches and
schema corruption are hard failures.

Outputs are written atomically and sorted deterministically by configured
chromosome order, coordinate, and stable event ID.

## Rejected Alternatives

### Greedy nearest-neighbor matching as the primary algorithm

Reciprocal nearest-neighbor matching is transparent and remains useful for
calibration and uniqueness checks. It is insufficient as the primary matcher
in crowded regions because an early local choice can force later boundaries
into incorrect gain/loss labels.

### Connected-component merging across samples

Pure graph-component merging can chain nearby boundaries into an increasingly
broad event. It is retained only for detecting and reporting complex local
geometry, not for assigning simple shifted events.

### Complete TAD alignment or differential-domain calling

Domain-based methods can identify split, merge, and topology changes, but the
existing multi-resolution 3DNetMod output contains overlapping hierarchical
TADs. Full-domain alignment would mix boundary evidence with an additional
domain-selection policy that is outside the requested scope.

### TADCompare as a drop-in replacement

TADCompare provides useful definitions for shifted and strength-change
boundaries, but it calculates its own boundary score from contact matrices.
Running it directly would replace rather than annotate the current 3DNetMod
boundary calls. This design adopts its interpretable event concepts while
keeping the established caller and multi-resolution evidence.

### Per-locus statistical significance

With one sample per condition and loci selected from the same IS curves later
used for strength annotation, per-locus P values would be circular and would
not estimate biological variance. The optional permutation is therefore
limited to set-level spatial association.

## Literature Basis

1. Cresswell, K. G., & Dozmorov, M. G. (2020). TADCompare: An R package for
   differential and temporal analysis of topologically associated domains.
   *Frontiers in Genetics, 11*, 158.
   https://doi.org/10.3389/fgene.2020.00158
   - Defines shifted boundaries using a configurable nearby-bin threshold and
     separates position shift from boundary-strength change.
   - Reports nonzero shifted and strength-change calls even between technical
     replicates, motivating explicit localization-jitter calibration.

2. Serra, F., Baù, D., Goodstadt, M., Castillo, D., Filion, G. J., &
   Marti-Renom, M. A. (2017). Automatic analysis and 3D-modelling of Hi-C data
   using TADbit reveals structural features of the fly chromatin colors.
   *PLoS Computational Biology, 13*(7), e1005665.
   https://doi.org/10.1371/journal.pcbi.1005665
   - TADbit documents both reciprocal-nearest and global
     Needleman-Wunsch-style alignment of ordered TAD-boundary lists.

3. Xu, J., Xu, X., Huang, D., et al. (2024). A comprehensive benchmarking with
   interpretation and operational guidance for the hierarchy of topologically
   associating domains. *Nature Communications, 15*, 4376.
   https://doi.org/10.1038/s41467-024-48593-7
   - Shows broadly that resolution materially changes TAD and boundary
     outputs, supporting multi-resolution evidence retention and
     leave-one-resolution stability checks. The paper catalogues 3DNetMod but
     excludes it from the final comparative benchmark because of parameter
     reproducibility concerns, so this evidence supports the validation
     strategy rather than directly validating the present 3DNetMod calls.

4. Gong, Y., Lazaris, C., Sakellaropoulos, T., et al. (2018). Stratification of
   TAD boundaries reveals preferential insulation of super-enhancers by strong
   boundaries. *Nature Communications, 9*, 542.
   https://doi.org/10.1038/s41467-018-03017-1
   - Demonstrates that shared boundary positions can vary in insulation
     strength and that strength carries biological information.

5. Open2C. cooltools insulation and boundary implementation.
   https://cooltools.readthedocs.io/en/stable/_modules/cooltools/api/insulation.html
   - Calculates boundary strength as peak prominence on minima of an
     insulation track, motivating the local `-IS` prominence statistic.

6. Heger, A., Webber, C., Goodson, M., Ponting, C. P., & Lunter, G. (2013).
   GAT: A simulation framework for testing the association of genomic
   intervals. *Bioinformatics, 29*(16), 2046-2048.
   https://doi.org/10.1093/bioinformatics/btt343
   - Motivates constrained genomic interval null models rather than uniform
     whole-genome placement.

7. Gel, B., Díez-Villanueva, A., Serra, E., Buschbeck, M., Peinado, M. A., &
   Malinverni, R. (2016). regioneR: An R/Bioconductor package for the
   association analysis of genomic regions based on permutation tests.
   *Bioinformatics, 32*(2), 289-291.
   https://doi.org/10.1093/bioinformatics/btv562
   - Supports chromosome-aware customizable randomization and evaluation of
     the local specificity of genomic associations.

## Final Design Decision

The approved implementation target is therefore:

```text
boundary-first genome-wide event map
+ raw/balanced provenance-preserving union
+ two-stage order-constrained one-to-one matching
+ data-calibrated jitter threshold
+ 100 kb primary shifted radius with 60/80/100 kb sensitivity
+ local IS-valley prominence across the existing three scale bands
+ separate position, strength, IS-consistency, and evidence-tier fields
+ descriptive per-locus annotation and optional set-level permutation only
```
