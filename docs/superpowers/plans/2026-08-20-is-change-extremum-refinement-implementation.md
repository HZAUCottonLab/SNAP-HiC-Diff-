# IS-change Extremum Refinement Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (- [ ]) syntax for tracking.

**Goal:** Within existing KNN-plus-Meltron candidate cores, move initially called IS-change point events to their nearest supported multi-scale Delta-IS peak or valley, then recover strong, unoccupied local extrema missed by the existing amplitude/transition/pulse conjunction.

**Architecture:** Preserve the current per-scale scorer, classifier, clusterer, and primary/preliminary assignment as Stage 1. Retain each candidate core's existing per-scale score tracks in memory, derive one conservative signed consensus Delta-IS z-score curve, then run: (1) coordinate refinement for existing non-broad clusters and (2) comparison-wide percentile recovery for residual extrema. The caller serializes provenance in its existing tables; the plotter continues to mark final summit coordinates.

**Tech Stack:** Python 3.12, NumPy, pandas, cooler, pytest, Snakemake, Donau dsub, existing SNAP-HiC-Diff Conda environment.

## Global Constraints

- Do not change KNN empirical-P-value selection, Meltron filtering, candidate merging, candidate-core coordinates, IS values, IS scales, score calculations, current thresholds, or broad-shift semantics.
- Search only the half-open candidate core [region_start, region_end). The 300 kb context supports Stage-1 score calculation and plotting only; it is never searched.
- A searchable consensus segment is a run of adjacent 4 kb bins where every configured scale has valid, finite z. Empty balanced bins, missing IS values, non-finite z, coordinate discontinuities, and candidate-core boundaries terminate a segment. Never bridge, interpolate, or cross a break.
- Consensus z at a searchable bin is the median of signed per-scale robust z-scores, where Delta-IS = IS_B - IS_A. Positive values seek peaks, negative values seek valleys, and zero/missing values are ambiguous.
- Stage 2 applies to existing point events in both primary and preliminary tiers. Broad shifts retain their current interval and are not coordinate-refined.
- Stage 3 may create only primary recovered_peak/recovered_valley rows. It must not create candidate windows, demote original rows, or relax Stage-1 calling criteria.
- Require settings.min_primary_scale_bands distinct directional bands (currently two) and reuse settings.high_percentile for both residual-prominence and per-scale amplitude gates. Add no new threshold.
- Preserve old output columns, appending audit columns in this exact order: discovery_mode, original_summit, refinement_status, refinement_direction, refinement_support_scale_bands.
- All final starts, ends, and summits remain 4 kb aligned. Moved and recovered point events use exactly [summit, summit + resolution_bp).
- All tests, caller runs, workflow refreshes, and real-data checks run in Donau dsub jobs, never on the login node.
- Do not commit results, PNGs, cluster logs, temporary files, or __pycache__.

---

## File structure

- Modify: scripts/is_change_detection.py — pure consensus, plateau/extremum, support, refinement, collision, and recovery helpers.
- Modify: scripts/call_is_change_loci.py — retain per-core tracks, orchestrate both stages, serialize audit fields.
- Modify: test/is_change/test_is_change_detection.py — numerical geometry and recovery unit tests.
- Modify: test/is_change/test_call_is_change_loci.py — output schemas and caller aggregation tests.
- Modify: test/is_change/test_plot_is_change_loci.py — provenance-compatible marker tests.
- Modify: test/is_change/run_is_change_realdata.dsub — new real-data invariants.
- Create: test/is_change/run_is_change_refresh.dsub — post-merge normal-results refresh on a compute node.

## Shared contracts

### Candidate-core score tracks

Change the private region analyser to return its already calculated per-scale core tracks:

~~~
clusters, masks, core_score_tracks = _event_rows_for_region(\n    region, settings, track_a, track_b, bins_a, bins_b,\n)

# Key: configured scale_bp. Value: score_scale output, clipped exactly to
# [region_start, region_end), retaining chrom, coordinate, valid, z, and
# evidence-percentile columns.
core_score_tracks: dict[int, pd.DataFrame]
~~~

These data stay in memory only; no raw-score output file is added. Update caller
test doubles from (clusters, masks) to (clusters, masks, {}).

### Detection-module helpers

Add these pure interfaces to scripts/is_change_detection.py:

~~~
CONSENSUS_COLUMNS = [
    "chrom", "coordinate", "consensus_z", "valid", "segment_id",
]
EXTREMUM_COLUMNS = [
    "chrom", "coordinate", "plateau_start", "plateau_end", "direction",
    "segment_id", "prominence", "left_flank_coordinate",
    "right_flank_coordinate",
]

def build_consensus_track(
    score_tracks: dict[int, pd.DataFrame], settings
) -> pd.DataFrame

def local_consensus_extrema(
    consensus: pd.DataFrame, settings
) -> pd.DataFrame

def directional_support(
    score_tracks: dict[int, pd.DataFrame], coordinate: int, direction: str,
    settings, require_high_amplitude: bool = False,
) -> tuple[list[str], list[int]]
~~~

build_consensus_track validates that all configured scales exist and share the
same chromosome and 4 kb coordinate grid. A row is valid only when every scale
is valid and its z is finite. It computes the signed median only on valid rows.
It labels every maximal run of valid, grid-adjacent bins with a segment id;
invalid rows carry NaN consensus and segment id.

local_consensus_extrema groups adjacent numerically equal values as one plateau
(use a documented strict tolerance: np.isclose with rtol=0.0, atol=1e-12). A
peak plateau is strictly above both immediately adjacent plateaus in its valid
segment; a valley is strictly below both. Segment-edge plateaus are not
extrema. For a peak, prominence is:

~~~
min(peak_z - nearest_left_valley_z, peak_z - nearest_right_valley_z)
~~~

For a valley the analogous nearest flanking peaks are used. A missing flank
makes prominence NaN and therefore ineligible for recovery. The emitted
coordinate is the plateau's left edge; refinement selects the member nearest
an original summit.

directional_support returns configuration-ordered bands and ascending scales
that are valid, finite, and sign-consistent at one exact coordinate. Peak
requires z > 0 and valley requires z < 0. With require_high_amplitude=True it
also requires amplitude_percentile >= settings.high_percentile. Distinct bands,
not scale count, determine acceptance.

### Provenance and priority

In scripts/call_is_change_loci.py define:

~~~
REFINEMENT_COLUMNS = [
    "discovery_mode", "original_summit", "refinement_status",
    "refinement_direction", "refinement_support_scale_bands",
]
~~~

Append this list to PRIMARY_COLUMNS. Use these values:

| Outcome | discovery_mode | original_summit | refinement_status | refinement_direction |
|---|---|---:|---|---|
| Original broad shift | original | missing | unchanged | missing |
| Clear but unchanged point | original | initial summit | unchanged | peak/valley |
| Ambiguous point | original | initial summit | ambiguous | missing |
| Moved point | refined | initial summit | refined | peak/valley |
| Recovered residual | recovered | missing | recovered | peak/valley |

refinement_support_scale_bands records comma-separated directional bands at the
final coordinate without the recovery amplitude condition; it is empty where
no direction exists. It is separate from current Stage-1 support_scale_bands.

Use this deterministic collision key whenever two point rows share a final
plateau: descending distinct Stage-1 support-band count; descending maximum
transition, pulse, then amplitude percentile; ascending original summit
(or final summit when original is missing). This matches current caller ranking
precedence and gives a reproducible winner.

### Two-pass caller orchestration

output_rows first analyses all candidate regions and retains each region,
Stage-1 clusters, masks, and core score tracks. It then:

1. makes a consensus and refines existing point events independently per core;
2. retains all initial and final point coordinates as recovery occupancy,
   including a Stage-1 row later removed by refinement collision;
3. pools residual extrema from every candidate core in this comparison,
   thresholds their prominence once, and creates recovery rows; and
4. appends recovery rows to the correct region's primary list, then performs
   the existing caller sort/rank procedure unchanged.

The global pool is essential: percentile thresholding separately per candidate
core is not the approved criterion.

## Task 1: Build consensus, plateaus, prominence, and directional support

**Files:**
- Modify: scripts/is_change_detection.py
- Modify: test/is_change/test_is_change_detection.py

**Interfaces:** Implement the three pure helpers in the shared contract.
They consume only score_scale frames and settings; they do not read files,
coolers, classifiers, or caller state.

- [ ] **Step 1: Write failing geometry tests**

Add a helper making aligned score frames for real configured small, medium, and
large scales. Add these tests before implementing the helper functions:

~~~
def test_consensus_median_uses_signed_z_and_common_valid_segments():
    # median([-3, -1, 1]) == -1; a false/NaN row creates a segment break.

def test_plateau_aware_extrema_emit_one_peak_and_valley_with_prominence():
    # valley -> two-bin peak plateau -> valley; check plateau bounds and
    # limiting-flank prominence.

def test_segment_edges_and_invalid_gaps_have_no_recovery_prominence():
    # neither condition can supply both opposite-polarity flanks.

def test_directional_support_requires_distinct_bands_and_high_amplitude_when_requested():
    # two small-scale rows fail; small+medium succeeds; low amplitude rejects
    # only the recovery-style query.
~~~

Use 4 kb coordinates and assert that no extrema appears in an invalid bin.

- [ ] **Step 2: Verify RED on a compute node**

~~~
dsub -s test/is_change/run_is_change_tests.dsub
~~~

Inspect djob -ll <jobid> and logs/cluster/is_change_tests_<jobid>.out.
Expected: the added tests fail because the helpers do not exist. Do not run
pytest locally.

- [ ] **Step 3: Implement the pure helpers**

Order tracks by settings.scales_bp; reject missing scales or mismatched
chromosome/coordinate series with ValueError. Construct common validity from
valid & finite z for every scale, calculate median z only on common-valid
rows, and scan actual coordinates to assign segments.

Extract one plateau for every contiguous equal-valued run within a segment.
Compare only immediate plateau neighbours and fill opposite-polarity flanks
after all extrema are known. Keep extrema with NaN prominence for refinement
logic, but reject them in recovery. directional_support must return empty
lists, not a neighbour value or an exception, for a missing coordinate.

- [ ] **Step 4: Verify GREEN and commit**

~~~
dsub -s test/is_change/run_is_change_tests.dsub
~~~

Expected: existing tests plus the four new tests pass. Then:

~~~
git add scripts/is_change_detection.py test/is_change/test_is_change_detection.py
git commit -m "feat: add multiscale IS consensus extrema"
~~~

## Task 2: Refine existing point summits and serialize provenance

**Files:**
- Modify: scripts/is_change_detection.py
- Modify: scripts/call_is_change_loci.py
- Modify: test/is_change/test_is_change_detection.py
- Modify: test/is_change/test_call_is_change_loci.py

**Interfaces:**

~~~
def refine_point_clusters(
    clusters: pd.DataFrame, consensus: pd.DataFrame, extrema: pd.DataFrame,
    score_tracks: dict[int, pd.DataFrame], settings,
) -> pd.DataFrame
~~~

The result retains all current cluster fields plus audit fields. It can contain
fewer point rows after collision resolution. Broad rows are retained unchanged.

- [ ] **Step 1: Write failing refinement and schema tests**

Add these detection tests:

~~~
def test_refinement_moves_a_positive_event_to_nearest_supported_peak():
    # A farther peak has larger absolute z; the nearer peak must win.

def test_refinement_moves_a_negative_event_to_nearest_supported_valley():
    # Assert exactly one final 4 kb bin.

def test_refinement_chooses_nearest_plateau_member_and_stops_at_invalid_gap():
    # No arbitrary radius and no gap crossing.

def test_refinement_preserves_ambiguous_or_single_band_event_and_broad_shift():
    # zero z -> ambiguous; one-band candidate -> unchanged; broad immutable.

def test_refinement_collision_keeps_stronger_multiband_event():
    # Two calls converge on one plateau; only evidence-key winner remains.
~~~

In the caller tests, append REFINEMENT_COLUMNS to expected primary and
preliminary schema, update the patched region analyser to return an empty
track mapping, and add a caller integration test asserting:

- old columns retain their exact order and the five fields are appended;
- a refined primary reports refined status, its original summit, and at least
  two refinement-support bands;
- an unchanged preliminary row remains preliminary and does not lose its
  existing Stage-1 tier/support fields.

- [ ] **Step 2: Verify RED on a compute node**

~~~
dsub -s test/is_change/run_is_change_tests.dsub
~~~

Expected: new refinement imports/schema assertions fail before implementation.

- [ ] **Step 3: Implement Stage 2 while leaving Stage 1 intact**

Add a pure selector:

~~~
def select_nearest_extremum(
    original_summit, direction, extrema, consensus, resolution_bp
)
~~~

For each non-broad cluster:

1. save its initial summit as original_summit;
2. read consensus at that coordinate; positive means peak, negative valley,
   otherwise mark original/ambiguous without moving;
3. consider same-direction extrema in the same segment, find the first
   eligible plateau to the left and first to the right;
4. convert each plateau to its member closest to original_summit, choose lower
   distance, then larger absolute consensus and lower coordinate on a tie;
5. call directional_support at that member without the amplitude condition and
   move only with the configured count of distinct bands.

On move, set summit/event_start/event_end to exactly the chosen bin. On a clear
but rejected search, retain the complete old Stage-1 geometry and set unchanged.
For clear rows record directional support at their final coordinate; for
ambiguous rows and broad shifts write empty support. Preserve all Stage-1
confidence, event type, evidence percentile, and band/scale fields.

After refinement, collapse only converged point rows using the shared
deterministic key. In the caller, retain each scale's core score table during
the existing mask loop, return it as the third value, and apply refinement
after all region analyses have been collected. Mask production and final sort/
rank behaviour otherwise remain unchanged.

- [ ] **Step 4: Verify GREEN and commit**

~~~
dsub -s test/is_change/run_is_change_tests.dsub
~~~

Expected: legacy tests and all new refinement/schema tests pass.

~~~
git add scripts/is_change_detection.py scripts/call_is_change_loci.py \
  test/is_change/test_is_change_detection.py test/is_change/test_call_is_change_loci.py
git commit -m "feat: refine IS-change point summits"
~~~

## Task 3: Recover residual extrema with comparison-wide prominence

**Files:**
- Modify: scripts/is_change_detection.py
- Modify: scripts/call_is_change_loci.py
- Modify: test/is_change/test_is_change_detection.py
- Modify: test/is_change/test_call_is_change_loci.py

**Interfaces:**

~~~
def residual_extrema(
    consensus: pd.DataFrame, occupied_coordinates: set[int], settings
) -> pd.DataFrame

def recover_residual_events(
    residuals_by_region, score_tracks_by_region, settings
) -> dict[object, pd.DataFrame]
~~~

Returned events use existing cluster fields plus audit fields. All are primary
and use recovered_peak or recovered_valley.

- [ ] **Step 1: Write failing recovery tests**

Create score tracks with at least three extrema and add:

~~~
def test_recovery_adds_supported_peak_and_valley_at_global_prominence_cutoff():
    # pool two cores and use a deterministic P95 cutoff.

def test_recovery_rejects_single_band_low_amplitude_and_low_prominence_extrema():
    # all three independent gates must hold.

def test_recovery_excludes_initial_final_occupied_plateaus_and_invalid_bins():
    # occupancy includes pre- and post-refinement positions.

def test_recovery_never_uses_context_or_segment_edge_extrema():
    # core score tracks only; no context search and no one-sided feature.

def test_recovery_collapses_duplicate_plateau_rows_by_evidence():
    # overlapping ownership/duplicate inputs retain one deterministic winner.
~~~

Add a subprocess caller test using at least two candidate regions. It must show
an eligible residual only in primary, with recovered provenance, a missing
original summit, a one-bin interval, a sign-matching type, and two or more
refinement-support bands. Assert no recovered row exists in preliminary.

- [ ] **Step 2: Verify RED on a compute node**

~~~
dsub -s test/is_change/run_is_change_tests.dsub
~~~

Expected: imports and recovery output assertions fail before residual logic
exists.

- [ ] **Step 3: Implement residual pooling and recovery**

residual_extrema starts from local_consensus_extrema. Exclude every plateau
containing an occupied coordinate. The occupancy union includes all Stage-1
point summits and final Stage-2 summits, even when a collision later drops a
row. Retain only finite, non-negative prominence; never derive occupancy from
the 300 kb context.

In output_rows, concatenate residual extrema from all candidate cores. When the
pool is empty, return no recovery rows. Otherwise calculate once:

~~~
prominence_cutoff = residual_pool["prominence"].quantile(settings.high_percentile)
~~~

Keep prominence >= cutoff. For a remaining plateau choose the member of largest
absolute consensus (then lower coordinate) and call directional_support with
require_high_amplitude=True. Require the configured number of distinct bands.
Construct each eligible result as:

~~~
{
    "event_type": "recovered_peak" if direction == "peak" else "recovered_valley",
    "event_start": coordinate,
    "event_end": coordinate + settings.resolution_bp,
    "summit": coordinate,
    "confidence_tier": "primary",
    "support_scale_bands": ",".join(bands),
    "support_scales": ",".join(map(str, scales)),
    "discovery_mode": "recovered",
    "original_summit": np.nan,
    "refinement_status": "recovered",
    "refinement_direction": direction,
    "refinement_support_scale_bands": ",".join(bands),
}
~~~

For current evidence-percentile fields, use finite maxima at this coordinate
across directional support scales, following existing _finite_max semantics.
At serialization broad_shift_length_bp is zero. Collapse duplicate recovery
plateaus via the shared evidence key before appending them to their region's
primary list. Do not rerun cluster_multiscale_events, alter original
preliminary rows, or recalculate Stage-1 calls.

- [ ] **Step 4: Verify GREEN and commit**

~~~
dsub -s test/is_change/run_is_change_tests.dsub
~~~

Expected: all recovery acceptance/rejection cases and the full existing suite
pass.

~~~
git add scripts/is_change_detection.py scripts/call_is_change_loci.py \
  test/is_change/test_is_change_detection.py test/is_change/test_call_is_change_loci.py
git commit -m "feat: recover residual IS-change extrema"
~~~

## Task 4: Preserve plotting and gather compute-node evidence

**Files:**
- Modify: test/is_change/test_plot_is_change_loci.py
- Modify: test/is_change/run_is_change_realdata.dsub
- Create: test/is_change/run_is_change_refresh.dsub

**Interfaces:** The plotter receives primary rows with additional fields but
continues to draw every non-broad type as a red marker at final summit. The
existing temporary real-data runner validates data only; the new refresh
runner is an explicit post-merge normal-result regeneration step.

- [ ] **Step 1: Add plot provenance/coordinate tests**

Extend a primary fixture with all five audit fields and a recovered_valley whose
final summit differs from its illustrative original summit. Assert:

1. the plotter accepts appended audit fields without relaxing its old required
   columns;
2. heatmap axes and all IS axes draw the red marker at final summit; and
3. no marker is drawn at original_summit.

Do not change YlOrRd/log1p heatmap behaviour, shared x axes, broad shading, or
top-N diagnostic selection.

- [ ] **Step 2: Strengthen the temporary real-data replay**

After caller completion in test/is_change/run_is_change_realdata.dsub, append
an inline Python check that reads temporary primary/preliminary tables and
prints:

~~~
primary_discovery_mode_counts=<counts>
preliminary_discovery_mode_counts=<counts>
refined_point_alignment=<true/false>
recovered_point_alignment=<true/false>
~~~

Fail if either table lacks audit columns; any non-broad final point is not
4 kb-aligned with a one-bin interval; or any recovered row is not primary.
Do not assert biological locus counts or add A05/region-56/region-911-specific
filters.

- [ ] **Step 3: Create a post-merge figure refresh runner**

Create test/is_change/run_is_change_refresh.dsub:

~~~
#!/bin/bash
#DSUB -n is_change_refresh
#DSUB -q arm
#DSUB -aa
#DSUB -R 'cpu=2'
#DSUB -o logs/cluster/is_change_refresh_%J.out
set -euo pipefail
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
ROOT=/share/home/xhhuang24/project/PanTE/SNAP-HiC-Diff
cd "$ROOT"
snakemake -s Snakefile --use-conda --conda-prefix "$ROOT/.snakemake/conda" \
  --cores 2 --rerun-incomplete \
  --forcerun call_is_change_loci plot_is_change_loci \
  results/is_change/H1398_vs_L134/H1398_vs_L134_is_change_loci.tsv \
  results/figures/is_change/H1398_vs_L134/manifest.tsv
~~~

Make it executable and run bash -n only before merge. Do not submit it from the
feature worktree because its hard-coded root is the merged main worktree.

- [ ] **Step 4: Run GREEN suite and temporary H1398/L134 replay on Donau**

~~~
dsub -s test/is_change/run_is_change_tests.dsub
dsub -s test/is_change/run_is_change_realdata.dsub
~~~

Inspect djob -ll <test-jobid>, djob -ll <replay-jobid>, and the associated
logs. Expected: all tests pass; the existing KNN/Meltron candidate evidence
remains intact; audit columns and both printed alignment invariants exist; and
there is no traceback, ValueError, or scheduler failure. Record MEM_MAX_SUM
and CPU use before increasing resources. Do not require a fixed recovered-locus
count.

- [ ] **Step 5: Commit validation assets**

~~~
git add test/is_change/test_plot_is_change_loci.py \
  test/is_change/run_is_change_realdata.dsub \
  test/is_change/run_is_change_refresh.dsub
git commit -m "test: validate IS-change extremum outputs"
~~~

## Task 5: Review, merge, and visual handoff

**Files:**
- Verify: all commits on codex/is-extremum-refinement
- Verify after a user-approved merge: normal H1398/L134 TSV and PNG manifest

- [ ] **Step 1: Prepare an independent review package**

~~~
BASE=$(git merge-base master HEAD)
bash /share/home/xhhuang24/.codex/skills/subagent-driven-development/scripts/review-package \
  "$BASE" HEAD .superpowers/sdd/is-extremum-refinement-review.diff
~~~

Request review of common-valid gaps, candidate-core bounds, sign direction,
plateau/tie handling, global P95 pooling, distinct-band support, collision
determinism, audit schema order, and final-summit plotting.

- [ ] **Step 2: Recheck branch and evidence**

~~~
git diff --check "$BASE"..HEAD
git status --short
djob -ll <test-jobid>
djob -ll <replay-jobid>
~~~

Expected: no whitespace errors, intended source/test/runner changes only,
successful compute jobs, and no generated artifacts staged.

- [ ] **Step 3: Merge only with user approval, then regenerate figures**

Report commits, reviewer outcome, and dsub evidence; do not merge or delete the
worktree without the user's choice. Once merged, submit:

~~~
dsub -s test/is_change/run_is_change_refresh.dsub
~~~

Inspect logs/cluster/is_change_refresh_<jobid>.out, the refreshed
results/is_change/H1398_vs_L134/H1398_vs_L134_is_change_loci.tsv, and
results/figures/is_change/H1398_vs_L134/manifest.tsv. The normal regenerated
figures, including the A05 candidate core containing 71.724 Mb, are the user's
visual review artifact. Do not add a special-case A05, region-56, or region-911
rule or manual replay filter.

