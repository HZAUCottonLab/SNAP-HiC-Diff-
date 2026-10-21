> **STATUS: SUPERSEDED — do not implement this plan.**
>
> The design below (containment-based split detection in `unify_boundaries`,
> per-locus boundary profiles with `profile_split_in_B` / `profile_merge_in_B`,
> density-region classification) was replaced before implementation. Two
> measurements killed it:
>
> 1. Boundary centres are derived quantities; the interval is the measurement.
>    Every rule built on centre distance produced artificial ties.
> 2. The tied components were reachable only through a transitive closure over
>    `max_shift_bp`, giving a median span of 1.39 Mb. At that scale "one
>    boundary became two" is not the question anyone asked.
>
> What was actually built is `scripts/is_boundary_point_state.py`: a locus is a
> point (its summit), a boundary is an interval, and the only question asked is
> whether the point falls inside one — per sample. See that module and
> `test/is_boundary/test_point_state.py` for the live behaviour.
>
> This file is kept for the record of what was considered and why it was
> rejected.

# Locus Boundary-Profile Classification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the boundary-event matching verdict at each IS-change locus with a direct per-sample boundary-profile classification, so that `complex` stops existing as an output state while every locus keeps a traceable list of contributing boundary coordinates.

**Architecture:** Three layers. (1) Unification stops calling a coarse/fine split a tie and emits a structured `split_within_sample_*` record instead. (2) A new pure-function module counts boundaries in a fixed window around each locus, per sample, and classifies the locus from the two counts and their offset — no cross-sample boundary pairing, therefore no tie. (3) A density rule absorbs dense regions into their own state. The existing event pipeline keeps running unchanged as a diagnostic view; the new classification is published in parallel columns first and only promoted to `position_state` after it agrees with the existing states on the loci where both should agree.

**Tech Stack:** Python 3.10+, pandas, NumPy, PyYAML, pytest, Snakemake 8+, existing `env/SNAP-HiC-Diff.yaml`. No new dependencies.

**Spec:** `docs/superpowers/specs/2026-09-03-is-change-tad-boundary-annotation-design.md` (the original boundary-annotation design; this plan amends its `complex` semantics, so read both).

## Global Constraints

- H1398 is sample A and L134 is sample B for `H1398_vs_L134`; all directional labels are A-to-B.
- Scale bands are fixed: small = 32/60/92 kb, medium = 120/152/180 kb, large = 212/240/272/300 kb.
- Base bin is 4,000 bp; TAD resolutions are 4,000 / 8,000 / 20,000 bp.
- A boundary interval is **not** a point. Never compare boundary centres as integers when deciding identity; compare intervals or sub-interval counts.
- `unevaluable_coverage_*` boundaries never contribute to a profile count.
- Every profile row must carry `chrom`, `start`, `end`, `center_bp`, and `boundary_id` so that a locus can be joined to an SV/CNV breakpoint table by coordinate.
- Outputs are deterministic, sorted, schema-checked, and atomically replaced.
- Missing coverage still produces `ambiguous`/`unevaluable`; it never produces gain or loss.
- Stable leading error codes are unchanged: `IS_BOUNDARY_CONFIG_INVALID`, `IS_BOUNDARY_INPUT_SCHEMA_INVALID`, `IS_BOUNDARY_MATCH_FAILED`, `IS_BOUNDARY_STRENGTH_FAILED`, `IS_BOUNDARY_ANNOTATION_FAILED`.
- Do not delete the existing `boundary_events.tsv` pipeline or its outputs; it remains the diagnostic view.
- Do not change `boundary_strength_by_scale.tsv` / `boundary_strength_summary.tsv` schemas in this plan.

## Review Focus

- **A locus whose window contains zero boundaries on both sides.** Expected: `profile_no_boundary`, distinct from `ambiguous`. A reasonable person would not call an empty window "complex".
- **A dense region with 40 boundaries in 100 kb.** Expected: exactly one `profile_dense_region` call with density numbers, not 40 pairwise judgements.
- **A coarse boundary that strictly contains two fine ones.** Expected: one `profile_split_in_B`-style call whose long-table rows list all three intervals with their own coordinates.
- **A locus whose anchor falls exactly on a boundary end coordinate.** Half-open interval semantics (`start <= x < end`) must be applied consistently in the window scan and the containment test, or the same boundary counts twice in one place and zero times in another.
- **A locus on a chromosome with an `unevaluable_coverage` boundary inside its window.** Expected: that boundary is invisible to the profile, and if it was the only one, the locus becomes `ambiguous`, never `lost_in_B`.

---

## File Structure

### Create

- `scripts/is_boundary_profile.py` — window construction, per-sample boundary collection, long-table row building, profile classification, and density classification. Pure functions only; no I/O.
- `test/is_boundary/test_profile.py` — unit tests for every classification branch, the half-open boundary semantics, and the unevaluable exclusion.

### Modify

- `scripts/is_boundary_io.py` — add `PROFILE_COLUMNS` and `PROFILE_SUMMARY_COLUMNS`.
- `scripts/is_boundary_unification.py` — replace the `_has_tie` short-circuit with the containment rule.
- `scripts/is_boundary_config.py` — add `ProfileSettings` to `ISBoundarySettings`.
- `scripts/is_boundary_annotation.py` — append profile columns to the annotated loci.
- `scripts/annotate_is_change_boundaries.py` — accept the two unified boundary tables and write the profile tables.
- `rules/is_boundary_annotation.smk` — pass boundaries into the annotate rule and declare the new outputs.
- `config/config.yaml` — add the `is_boundary_annotation.profile` block.
- `test/is_boundary/test_unification.py` — update the two tie tests to the containment semantics.
- `test/is_boundary/test_annotation.py` — extend the schema assertions for the new columns.

---

### Task 1: Containment replaces tie in boundary unification

**Files:**
- Modify: `scripts/is_boundary_unification.py:167-212` (`_has_tie`, `unify_boundaries` component loop)
- Test: `test/is_boundary/test_unification.py:89-132`

**Why:** 97% of measured ties are one wide interval containing two narrow ones, with the union's start and end each 4,000 bp inside the wide interval. That is a split, not an ambiguity. Encoding it as `complex_within_sample` throws away the structure and excludes the whole component downstream via `_is_excluded`.

- [ ] **Step 1: Write the failing test for a contained split**

Add to `test/is_boundary/test_unification.py`:

```python
def test_wide_boundary_containing_two_narrow_ones_is_a_split_not_a_tie():
    raw = boundary_frame("raw", [("A02", 0, 40000, "wide", "8000,20000", 2)])
    balanced = boundary_frame("balanced", [
        ("A02", 4000, 20000, "narrow_1", "8000", 1),
        ("A02", 20000, 36000, "narrow_2", "8000", 1),
    ])

    observed = unify_boundaries(raw, balanced, tie_tolerance_bp=4000)

    row = observed.iloc[0]
    assert row["component_state"] == "split_within_sample_fine_balanced"
    assert row["signal_support"] == "both"
    assert row["raw_boundary_ids"] == "wide"
    assert row["balanced_boundary_ids"] == "narrow_1,narrow_2"
    assert len(observed) == 1
```

- [ ] **Step 2: Run it and confirm it fails**

Run: `python3 -m pytest test/is_boundary/test_unification.py::test_wide_boundary_containing_two_narrow_ones_is_a_split_not_a_tie -q`
Expected: FAIL with `AssertionError` on `component_state == 'complex_within_sample'`.

- [ ] **Step 3: Implement the containment detector**

In `scripts/is_boundary_unification.py`, add above `unify_boundaries`:

```python
def _contained_children(parent, candidates):
    """Indices of candidates whose interval lies inside ``parent``.

    Uses half-open semantics: a child must satisfy
    ``parent.start <= child.start`` and ``child.end <= parent.end``.
    """
    return [
        index for index, item in enumerate(candidates)
        if int(parent.start) <= int(item.start) and int(item.end) <= int(parent.end)
    ]


def _containment_split(left, right, min_children=2):
    """Describe a coarse/fine split, or return ``None``.

    Returns ``(fine_signal, parent_position)`` where ``parent_position`` is one
    of ``"raw"``/``"balanced"`` and ``fine_signal`` the opposite side.
    """
    for signal, parent_frame, child_frame in (
        ("balanced", left, right), ("raw", right, left),
    ):
        for position in range(len(parent_frame)):
            children = _contained_children(parent_frame.iloc[position], child_frame)
            if len(children) >= min_children:
                return signal, position
    return None
```

- [ ] **Step 4: Replace the tie branch in `unify_boundaries`**

Replace the `if _has_tie(...)` block with:

```python
            split = _containment_split(component_left, component_right)
            if split is not None:
                fine_signal, _ = split
                output.append(_row(
                    component_left.itertuples(index=False),
                    component_right.itertuples(index=False),
                    f"split_within_sample_fine_{fine_signal}",
                ))
                continue
            was_tied = _has_tie(component_left, component_right, tie_tolerance_bp)
            pairs = balanced_priority_pairs(component_left, component_right)
```

The existing `was_tied` / `tie_resolution` block below stays as-is, so a genuine remaining tie still gets `balanced_priority_tied_component`.

- [ ] **Step 5: Run the test**

Run: `python3 -m pytest test/is_boundary/test_unification.py -q`
Expected: the new test passes. `test_adjacency_does_not_unify_and_many_to_many_is_complex` now sees `split_within_sample_fine_*` instead of `complex_within_sample`; update its assertion from

```python
    assert "complex_within_sample" not in set(observed["component_state"])
```

to

```python
    assert not set(observed["component_state"]).intersection({"complex_within_sample"})
```

and change the `tied = observed[...]` lookup to match `split_within_sample_fine_*` if that fixture triggers containment. Re-run until green.

- [ ] **Step 6: Commit**

```bash
git add scripts/is_boundary_unification.py test/is_boundary/test_unification.py
git commit -m "fix: treat a coarse boundary containing two fine ones as a split, not a tie"
```

---

### Task 2: Profile window and long-table rows

**Files:**
- Create: `scripts/is_boundary_profile.py`
- Create: `test/is_boundary/test_profile.py`
- Modify: `scripts/is_boundary_io.py` (add `PROFILE_COLUMNS`, `PROFILE_SUMMARY_COLUMNS`)

**Why:** Classification needs a per-sample, coordinate-bearing view of what boundaries exist near a locus. The long table is what makes an SV/CNV join possible later.

- [ ] **Step 1: Define the column schemas**

In `scripts/is_boundary_io.py`, after `ANNOTATION_COLUMNS`:

```python
PROFILE_COLUMNS = (
    "comparison", "region_id", "discovery_mode", "chrom", "locus_summit",
    "window_start", "window_end", "sample", "boundary_id", "start", "end",
    "center_bp", "center_offset_bp", "n_supporting_resolutions",
    "signal_support", "support_resolutions",
)
PROFILE_SUMMARY_COLUMNS = (
    "comparison", "region_id", "discovery_mode", "chrom", "locus_summit",
    "window_start", "window_end", "n_boundaries_a", "n_boundaries_b",
    "nearest_offset_bp", "boundaries_per_100kb_a", "boundaries_per_100kb_b",
    "profile_position_state", "profile_reason",
)
```

- [ ] **Step 2: Write the failing profile-row test**

Create `test/is_boundary/test_profile.py`:

```python
import pandas as pd

from is_boundary_profile import (
    locus_window, profile_rows, PROFILE_COLUMNS,
)


def _boundary(start, end, boundary_id, signal_support="both",
              state="simple", resolutions="4000", nsup=1, sample="S"):
    return {
        "chrom": "A02", "start": start, "end": end, "boundary_id": boundary_id,
        "center_bp": (start + end) // 2, "component_state": state,
        "signal_support": signal_support, "support_resolutions": resolutions,
        "n_supporting_resolutions": nsup, "n_source_boundaries": 1, "sample": sample,
    }


def test_locus_window_centres_on_the_anchor():
    locus = pd.Series({"is_change_type": "transition", "summit": 100000,
                       "event_start": 96000, "event_end": 104000})
    assert locus_window(locus, window_bp=100000, base_bin_bp=4000) == (52000, 152000)


def test_profile_rows_carry_coordinates_and_skip_unevaluable():
    boundaries = pd.DataFrame([
        _boundary(90000, 98000, "b1"),
        _boundary(110000, 118000, "b2"),
        _boundary(60000, 68000, "b3", state="unevaluable_coverage_raw_nonfinite_base_bins"),
    ])
    locus = pd.Series({"region_id": 7, "discovery_mode": "refined",
                       "is_change_type": "transition", "chrom": "A02",
                       "summit": 100000, "event_start": 96000, "event_end": 104000})

    rows = profile_rows(locus, boundaries, sample="H1398", comparison="A_vs_B",
                        window_bp=100000, base_bin_bp=4000)

    assert list(rows.columns) == list(PROFILE_COLUMNS)
    assert rows["boundary_id"].tolist() == ["b1", "b2"]
    assert rows["center_offset_bp"].tolist() == [-6000, 14000]
    assert rows["sample"].unique().tolist() == ["H1398"]
```

- [ ] **Step 3: Run it and confirm it fails**

Run: `python3 -m pytest test/is_boundary/test_profile.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'is_boundary_profile'`.

- [ ] **Step 4: Implement the module skeleton**

Create `scripts/is_boundary_profile.py`:

```python
"""Per-locus boundary profiles.

A boundary is an interval whose width is set by the resolution that called it,
not a point. Comparing boundary centres therefore mis-reports a coarse/fine
split as a positional tie. This module sidesteps matching entirely: for each
IS-change locus it counts the boundaries that fall in a fixed window, per
sample, and keeps every interval so downstream joins stay possible.
"""

import pandas as pd

from is_boundary_io import PROFILE_COLUMNS

UNEQUIVOCAL_EXCLUDE_PREFIX = "unevaluable_coverage_"


def locus_anchor(locus, base_bin_bp):
    """Anchor interval of a locus, matching the annotation stage."""
    if locus["is_change_type"] == "broad_shift":
        return int(locus["event_start"]), int(locus["event_end"])
    return int(locus["summit"]), int(locus["summit"]) + base_bin_bp


def locus_window(locus, window_bp, base_bin_bp):
    """Half-open window centred on the locus anchor's midpoint."""
    start, end = locus_anchor(locus, base_bin_bp)
    centre = (start + end) // 2
    half = window_bp // 2
    return max(0, centre - half), centre + half


def _overlaps(start, end, other_start, other_end):
    return max(start, other_start) < min(end, other_end)


def profile_rows(locus, boundaries, sample, comparison, window_bp, base_bin_bp):
    window_start, window_end = locus_window(locus, window_bp, base_bin_bp)
    peak = int(locus["summit"])
    frame = boundaries[boundaries["chrom"].eq(locus["chrom"])]
    rows = []
    for item in frame.itertuples(index=False):
        state = str(getattr(item, "component_state", ""))
        if state.startswith(UNEQUIVOCAL_EXCLUDE_PREFIX):
            continue
        start, end = int(item.start), int(item.end)
        if not _overlaps(window_start, window_end, start, end):
            continue
        rows.append({
            "comparison": comparison,
            "region_id": int(locus["region_id"]),
            "discovery_mode": str(locus["discovery_mode"]),
            "chrom": locus["chrom"],
            "locus_summit": peak,
            "window_start": window_start,
            "window_end": window_end,
            "sample": sample,
            "boundary_id": item.boundary_id,
            "start": start,
            "end": end,
            "center_bp": int(item.center_bp),
            "center_offset_bp": int(item.center_bp) - peak,
            "n_supporting_resolutions": int(item.n_supporting_resolutions),
            "signal_support": item.signal_support,
            "support_resolutions": item.support_resolutions,
        })
    frame = pd.DataFrame(rows, columns=PROFILE_COLUMNS)
    if frame.empty:
        return frame
    return frame.sort_values(
        ["chrom", "window_start", "center_bp", "boundary_id"], kind="stable"
    ).reset_index(drop=True)
```

- [ ] **Step 5: Run the tests**

Run: `python3 -m pytest test/is_boundary/test_profile.py -q`
Expected: 2 passed.

- [ ] **Step 6: Commit**

```bash
git add scripts/is_boundary_profile.py scripts/is_boundary_io.py test/is_boundary/test_profile.py
git commit -m "feat: add per-locus boundary profile rows with coordinates"
```

---

### Task 3: Classify a locus from its two boundary counts

**Files:**
- Modify: `scripts/is_boundary_profile.py`
- Test: `test/is_boundary/test_profile.py`

**Why:** This is the classification that replaces the matching verdict. It never needs to decide which A boundary corresponds to which B boundary, so it cannot produce a tie.

- [ ] **Step 1: Write the failing classification test**

Append to `test/is_boundary/test_profile.py`:

```python
from is_boundary_profile import classify_profile, summarize_profile


def test_one_boundary_becomes_two_is_split_in_b():
    state, reason = classify_profile(
        n_a=1, n_b=2, nearest_offset_bp=15000,
        shift_tolerance_bp=20000, dense_threshold=20.0, window_bp=100000,
    )
    assert state == "profile_split_in_B"
    assert reason == "boundary_count_1_to_2"


def test_two_become_one_is_merge_in_b():
    state, _ = classify_profile(
        n_a=2, n_b=1, nearest_offset_bp=15000,
        shift_tolerance_bp=20000, dense_threshold=20.0, window_bp=100000,
    )
    assert state == "profile_merge_in_B"


def test_offset_within_tolerance_is_conserved():
    state, _ = classify_profile(
        n_a=1, n_b=1, nearest_offset_bp=12000,
        shift_tolerance_bp=20000, dense_threshold=20.0, window_bp=100000,
    )
    assert state == "profile_conserved"


def test_offset_beyond_tolerance_is_shifted():
    state, _ = classify_profile(
        n_a=1, n_b=1, nearest_offset_bp=34000,
        shift_tolerance_bp=20000, dense_threshold=20.0, window_bp=100000,
    )
    assert state == "profile_shifted"


def test_empty_window_is_no_boundary_not_complex():
    state, _ = classify_profile(
        n_a=0, n_b=0, nearest_offset_bp=float("nan"),
        shift_tolerance_bp=20000, dense_threshold=20.0, window_bp=100000,
    )
    assert state == "profile_no_boundary"


def test_dense_region_overrides_count_arithmetic():
    state, reason = classify_profile(
        n_a=41, n_b=3, nearest_offset_bp=2000,
        shift_tolerance_bp=20000, dense_threshold=20.0, window_bp=100000,
    )
    assert state == "profile_dense_region"
    assert reason.startswith("boundaries_per_100kb=")
```

- [ ] **Step 2: Run it and confirm it fails**

Run: `python3 -m pytest test/is_boundary/test_profile.py -q -k classify or dense or offset or empty`
Expected: FAIL with `ImportError: cannot import name 'classify_profile'`.

- [ ] **Step 3: Implement classification and summary**

Append to `scripts/is_boundary_profile.py`:

```python
DENSE_PRECEDENCE_STATES = ("profile_dense_region",)


def _per_100kb(count, window_bp):
    return count * 100000.0 / window_bp if window_bp else 0.0


def classify_profile(n_a, n_b, nearest_offset_bp, shift_tolerance_bp,
                     dense_threshold, window_bp):
    """Classify a locus from boundary counts alone.

    Density is checked before count arithmetic: a window packed with
    boundaries cannot support a statement about one of them.
    """
    density_a = _per_100kb(n_a, window_bp)
    density_b = _per_100kb(n_b, window_bp)
    if max(density_a, density_b) >= dense_threshold:
        return (
            "profile_dense_region",
            f"boundaries_per_100kb={max(density_a, density_b):.1f}",
        )
    if n_a == 0 and n_b == 0:
        return "profile_no_boundary", "empty_window_both_samples"
    if n_a == 0:
        return "profile_gained_in_B", "boundaries_0_to_{}".format(n_b)
    if n_b == 0:
        return "profile_lost_in_B", "boundaries_{}_to_0".format(n_a)
    if n_a == 1 and n_b == 1:
        if nearest_offset_bp <= shift_tolerance_bp:
            return "profile_conserved", "single_boundary_within_tolerance"
        return "profile_shifted", "single_boundary_beyond_tolerance"
    if n_a < n_b:
        return "profile_split_in_B", f"boundary_count_{n_a}_to_{n_b}"
    if n_a > n_b:
        return "profile_merge_in_B", f"boundary_count_{n_a}_to_{n_b}"
    return "profile_multi_boundary", f"boundary_count_{n_a}_to_{n_b}"


def summarize_profile(rows, window_bp, shift_tolerance_bp, dense_threshold):
    """Collapse long-table rows into one classification per locus."""
    from is_boundary_io import PROFILE_SUMMARY_COLUMNS

    keys = ["comparison", "region_id", "discovery_mode", "chrom", "locus_summit",
            "window_start", "window_end"]
    summary = []
    if rows.empty:
        return pd.DataFrame(columns=PROFILE_SUMMARY_COLUMNS)
    for key, group in rows.groupby(keys, sort=False, observed=True):
        counts = group["sample"].value_counts()
        n_a = int(counts.get("H1398", 0))
        n_b = int(counts.get("L134", 0))
        offsets = group["center_offset_bp"].abs()
        nearest = float(offsets.min()) if len(offsets) else float("nan")
        state, reason = classify_profile(
            n_a, n_b, nearest, shift_tolerance_bp, dense_threshold, window_bp
        )
        summary.append({
            **dict(zip(keys, key)),
            "n_boundaries_a": n_a,
            "n_boundaries_b": n_b,
            "nearest_offset_bp": nearest,
            "boundaries_per_100kb_a": _per_100kb(n_a, window_bp),
            "boundaries_per_100kb_b": _per_100kb(n_b, window_bp),
            "profile_position_state": state,
            "profile_reason": reason,
        })
    return pd.DataFrame(summary, columns=PROFILE_SUMMARY_COLUMNS)
```

Note: the `H1398`/`L134` literals are wrong for any other comparison. Replace them with `sample_a`/`sample_b` arguments in Step 4.

- [ ] **Step 4: Make the sample names parameters**

Change `summarize_profile(rows, window_bp, shift_tolerance_bp, dense_threshold)` to

```python
def summarize_profile(rows, window_bp, shift_tolerance_bp, dense_threshold,
                      sample_a, sample_b):
```

and replace the two literals with `counts.get(sample_a, 0)` / `counts.get(sample_b, 0)`. Add a test:

```python
def test_summary_uses_the_supplied_sample_names():
    rows = pd.DataFrame([
        {"comparison": "C", "region_id": 1, "discovery_mode": "refined",
         "chrom": "A02", "locus_summit": 100000, "window_start": 50000,
         "window_end": 150000, "sample": "S1", "center_offset_bp": 0},
        {"comparison": "C", "region_id": 1, "discovery_mode": "refined",
         "chrom": "A02", "locus_summit": 100000, "window_start": 50000,
         "window_end": 150000, "sample": "S2", "center_offset_bp": 1000},
        {"comparison": "C", "region_id": 1, "discovery_mode": "refined",
         "chrom": "A02", "locus_summit": 100000, "window_start": 50000,
         "window_end": 150000, "sample": "S2", "center_offset_bp": 9000},
    ])
    summary = summarize_profile(rows, 100000, 20000, 20.0, "S1", "S2")
    row = summary.iloc[0]
    assert (row["n_boundaries_a"], row["n_boundaries_b"]) == (1, 2)
    assert row["profile_position_state"] == "profile_split_in_B"
```

- [ ] **Step 5: Run the tests**

Run: `python3 -m pytest test/is_boundary/test_profile.py -q`
Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add scripts/is_boundary_profile.py test/is_boundary/test_profile.py
git commit -m "feat: classify loci from per-sample boundary counts with density precedence"
```

---

### Task 4: Wire the profile into the annotation CLI

**Files:**
- Modify: `scripts/is_boundary_config.py:27-40` (`ISBoundarySettings`, parser)
- Modify: `config/config.yaml` (`is_boundary_annotation.profile`)
- Modify: `scripts/is_boundary_annotation.py` (append profile columns)
- Modify: `scripts/annotate_is_change_boundaries.py:61-73` (`_settings`), `build_parser`, `run_from_args`
- Modify: `rules/is_boundary_annotation.smk` (annotate rule inputs/outputs/shell)

**Why:** A module nobody calls changes nothing. This task makes the profile reach the locus tables.

- [ ] **Step 1: Add configuration**

In `config/config.yaml` under `is_boundary_annotation:`:

```yaml
  profile:
    window_bp: 100000
    shift_tolerance_bp: 20000
    dense_boundaries_per_100kb: 20.0
```

In `scripts/is_boundary_config.py`, add:

```python
@dataclass(frozen=True)
class ProfileSettings:
    window_bp: int
    shift_tolerance_bp: int
    dense_boundaries_per_100kb: float
```

Add `profile: ProfileSettings` as the last field of `ISBoundarySettings`, and parse it in `parse_is_boundary_settings` with an explicit `_required` lookup per field. Add a config test asserting that a missing `profile` block raises `IS_BOUNDARY_CONFIG_INVALID`.

- [ ] **Step 2: Append profile columns to annotated loci**

In `scripts/is_boundary_annotation.py`, add to `ANNOTATION_COLUMNS`-adjacent constants a new tuple:

```python
PROFILE_ANNOTATION_COLUMNS = (
    "profile_position_state", "profile_reason", "profile_n_boundaries_a",
    "profile_n_boundaries_b", "profile_nearest_offset_bp",
    "profile_boundaries_per_100kb_a", "profile_boundaries_per_100kb_b",
)
```

Then add a function:

```python
def attach_profile(annotated, summary):
    """Join the per-locus profile classification onto annotated loci.

    The profile never overwrites ``position_state``; it is published beside it
    until the two are shown to agree where they should.
    """
    keys = ["region_id", "discovery_mode"]
    lookup = summary.set_index(keys)
    frame = annotated.copy()
    for column in PROFILE_ANNOTATION_COLUMNS:
        source = column.removeprefix("profile_")
        frame[column] = frame.set_index(keys).index.map(lookup[source])
    return frame.reset_index(drop=True)
```

Add a test that a locus with no profile summary row gets `NaN` in all seven columns and keeps its original `position_state`.

- [ ] **Step 3: Wire the CLI**

In `scripts/annotate_is_change_boundaries.py`:
- Add `--boundaries-a`, `--boundaries-b` required arguments.
- Add outputs `--out-profile` and `--out-profile-summary`.
- After `annotate_loci`, build the long table:

```python
        profile_frames = []
        for sample, path in ((sample_a, args.boundaries_a), (sample_b, args.boundaries_b)):
            boundaries = read_unified_boundaries(path)
            for frame in (primary, preliminary):
                for locus in frame.itertuples(index=False):
                    profile_frames.append(profile_rows(
                        locus._asdict(), boundaries, sample=sample,
                        comparison=args.comparison,
                        window_bp=boundary_settings.profile.window_bp,
                        base_bin_bp=boundary_settings.base_bin_bp,
                    ))
        profile = pd.concat(profile_frames, ignore_index=True)
        profile_summary = summarize_profile(
            profile, boundary_settings.profile.window_bp,
            boundary_settings.profile.shift_tolerance_bp,
            boundary_settings.profile.dense_boundaries_per_100kb,
            sample_a, sample_b,
        )
        annotated_primary = attach_profile(annotated_primary, profile_summary)
        annotated_preliminary = attach_profile(annotated_preliminary, profile_summary)
```

- [ ] **Step 4: Update the Snakemake rule**

In `rules/is_boundary_annotation.smk`, add to the annotate rule:

```yaml
        boundaries_a=lambda w: unified_boundary_path(get_pair_samples(w.comparison)[0]),
        boundaries_b=lambda w: unified_boundary_path(get_pair_samples(w.comparison)[1]),
```

Add to `output`:

```yaml
        profile=OUTDIR + "/tad_boundary_events/{comparison}.boundary_profile.tsv",
        profile_summary=OUTDIR + "/tad_boundary_events/{comparison}.boundary_profile_summary.tsv",
```

and the matching `--boundaries-a`, `--boundaries-b`, `--out-profile`, `--out-profile-summary` flags in the shell block.

- [ ] **Step 5: Run the boundary test suite**

Run: `python3 -m pytest test/is_boundary/ -q`
Expected: all pass except the two known pre-existing/flaky failures (`test_plot_rule_reverts_to_unannotated_primary_when_disabled` is pre-existing; the matching CLI test passes in isolation).

- [ ] **Step 6: Commit**

```bash
git add config/config.yaml scripts/is_boundary_config.py scripts/is_boundary_annotation.py \
        scripts/annotate_is_change_boundaries.py rules/is_boundary_annotation.smk \
        test/is_boundary/
git commit -m "feat: publish per-locus boundary profiles beside the event-based state"
```

---

### Task 5: Validate against the current states, then rerun

**Files:**
- Read: `results/tad_boundary_events/H1398_vs_L134.boundary_profile_summary.tsv` (produced by Task 4)
- Read: `results/multiscale/r4000/H1398_vs_L134/is_change/*_boundary_annotated.tsv`
- Create: `results/audits/<timestamp>/profile_vs_event_agreement.tsv`

**Why:** Do not promote the profile to `position_state` until it agrees with the event pipeline on the loci where both must agree. Where they disagree, the disagreement is evidence about which one is wrong.

- [ ] **Step 1: Submit the rerun**

Write `scripts/rerun_boundary_profile.dsub` following `scripts/rerun_balanced_priority_annotation.dsub`, with `--allowed-rules unify_tad_boundaries call_tad_boundary_events annotate_is_change_boundaries` and the annotated primary as the target. Preserve the previous outputs in an audit directory before running.

- [ ] **Step 2: Measure agreement**

Compare `profile_position_state` against `position_state` on the 1,572 primary loci:

| existing `position_state` | expected profile state | tolerance |
|---|---|---|
| `conserved` | `profile_conserved` | any agreement rate below 85% means the window or tolerance is wrong |
| `shifted` | `profile_shifted` | same |
| `gained_in_B` | `profile_gained_in_B` | same |
| `lost_in_B` | `profile_lost_in_B` | same |
| `no_boundary_association` | `profile_no_boundary` | same |
| `ambiguous` | any | masks dominate; report only |
| `complex` | any | this is the population being replaced |

- [ ] **Step 3: Report, do not promote yet**

Write the agreement table and hand it to the user before changing `position_state`. Promotion is a separate decision that depends on the agreement numbers.

- [ ] **Step 4: Commit the validation script**

```bash
git add results/audits/*/profile_vs_event_agreement.tsv scripts/rerun_boundary_profile.dsub
git commit -m "chore: add boundary profile rerun and agreement audit"
```

---

## Calibration notes

These values are starting points and must be checked against the Task 5 output:

- `window_bp: 100000` — at the genome-wide rate of roughly 1.8 unified boundaries per 100 kb, a 100 kb window holds about two boundaries by chance. Wider windows make every locus look dense; narrower windows make the 20 kb shift tolerance consume most of the window.
- `shift_tolerance_bp: 20000` — matches the measured raw-vs-balanced centre-offset fallback and is below the 22–24 kb median rival separation measured in tied components, so a genuine split stays a split.
- `dense_boundaries_per_100kb: 20.0` — roughly 11× the genome-wide rate. Ties observed so far reach 37 boundaries per component; this threshold should separate them from ordinary loci. Verify the resulting density distribution before fixing the value.
