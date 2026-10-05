# IS-change TAD Boundary Annotation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a provenance-aware boundary annotation stage that classifies each primary or preliminary IS-change locus as conserved, gained, lost, shifted, complex, ambiguous, or not boundary-associated, with a separate multi-scale IS boundary-strength state.

**Architecture:** Add three independently retryable workflow stages: unify raw/balanced multi-resolution boundaries once per sample, construct an ordered genome-wide boundary-event map once per comparison, and annotate primary/preliminary IS loci from that event map and the existing balanced IS tracks. Keep matching, strength, association, and optional permutation logic in small pure-function modules; CLIs own schema validation and atomic publication, while Snakemake owns dependencies and output paths.

**Tech Stack:** Python 3.10+, pandas, NumPy, SciPy `signal.find_peaks` and `signal.peak_prominences`, PyYAML, pytest, Snakemake 8+, existing `env/SNAP-HiC-Diff.yaml`.

## Global Constraints

- The approved specification is `docs/superpowers/specs/2026-09-03-is-change-tad-boundary-annotation-design.md`.
- H1398 is sample A and L134 is sample B for `H1398_vs_L134`; all directional labels are A-to-B.
- Both samples use the same reference genome and the same 4,000 bp IS grid.
- Analyze boundary positions only; do not infer complete-TAD split, fusion, merge, or nested topology.
- Preserve the existing raw and balanced 3DNetMod outputs. The unified boundary table is a downstream annotation view and must not replace either source family.
- Retain unmatched raw and balanced boundaries with `raw_only` or `balanced_only` provenance.
- Main summaries use `primary` loci. `single_band_preliminary` loci receive identical annotations in a separate output.
- Use `J = max(4,000 bp, Q95 of eligible technical offsets)` with at least 100 pooled offsets; otherwise use the documented 20,000 bp fallback.
- Use 100,000 bp as the primary maximum shift and publish 60,000/80,000/100,000 bp sensitivity results.
- Do not use the old-to-new shifted event span as an IS-overlap target; only endpoint support intervals count.
- Use local prominence of `-IS`, not a single raw IS value, as the boundary-strength statistic.
- Preserve the existing IS scale bands exactly: small = 32/60/92 kb, medium =
  120/152/180 kb, and large = 212/240/272/300 kb.
- Keep `position_state`, `strength_state`, `is_consistency_state`, and `evidence_tier` as separate fields.
- Emit no locus-level inferential P values. The optional 10,000-permutation analysis is set-level and disabled by default.
- All outputs are deterministic, sorted, schema-checked, and atomically replaced.
- Missing coverage produces `ambiguous` or `unevaluable`, never gain or loss.
- Use stable leading error codes exactly:
  `IS_BOUNDARY_CONFIG_INVALID`, `IS_BOUNDARY_INPUT_SCHEMA_INVALID`,
  `IS_BOUNDARY_REFERENCE_MISMATCH`, `IS_BOUNDARY_JITTER_CALIBRATION_FAILED`,
  `IS_BOUNDARY_MATCH_FAILED`, `IS_BOUNDARY_STRENGTH_FAILED`, and
  `IS_BOUNDARY_ANNOTATION_FAILED`; scientific insufficiency is written to QC
  and does not abort unaffected chromosomes.
- No new Conda dependencies are required; SciPy, pandas, NumPy, PyYAML, and pytest already exist in `env/SNAP-HiC-Diff.yaml`.

---

## File Structure

### Create

- `scripts/is_boundary_errors.py` — stable error-code exception type.
- `scripts/is_boundary_config.py` — immutable configuration objects and cross-checks against existing IS/TAD settings.
- `scripts/is_boundary_io.py` — BED/TSV parsing, output schemas, deterministic sorting, and atomic table publication.
- `scripts/is_boundary_unification.py` — raw/balanced overlap assignment, complex-component detection, technical-offset extraction, and source/leave-one-resolution-out views.
- `scripts/unify_tad_boundaries.py` — per-sample unification CLI.
- `scripts/is_boundary_matching.py` — jitter estimation, ordered dynamic programming, uniqueness diagnostics, event construction, and view-level stability comparison.
- `scripts/call_tad_boundary_events.py` — per-comparison boundary-event CLI.
- `scripts/is_boundary_strength.py` — local valley prominence, robust standardization, band aggregation, and direction support.
- `scripts/is_boundary_annotation.py` — IS-anchor association, evaluability rules, state inheritance, evidence tier, and summaries.
- `scripts/annotate_is_change_boundaries.py` — primary/preliminary annotation CLI.
- `scripts/is_boundary_permutation.py` — constrained candidate-region randomization and empirical set-level statistics.
- `scripts/permute_is_boundary_association.py` — optional permutation CLI.
- `rules/is_boundary_annotation.smk` — four Snakemake rules and their stable logs/benchmarks.
- `test/is_boundary/test_config.py` — configuration contract tests.
- `test/is_boundary/test_unification.py` — boundary parsing, union, provenance, and complex-component tests.
- `test/is_boundary/test_matching.py` — jitter, DP, shift, gain/loss, ambiguity, and direction-swap tests.
- `test/is_boundary/test_strength.py` — prominence and multi-band strength tests.
- `test/is_boundary/test_annotation.py` — endpoint association, missingness, and evidence-tier tests.
- `test/is_boundary/test_permutation.py` — constrained-null and empirical-P tests.
- `test/is_boundary/test_workflow.py` — Snakemake rule and dry-run contracts.
- `test/is_boundary/test_synthetic_benchmark.py` — deterministic and empirical-jitter truth benchmark.
- `docs/adr/0002-downstream-unified-tad-boundary-view.md` — compatibility decision relative to ADR-0001.

### Modify

- `config/config.yaml` — add the approved `is_boundary_annotation` section.
- `Snakefile` — parse settings, add path helpers, include rules, and add stable outputs to `rule all`.
- `scripts/plot_is_change_loci.py` — consume annotated primary loci and overlay A/B boundary endpoint intervals.
- `rules/is_change.smk` — make plotting consume the annotated primary table and boundary-event map.
- `test/is_change/test_plot_is_change_loci.py` — test boundary overlays and annotated schema.
- `test/is_change/test_is_change_workflow.py` — update the plotting dependency contract.
- `CONTEXT.md` — define unified boundary, boundary position state, shifted boundary, and IS boundary strength.
- `README.md` — document outputs, interpretation limits, targeted execution, and QC files.

## Cross-task Public Interfaces

These names and signatures are fixed for the rest of the plan:

| Module | Public interface |
|---|---|
| `is_boundary_io.py` | `read_multires_boundaries(path: str or Path, signal: str) -> pd.DataFrame` |
| `is_boundary_io.py` | `read_resolution_boundaries(paths: Sequence of resolution/path pairs, signal: str) -> pd.DataFrame` |
| `is_boundary_io.py` | `atomic_write_tsv(frame: pd.DataFrame, path: str or Path, columns: Sequence of str) -> None` |
| `is_boundary_unification.py` | `unify_boundaries(raw: pd.DataFrame, balanced: pd.DataFrame, tie_tolerance_bp: int, chrom_order: Sequence of str or None = None) -> pd.DataFrame` |
| `is_boundary_unification.py` | `build_jitter_offsets(raw: pd.DataFrame, balanced: pd.DataFrame, per_resolution: pd.DataFrame, base_bin_bp: int) -> pd.DataFrame` |
| `is_boundary_unification.py` | `build_sensitivity_views(sample: str, raw: pd.DataFrame, balanced: pd.DataFrame, per_resolution: pd.DataFrame, settings: ISBoundarySettings) -> pd.DataFrame` |
| `is_boundary_matching.py` | `estimate_jitter(offsets_a: pd.DataFrame, offsets_b: pd.DataFrame, settings: ISBoundarySettings) -> JitterEstimate` |
| `is_boundary_matching.py` | `align_ordered(boundaries_a: pd.DataFrame, boundaries_b: pd.DataFrame, max_distance_bp: int) -> list of integer index pairs` |
| `is_boundary_matching.py` | `call_boundary_events(boundaries_a: pd.DataFrame, boundaries_b: pd.DataFrame, jitter: JitterEstimate, settings: ISBoundarySettings, comparison: str, sample_a: str, sample_b: str) -> pd.DataFrame` |
| `is_boundary_matching.py` | `compare_event_views(reference_events: pd.DataFrame, views_a: pd.DataFrame, views_b: pd.DataFrame, jitter: JitterEstimate, settings: ISBoundarySettings, comparison: str, sample_a: str, sample_b: str) -> pd.DataFrame` |
| `is_boundary_strength.py` | `compute_strength_long(events: pd.DataFrame, track_a: pd.DataFrame, track_b: pd.DataFrame, boundary_settings: ISBoundarySettings, is_settings: ISChangeSettings) -> pd.DataFrame` |
| `is_boundary_strength.py` | `summarize_strength(strength_long: pd.DataFrame, events: pd.DataFrame, boundary_settings: ISBoundarySettings, is_settings: ISChangeSettings) -> pd.DataFrame` |
| `is_boundary_annotation.py` | `annotate_loci(loci: pd.DataFrame, events: pd.DataFrame, strength_summary: pd.DataFrame, masked: pd.DataFrame, settings: ISBoundarySettings) -> pd.DataFrame` |
| `is_boundary_annotation.py` | `summarize_annotations(primary: pd.DataFrame, preliminary: pd.DataFrame) -> pd.DataFrame` |
| `is_boundary_permutation.py` | `randomize_loci_once(loci: pd.DataFrame, candidate_regions: pd.DataFrame, masked: pd.DataFrame, settings: ISBoundarySettings, rng: np.random.Generator) -> pd.DataFrame` |
| `is_boundary_permutation.py` | `run_permutations(loci: pd.DataFrame, candidate_regions: pd.DataFrame, events: pd.DataFrame, masked: pd.DataFrame, settings: ISBoundarySettings) -> pd.DataFrame` |

The schema constants in `is_boundary_io.py` must contain exactly these fields:

```text
MULTIRES_INPUT_COLUMNS =
chrom,start,end,boundary_id,support_resolutions,
n_supporting_resolutions,n_source_boundaries

UNIFIED_BOUNDARY_COLUMNS =
chrom,start,end,unified_boundary_id,center_bp,component_state,
signal_support,raw_boundary_ids,balanced_boundary_ids,support_resolutions,
n_supporting_resolutions,n_source_boundaries

BOUNDARY_VIEW_COLUMNS =
sample,view,chrom,start,end,unified_boundary_id,center_bp,component_state,
signal_support,raw_boundary_ids,balanced_boundary_ids,support_resolutions,
n_supporting_resolutions,n_source_boundaries

BOUNDARY_EVENT_COLUMNS =
comparison,sample_a,sample_b,chrom,boundary_event_id,position_state,
position_substate,boundary_a_id,boundary_b_id,a_start,a_end,b_start,b_end,
event_span_start,event_span_end,shift_bp,match_stage,match_is_reciprocal,
second_best_margin_bp,n_local_candidates,jitter_threshold_bp,jitter_source,
max_shift_bp,ambiguity_reason,state_at_60000,state_at_80000,state_at_100000,
position_state_stable,signal_support_a,signal_support_b,resolution_support_a,
resolution_support_b

ANNOTATION_COLUMNS =
boundary_event_id,position_state,position_substate,strength_state,
strength_support_bands,is_consistency_state,evidence_tier,boundary_a_id,
boundary_b_id,shift_bp,signal_support_a,signal_support_b,
resolution_support_a,resolution_support_b,n_boundary_candidates,
jitter_threshold_bp,max_shift_bp,annotation_reason

STRENGTH_LONG_COLUMNS =
comparison,boundary_event_id,chrom,position_state,coordinate_role,IS_distance,
scale_band,a_start,a_end,b_start,b_end,strength_a,strength_b,delta_strength,
robust_z,valley_a_start,valley_b_start,detection_a,detection_b,
valid_fraction_a,valid_fraction_b,reference_n,reference_median,reference_mad

STRENGTH_SUMMARY_COLUMNS =
comparison,boundary_event_id,position_state,strength_state,
strength_support_bands,is_consistency_state,small_band_score,
medium_band_score,large_band_score,n_evaluable_bands,state_at_z1_5,
state_at_z2_0,state_at_z2_5,strength_state_stable

BOUNDARY_STABILITY_COLUMNS =
comparison,boundary_event_id,view,reference_position_state,
sensitivity_position_state,concordant,reason
```

`build_jitter_offsets()` returns exactly:

```text
sample,chrom,source,comparison,resolution_a,resolution_b,boundary_a_id,
boundary_b_id,center_a,center_b,offset_bp,eligible
```

No task may rename a function or schema field without updating every later
consumer and test in the same commit.

---

### Task 1: Configuration, Errors, and Domain Decision

**Files:**

- Create: `scripts/is_boundary_errors.py`
- Create: `scripts/is_boundary_config.py`
- Create: `test/is_boundary/test_config.py`
- Create: `docs/adr/0002-downstream-unified-tad-boundary-view.md`
- Modify: `config/config.yaml:147-220`
- Modify: `CONTEXT.md:1-24`

**Interfaces:**

- Consumes: `CworldISSettings` from `scripts/cworld_is_config.py` and `TadSettings` from `scripts/tad_config.py`.
- Produces: `ISBoundarySettings`, `StrengthSettings`, `PermutationSettings`, `parse_is_boundary_settings()`, and `ISBoundaryError`.

- [ ] **Step 1: Write the failing repository-contract test**

Create `test/is_boundary/test_config.py` with a repository-config test and
parameterized invalid cases:

```python
from copy import deepcopy
from pathlib import Path
import sys

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from cworld_is_config import parse_cworld_is_settings
from is_boundary_config import parse_is_boundary_settings
from tad_config import parse_tad_settings


def repository_config():
    return yaml.safe_load((REPO / "config" / "config.yaml").read_text())


def test_repository_boundary_annotation_contract_is_exact():
    config = repository_config()
    settings = parse_is_boundary_settings(
        config,
        parse_cworld_is_settings(config),
        parse_tad_settings(config),
    )
    assert settings.base_bin_bp == 4000
    assert settings.jitter_quantile == 0.95
    assert settings.jitter_min_pairs == 100
    assert settings.jitter_fallback_bp == 20000
    assert settings.max_shift_bp == 100000
    assert settings.shift_sensitivity_bp == (60000, 80000, 100000)
    assert settings.match_tie_tolerance_bp == 4000
    assert settings.strength.robust_z_threshold == 2.0
    assert settings.strength.stable_abs_z_threshold == 1.0
    assert settings.permutation.enabled is False
    assert settings.permutation.n_permutations == 10000


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("base_bin_bp", 8000, "base_bin_bp"),
        ("jitter_quantile", 1.0, "jitter_quantile"),
        ("jitter_min_pairs", 0, "jitter_min_pairs"),
        ("max_shift_bp", 90000, "shift_sensitivity_bp"),
        ("match_tie_tolerance_bp", 2000, "4 kb"),
    ],
)
def test_rejects_invalid_boundary_annotation_values(field, value, message):
    config = deepcopy(repository_config())
    config["is_boundary_annotation"][field] = value
    with pytest.raises(ValueError, match=message):
        parse_is_boundary_settings(
            config,
            parse_cworld_is_settings(config),
            parse_tad_settings(config),
        )


def test_permutation_cannot_be_enabled_when_annotation_is_disabled():
    config = deepcopy(repository_config())
    config["is_boundary_annotation"]["enabled"] = False
    config["is_boundary_annotation"]["permutation"]["enabled"] = True
    with pytest.raises(ValueError, match="permutation.enabled"):
        parse_is_boundary_settings(
            config,
            parse_cworld_is_settings(config),
            parse_tad_settings(config),
        )
```

- [ ] **Step 2: Run the config test and confirm the missing module failure**

Run:

```bash
python -m pytest -q test/is_boundary/test_config.py
```

Expected: collection fails with `ModuleNotFoundError: No module named 'is_boundary_config'`.

- [ ] **Step 3: Add the stable exception and immutable settings**

Create `scripts/is_boundary_errors.py`:

```python
class ISBoundaryError(RuntimeError):
    def __init__(self, code, **context):
        self.code = code
        self.context = context
        details = " ".join(
            f"{key}={value}" for key, value in context.items() if value is not None
        )
        super().__init__(f"ERROR[{code}] {details}".rstrip())
```

Create `scripts/is_boundary_config.py` with frozen dataclasses and these public
fields:

```python
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from cworld_is_config import CworldISSettings
from tad_config import TadSettings


@dataclass(frozen=True)
class StrengthSettings:
    local_radius_scale_factor: float
    min_reference_boundaries: int
    min_valid_scales_per_band: int
    robust_z_threshold: float
    stable_abs_z_threshold: float
    sensitivity_thresholds: Sequence[float]


@dataclass(frozen=True)
class PermutationSettings:
    enabled: bool
    n_permutations: int
    seed: int


@dataclass(frozen=True)
class ISBoundarySettings:
    enabled: bool
    association_padding_bp: int
    raw_balanced_overlap_only: bool
    base_bin_bp: int
    jitter_quantile: float
    jitter_min_pairs: int
    jitter_fallback_bp: int
    max_shift_bp: int
    shift_sensitivity_bp: Sequence[int]
    match_tie_tolerance_bp: int
    strength: StrengthSettings
    permutation: PermutationSettings


def _positive_int(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"ERROR[IS_BOUNDARY_CONFIG_INVALID] {name} must be positive")
    return value


def parse_is_boundary_settings(
    config: Mapping[str, Any],
    cworld: CworldISSettings,
    tad: TadSettings,
) -> ISBoundarySettings:
    try:
        raw = config["is_boundary_annotation"]
        strength_raw = raw["strength"]
        permutation_raw = raw["permutation"]
    except (KeyError, TypeError) as error:
        raise ValueError(
            f"ERROR[IS_BOUNDARY_CONFIG_INVALID] missing={error}"
        ) from error
    base = _positive_int(raw["base_bin_bp"], "base_bin_bp")
    if base != cworld.source_resolution_bp or base != min(tad.resolutions):
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] base_bin_bp must equal the 4 kb IS and finest TAD resolution"
        )
    if raw["association_padding_bp"] != 0:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] association_padding_bp must equal zero"
        )
    if raw["raw_balanced_overlap_only"] is not True:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] raw_balanced_overlap_only must be true"
        )
    quantile = raw["jitter_quantile"]
    if isinstance(quantile, bool) or not isinstance(quantile, (int, float)) or not 0 < quantile < 1:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] jitter_quantile must be between zero and one"
        )
    sensitivity = tuple(raw["shift_sensitivity_bp"])
    if sensitivity != tuple(sorted(set(sensitivity))):
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] shift_sensitivity_bp must be sorted and unique"
        )
    max_shift = _positive_int(raw["max_shift_bp"], "max_shift_bp")
    if max_shift not in sensitivity:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] shift_sensitivity_bp must contain max_shift_bp"
        )
    grid_values = sensitivity + (
        raw["jitter_fallback_bp"], raw["match_tie_tolerance_bp"], max_shift,
    )
    if any(value % base for value in grid_values):
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] distance values must align to the 4 kb grid"
        )
    strength = StrengthSettings(
        local_radius_scale_factor=float(strength_raw["local_radius_scale_factor"]),
        min_reference_boundaries=_positive_int(
            strength_raw["min_reference_boundaries"], "min_reference_boundaries"
        ),
        min_valid_scales_per_band=_positive_int(
            strength_raw["min_valid_scales_per_band"], "min_valid_scales_per_band"
        ),
        robust_z_threshold=float(strength_raw["robust_z_threshold"]),
        stable_abs_z_threshold=float(strength_raw["stable_abs_z_threshold"]),
        sensitivity_thresholds=tuple(strength_raw["sensitivity_thresholds"]),
    )
    if strength.local_radius_scale_factor <= 0:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] local_radius_scale_factor must be positive"
        )
    if strength.stable_abs_z_threshold >= strength.robust_z_threshold:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] stable_abs_z_threshold must be smaller than robust_z_threshold"
        )
    if (
        tuple(strength.sensitivity_thresholds)
        != tuple(sorted(set(strength.sensitivity_thresholds)))
        or strength.robust_z_threshold not in strength.sensitivity_thresholds
    ):
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] sensitivity_thresholds must be sorted, unique, and contain robust_z_threshold"
        )
    permutation = PermutationSettings(
        enabled=permutation_raw["enabled"],
        n_permutations=_positive_int(
            permutation_raw["n_permutations"], "n_permutations"
        ),
        seed=_positive_int(permutation_raw["seed"], "seed"),
    )
    if not isinstance(raw["enabled"], bool) or not isinstance(permutation.enabled, bool):
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] enabled fields must be boolean"
        )
    if permutation.enabled and not raw["enabled"]:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] permutation.enabled requires is_boundary_annotation.enabled"
        )
    jitter_fallback = _positive_int(raw["jitter_fallback_bp"], "jitter_fallback_bp")
    if jitter_fallback != max(tad.resolutions):
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] jitter_fallback_bp must equal the maximum TAD resolution"
        )
    return ISBoundarySettings(
        enabled=raw["enabled"],
        association_padding_bp=raw["association_padding_bp"],
        raw_balanced_overlap_only=True,
        base_bin_bp=base,
        jitter_quantile=float(quantile),
        jitter_min_pairs=_positive_int(raw["jitter_min_pairs"], "jitter_min_pairs"),
        jitter_fallback_bp=jitter_fallback,
        max_shift_bp=max_shift,
        shift_sensitivity_bp=sensitivity,
        match_tie_tolerance_bp=_positive_int(
            raw["match_tie_tolerance_bp"], "match_tie_tolerance_bp"
        ),
        strength=strength,
        permutation=permutation,
    )
```

- [ ] **Step 4: Add the exact approved configuration**

Append this section to `config/config.yaml` after `is_change` and before `tad`:

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

- [ ] **Step 5: Resolve the ADR-0001 compatibility boundary**

Create `docs/adr/0002-downstream-unified-tad-boundary-view.md` stating:

```markdown
# ADR-0002: Build a downstream unified TAD-boundary view for IS annotation

## Status

Accepted 2026-09-03.

## Decision

Keep the raw and balanced 3DNetMod TAD and multi-resolution boundary outputs
separate as required by ADR-0001. Add a downstream, provenance-preserving union
view used only for IS-change annotation. Overlapping raw/balanced boundaries
may become `both`; unmatched boundaries remain `raw_only` or `balanced_only`.
The union view never overwrites or aliases either source output.

## Consequences

Scientific users can audit normalization sensitivity while the comparison
stage receives one nonduplicated boundary collection per sample. Downstream
evidence tiers retain signal and resolution support, and complex overlap
components remain explicit rather than being forced into one boundary.
```

Add these terms to `CONTEXT.md`: `Unified boundary region`, `Boundary position
state`, `Shifted boundary`, and `IS boundary strength`. Define them with the
same vocabulary as the approved spec and state that the unified view is
downstream-only.

- [ ] **Step 6: Run the config test and existing config suites**

Run:

```bash
python -m pytest -q test/is_boundary/test_config.py test/is_change/test_is_change_config.py test/tad/test_tad_config.py
```

Expected: all selected tests pass with exit code 0.

- [ ] **Step 7: Commit the configuration and decision contract**

```bash
git add config/config.yaml CONTEXT.md docs/adr/0002-downstream-unified-tad-boundary-view.md scripts/is_boundary_errors.py scripts/is_boundary_config.py test/is_boundary/test_config.py
git commit -m "feat: define IS boundary annotation contract"
```

---

### Task 2: Boundary I/O and Per-sample Raw/Balanced Unification

**Files:**

- Create: `scripts/is_boundary_io.py`
- Create: `scripts/is_boundary_unification.py`
- Create: `scripts/unify_tad_boundaries.py`
- Create: `test/is_boundary/test_unification.py`

**Interfaces:**

- Consumes: existing BED7 multi-resolution files and BED6 per-resolution files.
- Produces: `read_multires_boundaries()`, `read_resolution_boundaries()`,
  `unify_boundaries()`, `build_jitter_offsets()`,
  `build_sensitivity_views()`, and a CLI that writes one unified TSV, one
  technical-offset TSV, and one sensitivity-view TSV per sample.

- [ ] **Step 1: Write failing parser and one-to-one union tests**

Create fixtures with exact BED rows and assert deterministic provenance:

```python
from pathlib import Path
import sys

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from is_boundary_io import read_multires_boundaries
from is_boundary_unification import unify_boundaries


def write_bed7(path, rows):
    path.write_text("".join("\t".join(map(str, row)) + "\n" for row in rows))
    return path


def test_reads_bed7_and_unifies_overlap_with_provenance(tmp_path):
    raw_path = write_bed7(tmp_path / "S.multires.raw.boundaries.bed", [
        ("A02", 0, 12000, "raw_1", "4000,8000", 2, 3),
        ("A02", 40000, 48000, "raw_2", "4000", 1, 1),
    ])
    balanced_path = write_bed7(tmp_path / "S.multires.balanced.boundaries.bed", [
        ("A02", 4000, 16000, "bal_1", "4000,20000", 2, 2),
    ])
    raw = read_multires_boundaries(raw_path, "raw")
    balanced = read_multires_boundaries(balanced_path, "balanced")

    observed = unify_boundaries(raw, balanced, tie_tolerance_bp=4000)

    assert observed[[
        "chrom", "start", "end", "signal_support", "support_resolutions",
        "n_source_boundaries", "component_state",
    ]].to_dict("records") == [
        {
            "chrom": "A02", "start": 0, "end": 16000,
            "signal_support": "both", "support_resolutions": "4000,8000,20000",
            "n_source_boundaries": 5, "component_state": "simple",
        },
        {
            "chrom": "A02", "start": 40000, "end": 48000,
            "signal_support": "raw_only", "support_resolutions": "4000",
            "n_source_boundaries": 1, "component_state": "simple",
        },
    ]
```

- [ ] **Step 2: Run the focused test and confirm the missing-module failure**

Run:

```bash
python -m pytest -q test/is_boundary/test_unification.py::test_reads_bed7_and_unifies_overlap_with_provenance
```

Expected: collection fails because `is_boundary_io` does not exist.

- [ ] **Step 3: Implement strict BED parsing and atomic TSV output**

In `scripts/is_boundary_io.py`, define exact schemas and reject malformed,
unsorted, duplicate-ID, nonpositive, or non-half-open coordinates:

```python
import os
from pathlib import Path
import tempfile

import pandas as pd

from is_boundary_errors import ISBoundaryError

MULTIRES_INPUT_COLUMNS = (
    "chrom", "start", "end", "boundary_id", "support_resolutions",
    "n_supporting_resolutions", "n_source_boundaries",
)
UNIFIED_BOUNDARY_COLUMNS = (
    "chrom", "start", "end", "unified_boundary_id", "center_bp",
    "component_state", "signal_support", "raw_boundary_ids",
    "balanced_boundary_ids", "support_resolutions",
    "n_supporting_resolutions", "n_source_boundaries",
)


def _resolution_string(value):
    values = tuple(sorted({int(item) for item in str(value).split(",") if item}))
    if not values or any(item <= 0 for item in values):
        raise ValueError("resolution list must contain positive integers")
    return ",".join(map(str, values))


def read_multires_boundaries(path, signal: str) -> pd.DataFrame:
    if signal not in {"raw", "balanced"}:
        raise ValueError("signal must be raw or balanced")
    frame = pd.read_csv(path, sep="\t", header=None, names=MULTIRES_INPUT_COLUMNS)
    if frame.empty:
        return frame.assign(signal=pd.Series(dtype=str))
    frame["start"] = pd.to_numeric(frame["start"], errors="raise").astype(int)
    frame["end"] = pd.to_numeric(frame["end"], errors="raise").astype(int)
    if (frame["start"] < 0).any() or (frame["end"] <= frame["start"]).any():
        raise ISBoundaryError("IS_BOUNDARY_INPUT_SCHEMA_INVALID", path=path)
    if frame["boundary_id"].duplicated().any():
        raise ISBoundaryError(
            "IS_BOUNDARY_INPUT_SCHEMA_INVALID", path=path, reason="duplicate_boundary_id"
        )
    frame["support_resolutions"] = frame["support_resolutions"].map(_resolution_string)
    frame["signal"] = signal
    chrom_order = {
        chrom: index for index, chrom in enumerate(dict.fromkeys(frame["chrom"]))
    }
    ordered = (
        frame.assign(_chrom_rank=frame["chrom"].map(chrom_order))
        .sort_values(["_chrom_rank", "start", "end", "boundary_id"], kind="stable")
        .drop(columns="_chrom_rank")
    )
    if ordered.index.tolist() != frame.index.tolist():
        raise ISBoundaryError(
            "IS_BOUNDARY_INPUT_SCHEMA_INVALID", path=path, reason="unsorted_rows"
        )
    return frame.reset_index(drop=True)


def atomic_write_tsv(frame: pd.DataFrame, path, columns) -> None:
    if tuple(frame.columns) != tuple(columns):
        raise ISBoundaryError(
            "IS_BOUNDARY_INPUT_SCHEMA_INVALID", path=path, reason="output_schema"
        )
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=target.parent, prefix=f".{target.name}.", suffix=".tmp"
    )
    os.close(descriptor)
    try:
        frame.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, target)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise
```

Add `read_resolution_boundaries(paths, signal)` in the same module. It accepts
a sequence of `(resolution, path)` pairs, validates BED6 rows, and returns the
same normalized coordinate fields plus integer `resolution` and `signal`.

- [ ] **Step 4: Implement overlap assignment and deterministic unified IDs**

In `scripts/is_boundary_unification.py`, use half-open overlap and an ordered
maximum-cardinality assignment. The comparison key is `(number of matches,
negative total center distance)` so distance cannot reduce the number of
eligible overlap matches:

```python
import pandas as pd

from is_boundary_io import UNIFIED_BOUNDARY_COLUMNS


def overlaps(left, right):
    return max(int(left.start), int(right.start)) < min(int(left.end), int(right.end))


def maximum_cardinality_pairs(left, right, eligible):
    n, m = len(left), len(right)
    score = [[(0, 0) for _ in range(m + 1)] for _ in range(n + 1)]
    trace = [[None for _ in range(m + 1)] for _ in range(n + 1)]
    for i in range(n + 1):
        for j in range(m + 1):
            if i < n and score[i][j] > score[i + 1][j]:
                score[i + 1][j] = score[i][j]
                trace[i + 1][j] = (i, j, "skip_left")
            if j < m and score[i][j] > score[i][j + 1]:
                score[i][j + 1] = score[i][j]
                trace[i][j + 1] = (i, j, "skip_right")
            if i < n and j < m and eligible(left.iloc[i], right.iloc[j]):
                distance = abs(int(left.iloc[i].center_bp) - int(right.iloc[j].center_bp))
                candidate = (score[i][j][0] + 1, score[i][j][1] - distance)
                if candidate > score[i + 1][j + 1]:
                    score[i + 1][j + 1] = candidate
                    trace[i + 1][j + 1] = (i, j, "match")
    pairs = []
    i, j = n, m
    while i or j:
        previous = trace[i][j]
        if previous is None:
            if i:
                i -= 1
            else:
                j -= 1
            continue
        previous_i, previous_j, action = previous
        if action == "match":
            pairs.append((previous_i, previous_j))
        i, j = previous_i, previous_j
    return list(reversed(pairs))
```

Before matching, add `center_bp=(start+end)//2` to both frames. For every
chromosome, detect overlap components. Emit a component union as
`complex_within_sample` when it contains multiple members on both sides, or a
matched member's closest and second-closest center distances differ by less
than `tie_tolerance_bp`. Otherwise merge assigned pairs, retain unmatched rows,
and use `start=min(source starts)`, `end=max(source ends)`, and
`center_bp=(start+end)//2`. Union resolution IDs, sum source-boundary counts,
sort using the mcool chromosome order supplied by the CLI, and assign
`unified_boundary_000001`-style IDs in that order. The comparison CLI verifies
that the relative order of chromosomes shared by A and B is identical before
matching.

- [ ] **Step 5: Write and pass complex-component and adjacency tests**

Add:

```python
def test_adjacency_does_not_unify_and_many_to_many_is_complex(tmp_path):
    raw = pd.DataFrame([
        {"chrom": "A02", "start": 0, "end": 8000, "boundary_id": "r1", "support_resolutions": "4000", "n_supporting_resolutions": 1, "n_source_boundaries": 1, "signal": "raw"},
        {"chrom": "A02", "start": 6000, "end": 14000, "boundary_id": "r2", "support_resolutions": "4000", "n_supporting_resolutions": 1, "n_source_boundaries": 1, "signal": "raw"},
        {"chrom": "A02", "start": 30000, "end": 38000, "boundary_id": "r3", "support_resolutions": "4000", "n_supporting_resolutions": 1, "n_source_boundaries": 1, "signal": "raw"},
    ])
    balanced = pd.DataFrame([
        {"chrom": "A02", "start": 4000, "end": 10000, "boundary_id": "b1", "support_resolutions": "4000", "n_supporting_resolutions": 1, "n_source_boundaries": 1, "signal": "balanced"},
        {"chrom": "A02", "start": 8000, "end": 16000, "boundary_id": "b2", "support_resolutions": "4000", "n_supporting_resolutions": 1, "n_source_boundaries": 1, "signal": "balanced"},
        {"chrom": "A02", "start": 38000, "end": 46000, "boundary_id": "b3", "support_resolutions": "4000", "n_supporting_resolutions": 1, "n_source_boundaries": 1, "signal": "balanced"},
    ])

    observed = unify_boundaries(raw, balanced, tie_tolerance_bp=4000)

    assert (observed["component_state"] == "complex_within_sample").sum() == 1
    assert set(observed.loc[observed.start >= 30000, "signal_support"]) == {"raw_only", "balanced_only"}
```

Run:

```bash
python -m pytest -q test/is_boundary/test_unification.py
```

Expected: all unification tests pass.

- [ ] **Step 6: Implement technical offsets and sensitivity views**

`build_jitter_offsets()` emits:

```text
sample chrom source comparison resolution_a resolution_b
boundary_a_id boundary_b_id center_a center_b offset_bp eligible
```

Include raw-versus-balanced one-to-one overlapping pairs, reciprocal
full-versus-leave-one-resolution-out pairs, and reciprocal pairs among the
4/8/20 kb resolution-specific sets. Mark complex components, nonreciprocal
ties, intervals touching a known chromosome end, nonfinite 4 kb cooler bins,
and intervals overlapping `TAD.bad_region_bed` as `eligible=false`; only
`eligible=true` offsets enter the Q95 calculation. The CLI reads chromosome
lengths and 4 kb bin weights from the sample mcool rather than inferring a
chromosome end from the boundary table.

Add `build_sensitivity_views(sample, raw, balanced, per_resolution, settings)`.
It emits `BOUNDARY_VIEW_COLUMNS` for exactly these view names:

```text
unified
raw_only_input
balanced_only_input
leave_r4000_out
leave_r8000_out
leave_r20000_out
```

`raw_only_input` and `balanced_only_input` retain the corresponding original
multi-resolution coordinates and mark the absent signal-ID column as empty;
they are not filters of the already merged unified view. For each leave-one-out
view, remove that resolution from the BED6 input rows, rebuild each signal's
multi-resolution intervals using the same strict half-open transitive-overlap
semantics as `tad_boundaries.merge_intervals()`, then call
`unify_boundaries()`. Add a regression fixture proving that a boundary
supported only at 4 kb disappears from `leave_r4000_out` but remains in the
other two leave-one-out views.

- [ ] **Step 7: Implement the per-sample CLI**

Create `scripts/unify_tad_boundaries.py` with explicit arguments:

```text
--sample
--raw-multires
--balanced-multires
--raw-resolution-boundaries PATH PATH PATH
--balanced-resolution-boundaries PATH PATH PATH
--mcool
--config
--output-unified
--output-offsets
--output-views
--log
```

Map the three resolution paths against the sorted configured TAD resolutions,
call the pure functions, assert exact output schemas, and publish all three
tables atomically. Catch `ISBoundaryError` only to print its stable message to stderr
and return exit code 1; unexpected exceptions are wrapped as
`IS_BOUNDARY_ANNOTATION_FAILED` with `stage=unification`.

- [ ] **Step 8: Add a subprocess CLI test and run the task suite**

The subprocess test supplies three small BED6 files per signal, runs the CLI,
and asserts that all three output tables have headers, stable IDs, exact view
names, and no temporary files. Run:

```bash
python -m pytest -q test/is_boundary/test_unification.py
```

Expected: all tests pass with exit code 0.

- [ ] **Step 9: Commit the per-sample unification stage**

```bash
git add scripts/is_boundary_io.py scripts/is_boundary_unification.py scripts/unify_tad_boundaries.py test/is_boundary/test_unification.py
git commit -m "feat: unify raw and balanced TAD boundaries"
```

---

### Task 3: Jitter Calibration and Ordered A/B Boundary Events

**Files:**

- Create: `scripts/is_boundary_matching.py`
- Create: `scripts/call_tad_boundary_events.py`
- Create: `test/is_boundary/test_matching.py`
- Modify: `scripts/is_boundary_io.py`

**Interfaces:**

- Consumes: two unified-boundary TSVs, two technical-offset TSVs, and two
  sensitivity-view TSVs.
- Produces: `JitterEstimate`, `align_ordered()`, `call_boundary_events()`,
  `compare_event_views()`, a genome-wide boundary-event table, a long
  stability table, and matching QC.

- [ ] **Step 1: Write failing jitter tests**

```python
import pandas as pd

from is_boundary_matching import estimate_jitter


def test_jitter_uses_q95_and_falls_back_below_minimum(settings):
    enough = pd.DataFrame({
        "offset_bp": list(range(0, 400000, 4000)),
        "eligible": True,
    })
    estimate = estimate_jitter(enough, enough, settings)
    assert estimate.source == "empirical_q95"
    assert estimate.n_offsets == 200
    assert estimate.threshold_bp % 4000 == 0
    assert estimate.threshold_bp >= 4000

    too_few = enough.head(10)
    fallback = estimate_jitter(too_few, too_few, settings)
    assert fallback.threshold_bp == 20000
    assert fallback.source == "fallback_max_tad_resolution"
```

The fixture `settings` is parsed from the repository configuration.

- [ ] **Step 2: Run the jitter test and confirm the import failure**

Run:

```bash
python -m pytest -q test/is_boundary/test_matching.py::test_jitter_uses_q95_and_falls_back_below_minimum
```

Expected: collection fails because `is_boundary_matching` does not exist.

- [ ] **Step 3: Implement grid-aligned empirical jitter**

```python
from dataclasses import dataclass
import math

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class JitterEstimate:
    threshold_bp: int
    source: str
    n_offsets: int
    raw_quantile_bp: float


def estimate_jitter(offsets_a, offsets_b, settings) -> JitterEstimate:
    values = pd.concat(
        [offsets_a.loc[offsets_a.eligible, "offset_bp"], offsets_b.loc[offsets_b.eligible, "offset_bp"]],
        ignore_index=True,
    ).dropna().astype(float)
    if len(values) < settings.jitter_min_pairs:
        return JitterEstimate(
            settings.jitter_fallback_bp,
            "fallback_max_tad_resolution",
            len(values),
            float("nan"),
        )
    raw = float(np.quantile(values, settings.jitter_quantile))
    aligned = int(math.ceil(raw / settings.base_bin_bp) * settings.base_bin_bp)
    return JitterEstimate(
        max(settings.base_bin_bp, aligned),
        "empirical_q95",
        len(values),
        raw,
    )
```

- [ ] **Step 4: Write failing ordered-alignment and event-state tests**

Construct two chromosomes with overlapping, jitter-close, shifted, A-only, and
B-only boundaries. Assert exact results and noncrossing pairs:

```python
from is_boundary_io import UNIFIED_BOUNDARY_COLUMNS
from is_boundary_matching import JitterEstimate, call_boundary_events


def unified_frame(rows):
    records = []
    for chrom, start, end, boundary_id in rows:
        records.append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "unified_boundary_id": boundary_id,
            "center_bp": (start + end) // 2,
            "component_state": "simple",
            "signal_support": "both",
            "raw_boundary_ids": f"raw_{boundary_id}",
            "balanced_boundary_ids": f"balanced_{boundary_id}",
            "support_resolutions": "4000,8000",
            "n_supporting_resolutions": 2,
            "n_source_boundaries": 2,
        })
    return pd.DataFrame(records, columns=UNIFIED_BOUNDARY_COLUMNS)


def test_calls_fixed_shift_gain_and_loss_without_crossing(settings):
    a = unified_frame([
        ("A02", 0, 8000, "a1"),
        ("A02", 40000, 48000, "a2"),
        ("A02", 200000, 208000, "a3"),
        ("A02", 500000, 508000, "a4"),
    ])
    b = unified_frame([
        ("A02", 4000, 12000, "b1"),
        ("A02", 52000, 60000, "b2"),
        ("A02", 272000, 280000, "b3"),
        ("A02", 700000, 708000, "b4"),
    ])
    jitter = JitterEstimate(20000, "empirical_q95", 200, 17200.0)

    events = call_boundary_events(
        a, b, jitter, settings, "H1398_vs_L134", "H1398", "L134"
    )

    assert events["position_state"].tolist() == [
        "conserved", "conserved", "shifted", "lost_in_B", "gained_in_B"
    ]
    shifted = events.loc[events.position_state == "shifted"].iloc[0]
    assert shifted.boundary_a_id == "a3"
    assert shifted.boundary_b_id == "b3"
    assert shifted.shift_bp == 72000
```

- [ ] **Step 5: Implement banded dynamic programming**

Use a scalar cost matrix with deterministic tie priority `match`, `delete`,
then `insert`:

```python
def align_ordered(boundaries_a, boundaries_b, max_distance_bp: int):
    n, m = len(boundaries_a), len(boundaries_b)
    cost = np.full((n + 1, m + 1), np.inf)
    trace = np.full((n + 1, m + 1), "", dtype=object)
    cost[:, 0] = np.arange(n + 1) * 0.5
    cost[0, :] = np.arange(m + 1) * 0.5
    trace[1:, 0] = "delete"
    trace[0, 1:] = "insert"
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            distance = abs(
                int(boundaries_a.iloc[i - 1].center_bp)
                - int(boundaries_b.iloc[j - 1].center_bp)
            )
            candidates = [
                (cost[i - 1, j] + 0.5, 1, "delete"),
                (cost[i, j - 1] + 0.5, 2, "insert"),
            ]
            if distance <= max_distance_bp:
                candidates.append(
                    (cost[i - 1, j - 1] + distance / max_distance_bp, 0, "match")
                )
            best_cost, _, action = min(candidates)
            cost[i, j] = best_cost
            trace[i, j] = action
    pairs = []
    i, j = n, m
    while i or j:
        action = trace[i, j]
        if action == "match":
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif action == "delete":
            i -= 1
        else:
            j -= 1
    return list(reversed(pairs))
```

Run this per chromosome only. Stage 1 first assigns interval overlaps and then
mutually closest pairs within `J`. An assigned interval overlap is
`position_state=conserved, position_substate=conserved_exact`; a nonoverlapping
mutually closest pair within `J` is
`position_state=conserved, position_substate=conserved_within_jitter`. Stage 2
calls `align_ordered()` on the remaining boundaries with 100 kb. Assign
unmatched A as `lost_in_B` and unmatched B as `gained_in_B`.
For every matched pair, define signed `shift_bp=center_b-center_a`; unmatched
events store a missing shift rather than zero.

- [ ] **Step 6: Implement uniqueness, complex, and ambiguity rules**

For each candidate pair, calculate closest and second-closest distances on
both sides. Set:

```python
is_reciprocal = closest_b_for_a == b_index and closest_a_for_b == a_index
second_best_margin_bp = min(second_a - first_a, second_b - first_b)
is_tied = second_best_margin_bp < settings.match_tie_tolerance_bp
```

If a local candidate component contains multiple A and B boundaries and is
tied, emit one `complex` component event. A nonreciprocal or tied one-to-one
pair outside a many-to-many component is `ambiguous`. If
`jitter.threshold_bp >= settings.max_shift_bp`, convert otherwise eligible
shifted pairs to `ambiguous` with
`ambiguity_reason=jitter_not_smaller_than_shift_radius`.

Populate the exact `BOUNDARY_EVENT_COLUMNS` from the design, including separate
`a_start/a_end` and `b_start/b_end`. `event_span_start/end` are diagnostic
columns only.

- [ ] **Step 7: Add sensitivity and direction-swap tests**

Call events at 60, 80, and 100 kb and record one column per radius:

```text
state_at_60000
state_at_80000
state_at_100000
position_state_stable
```

Add a test that swaps A/B and verifies gain/loss reversal, conserved and
shifted preservation, equal absolute shift, and opposite signed shift.

Add `compare_event_views()`. For each shared non-`unified` view, call the same
event classifier and associate each reference event to sensitivity events by
its A and B endpoint intervals, never by the shifted corridor. Zero associated
events means `sensitivity_position_state=no_boundary_in_view`; more than one
or a complex association means `complex`; one means the associated state.
Emit one row per reference event and view using `BOUNDARY_STABILITY_COLUMNS`.
Set `position_state_stable=True` in the reference event table only when its
state agrees across all three shift-radius columns and every applicable source
and leave-one-resolution-out view. Unit tests must cover one stable event and
one event lost only after removing 4 kb support.

Run:

```bash
python -m pytest -q test/is_boundary/test_matching.py
```

Expected: all matching tests pass with exit code 0.

- [ ] **Step 8: Implement the boundary-event CLI and atomic QC output**

Create `scripts/call_tad_boundary_events.py` with:

```text
--comparison
--sample-a
--sample-b
--boundaries-a
--boundaries-b
--offsets-a
--offsets-b
--views-a
--views-b
--config
--output-events
--output-stability
--output-qc
--log
```

The QC table uses exact columns `metric,value,detail` and contains
`jitter_threshold_bp`, `jitter_source`, `n_offsets`,
`raw_quantile_bp`, counts by position state, shift-radius concordance,
raw-only/balanced-only concordance, and leave-one-resolution-out concordance.
Reject a comparison if the supplied sample names equal one another or the
event input chromosomes are not present on the same coordinate grid, using
`IS_BOUNDARY_REFERENCE_MISMATCH`. Wrap an unrecoverable calibration error with
`IS_BOUNDARY_JITTER_CALIBRATION_FAILED` and an event-construction error with
`IS_BOUNDARY_MATCH_FAILED`; insufficient but nonzero calibration support uses
the fallback and QC instead of raising.

- [ ] **Step 9: Run module and CLI tests, then commit**

Run:

```bash
python -m pytest -q test/is_boundary/test_matching.py
```

Expected: all tests pass.

Commit:

```bash
git add scripts/is_boundary_io.py scripts/is_boundary_matching.py scripts/call_tad_boundary_events.py test/is_boundary/test_matching.py
git commit -m "feat: classify ordered TAD boundary events"
```

---

### Task 4: Multi-scale IS Valley Strength

**Files:**

- Create: `scripts/is_boundary_strength.py`
- Create: `test/is_boundary/test_strength.py`

**Interfaces:**

- Consumes: boundary events, validated balanced IS tracks from both samples, `ISBoundarySettings`, and `ISChangeSettings`.
- Produces: a per-event/per-coordinate-role/per-scale long table and an
  event-level strength summary using `STRENGTH_SUMMARY_COLUMNS`.

- [ ] **Step 1: Write failing valley-prominence tests**

Use a 4 kb track with a clear valley, a fully valid flat track, and a track
with an invalid flank:

```python
import numpy as np
import pandas as pd

from is_boundary_strength import local_valley_strength


def test_local_valley_strength_distinguishes_valley_flat_and_missing():
    coordinates = np.arange(9) * 4000
    valley = pd.DataFrame({
        "start": coordinates,
        "IS": [1.0, 0.8, 0.5, -2.0, 0.4, 0.7, 1.0, 0.9, 0.8],
    })
    flat = valley.assign(IS=0.0)
    missing = valley.copy()
    missing.loc[0:3, "IS"] = np.nan

    strong = local_valley_strength(valley, 12000, 20000, 32000, 0.80)
    no_valley = local_valley_strength(flat, 12000, 20000, 32000, 0.80)
    invalid = local_valley_strength(missing, 12000, 20000, 32000, 0.80)

    assert strong.strength > 0
    assert strong.detection == "valley"
    assert no_valley.strength == 0
    assert no_valley.detection == "no_valley"
    assert np.isnan(invalid.strength)
    assert invalid.detection == "insufficient_coverage"
```

- [ ] **Step 2: Run the focused test and confirm the missing-module failure**

Run:

```bash
python -m pytest -q test/is_boundary/test_strength.py::test_local_valley_strength_distinguishes_valley_flat_and_missing
```

Expected: collection fails because `is_boundary_strength` does not exist.

- [ ] **Step 3: Implement finite-run, plateau-aware local prominence**

```python
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_prominences


@dataclass(frozen=True)
class StrengthObservation:
    strength: float
    valley_start: float
    detection: str
    valid_fraction: float


def _finite_runs(mask):
    indices = np.flatnonzero(mask)
    if not len(indices):
        return []
    breaks = np.flatnonzero(np.diff(indices) > 1) + 1
    return [part for part in np.split(indices, breaks) if len(part)]


def local_valley_strength(track, start, end, scale_bp, min_valid_fraction):
    radius = int(scale_bp)
    segment = track.loc[
        track.start.between(start - radius, end + radius - 1), ["start", "IS"]
    ].sort_values("start", kind="stable")
    finite = np.isfinite(segment.IS.to_numpy(dtype=float))
    valid_fraction = float(finite.mean()) if len(finite) else 0.0
    if not len(segment) or valid_fraction < min_valid_fraction:
        return StrengthObservation(float("nan"), float("nan"), "insufficient_coverage", valid_fraction)
    candidates = []
    values = segment.IS.to_numpy(dtype=float)
    starts = segment.start.to_numpy(dtype=int)
    for run in _finite_runs(finite):
        run_values = values[run]
        peaks, _ = find_peaks(-run_values, plateau_size=True)
        if not len(peaks):
            continue
        prominences = peak_prominences(-run_values, peaks)[0]
        for local_peak, prominence in zip(peaks, prominences):
            global_index = run[local_peak]
            coordinate = starts[global_index]
            if start <= coordinate < end:
                candidates.append((float(prominence), int(coordinate)))
    if not candidates:
        return StrengthObservation(0.0, float("nan"), "no_valley", valid_fraction)
    prominence, coordinate = max(candidates, key=lambda item: (item[0], -item[1]))
    return StrengthObservation(prominence, coordinate, "valley", valid_fraction)
```

Before accepting a candidate, separately require left and right scale-radius
flanks to meet `is_settings.min_flank_valid_fraction`. Preserve `NA` when either
flank fails.

- [ ] **Step 4: Write failing robust-Z and two-band classification tests**

Build 60 conserved events across all ten scales so the reference count exceeds
50. Inject a positive effect into small and medium bands for one target event
and assert:

```python
summary = summarize_strength(strength_long, events, settings, is_settings)
target = summary.loc[summary.boundary_event_id == "event_target"].iloc[0]
assert target.strength_state == "strengthened_in_B"
assert target.strength_support_bands == "small,medium"
```

Also assert negative effects become `weakened_in_B`, all band scores below one
become `stable_strength`, and zero MAD makes the affected scale unevaluable.
Run the same fixture at robust-Z thresholds 1.5, 2.0, and 2.5; assert that the
summary records `state_at_z1_5`, `state_at_z2_0`, `state_at_z2_5`, and
`strength_state_stable` without overwriting the primary 2.0 call.

- [ ] **Step 5: Implement coordinate roles and robust standardization**

`compute_strength_long()` emits roles exactly as follows:

```python
ROLE_BY_STATE = {
    "conserved": ("shared",),
    "gained_in_B": ("new_b",),
    "lost_in_B": ("old_a",),
    "shifted": ("old_a", "new_b"),
    "complex": ("diagnostic",),
    "ambiguous": ("diagnostic",),
}
```

For each role and scale, evaluate both tracks at the role's boundary interval
and write `strength_a`, `strength_b`, `delta_strength`, valley coordinates,
detection states, and valid fractions. Calculate reference median/MAD from
`position_state=conserved`, `coordinate_role=shared`, and finite deltas only.
Require `settings.strength.min_reference_boundaries` values and nonzero MAD.

Band scores are medians of valid robust-Z values and require
`min_valid_scales_per_band`. At the primary threshold, at least two band scores
at or above +2 yield `strengthened_in_B`; at least two at or below -2 yield
`weakened_in_B`; at least two evaluable bands with every absolute score below
1 yield `stable_strength`; all remaining cases are
`discordant_or_unevaluable`.
The event-level summary contains the primary `strength_state`, the three
threshold-sensitivity state columns, and a Boolean `strength_state_stable`.
For gained and lost events, the expected sign must cross the absolute primary
threshold in at least two bands to be `direction_concordant`. For shifted
events, require old-A weakening and new-B strengthening separately: both
endpoints supported gives `bidirectional_shift_support`, exactly one gives
`partial_shift_support`, an opposite or absent effect gives `discordant`, and
insufficient valid scales gives `unevaluable`. Thus rearrangements populate:

```text
direction_concordant
bidirectional_shift_support
partial_shift_support
discordant
unevaluable
```

- [ ] **Step 6: Run the strength suite and commit**

Run:

```bash
python -m pytest -q test/is_boundary/test_strength.py
```

Expected: all strength tests pass.

Commit:

```bash
git add scripts/is_boundary_strength.py test/is_boundary/test_strength.py
git commit -m "feat: score multiscale IS boundary strength"
```

---

### Task 5: IS-locus Association, Evidence Tier, and Annotation CLI

**Files:**

- Create: `scripts/is_boundary_annotation.py`
- Create: `scripts/annotate_is_change_boundaries.py`
- Create: `test/is_boundary/test_annotation.py`
- Modify: `scripts/is_boundary_io.py`

**Interfaces:**

- Consumes: primary/preliminary IS tables, event map, both IS tracks, and the existing IS mask.
- Produces: annotated primary/preliminary TSVs, per-scale strength TSV,
  event-level strength-summary TSV, annotation QC, and annotation summary.

- [ ] **Step 1: Write failing endpoint-only association tests**

```python
import pandas as pd

from is_boundary_annotation import associated_event_ids


def test_shifted_span_is_not_an_overlap_target():
    event = pd.Series({
        "boundary_event_id": "e1",
        "position_state": "shifted",
        "a_start": 100000,
        "a_end": 108000,
        "b_start": 180000,
        "b_end": 188000,
        "event_span_start": 100000,
        "event_span_end": 188000,
    })
    between = pd.Series({
        "is_change_type": "transition", "summit": 140000,
        "event_start": 140000, "event_end": 144000,
    })
    old_endpoint = between.copy()
    old_endpoint["summit"] = 104000
    old_endpoint["event_start"] = 104000
    old_endpoint["event_end"] = 108000

    assert associated_event_ids(between, pd.DataFrame([event]), 4000) == []
    assert associated_event_ids(old_endpoint, pd.DataFrame([event]), 4000) == ["e1"]
```

- [ ] **Step 2: Run the focused test and confirm the missing-module failure**

Run:

```bash
python -m pytest -q test/is_boundary/test_annotation.py::test_shifted_span_is_not_an_overlap_target
```

Expected: collection fails because `is_boundary_annotation` does not exist.

- [ ] **Step 3: Implement anchor creation and candidate association**

```python
SIMPLE_POSITION_STATES = {"conserved", "gained_in_B", "lost_in_B", "shifted"}


def locus_anchor(locus, base_bin_bp):
    if locus.is_change_type == "broad_shift":
        return int(locus.event_start), int(locus.event_end)
    return int(locus.summit), int(locus.summit) + base_bin_bp


def _interval_overlap(start_a, end_a, start_b, end_b):
    return max(start_a, start_b) < min(end_a, end_b)


def associated_event_ids(locus, events, base_bin_bp):
    anchor_start, anchor_end = locus_anchor(locus, base_bin_bp)
    identifiers = []
    for event in events.itertuples(index=False):
        supports = []
        if pd.notna(event.a_start) and pd.notna(event.a_end):
            supports.append((int(event.a_start), int(event.a_end)))
        if pd.notna(event.b_start) and pd.notna(event.b_end):
            supports.append((int(event.b_start), int(event.b_end)))
        if any(_interval_overlap(anchor_start, anchor_end, start, end) for start, end in supports):
            identifiers.append(event.boundary_event_id)
    return identifiers
```

Filter event candidates to the locus chromosome before calling this function.
Never inspect `event_span_start/end` in association code.
Store zero candidates as an empty `boundary_event_id`, one candidate as its
ID, and multiple candidates as comma-separated IDs in deterministic event
order; `n_boundary_candidates` always stores the integer count.

- [ ] **Step 4: Write failing state, missingness, and evidence-tier tests**

Test these exact outcomes:

```text
zero candidates + clean mask       -> no_boundary_association
zero candidates + overlapping mask -> ambiguous
one simple candidate               -> inherit event state
two distinct simple candidates     -> complex
one complex event                  -> complex
gain/loss + overlapping mask       -> ambiguous
primary + unique + robust support  -> high
primary + one missing support axis -> moderate
primary + position instability    -> moderate
single-signal/resolution support  -> low
fallback jitter                   -> low
preliminary                        -> low
complex or ambiguous               -> not_applicable
```

For `conserved`, assert that the independently supplied strength summary can
produce `strength_state=weakened_in_B` without changing `position_state`.

- [ ] **Step 5: Implement annotation and evidence tiers**

`annotate_loci()` preserves all original input columns in original order, then
adds the exact `ANNOTATION_COLUMNS`. It joins strength only by
`boundary_event_id` when `n_boundary_candidates == 1`; multi-candidate loci
receive `strength_state=not_applicable` and retain the comma-separated IDs only
as provenance. Apply this tier function:

```python
def resolution_count(value):
    if pd.isna(value) or str(value).strip() == "":
        return 0
    return len({int(item) for item in str(value).split(",")})


def evidence_tier(row):
    if row.position_state in {"complex", "ambiguous"}:
        return "not_applicable"
    if row.confidence_tier == "single_band_preliminary":
        return "low"
    if row.jitter_source != "empirical_q95":
        return "low"
    support_a = row.signal_support_a == "both" or resolution_count(row.resolution_support_a) >= 2
    support_b = row.signal_support_b == "both" or resolution_count(row.resolution_support_b) >= 2
    present_support = {
        "conserved": support_a and support_b,
        "shifted": support_a and support_b,
        "lost_in_B": support_a,
        "gained_in_B": support_b,
        "no_boundary_association": True,
    }[row.position_state]
    if not present_support:
        return "low"
    direction_ok = row.is_consistency_state in {
        "direction_concordant", "bidirectional_shift_support", "not_applicable"
    }
    support_bands = {
        value for value in str(row.strength_support_bands).split(",")
        if value and value != "nan"
    }
    if row.is_consistency_state != "not_applicable" and len(support_bands) < 2:
        return "low"
    strength_stable = (
        row.position_state == "no_boundary_association"
        or (pd.notna(row.strength_state_stable) and bool(row.strength_state_stable))
    )
    if row.position_state_stable and strength_stable and direction_ok:
        return "high"
    return "moderate"
```

For conserved strength changes, treat a two-band strengthened or weakened state
as direction support. A stable conserved boundary uses
`is_consistency_state=not_applicable`.

- [ ] **Step 6: Implement summary tables**

`summarize_annotations()` returns tidy counts and fractions grouped by:

```text
confidence_tier
position_state
strength_state
is_consistency_state
evidence_tier
```

The exact summary columns are
`confidence_tier,position_state,strength_state,is_consistency_state,evidence_tier,n_loci,denominator,fraction`.
The QC table uses `metric,value,detail` and records missing inputs, unevaluable
scales, candidate multiplicity, fallback jitter use, and label changes across
sensitivity thresholds. Sort summary states by explicit categorical order,
not lexical order.

- [ ] **Step 7: Implement the annotation CLI**

Create `scripts/annotate_is_change_boundaries.py` with:

```text
--comparison
--events
--primary
--preliminary
--masked
--is-a
--is-b
--config
--out-primary
--out-preliminary
--out-strength
--out-strength-summary
--out-qc
--out-summary
--log
```

Reuse `read_is_track()` and `index_is_track()` from
`scripts/call_is_change_loci.py`; do not create a third IS parser. Parse both
existing settings objects, calculate strength once for the genome-wide event
map, annotate primary and preliminary separately, verify exact schemas, and
atomically publish all six tables. Wrap failures before strength-summary
publication as `IS_BOUNDARY_STRENGTH_FAILED`; wrap association, evidence-tier,
or annotation-output failures as `IS_BOUNDARY_ANNOTATION_FAILED`.

- [ ] **Step 8: Add subprocess tests for empty and populated inputs**

Verify that header-only primary or preliminary inputs produce valid header-only
annotated tables, and that a populated fixture preserves every original column
before appending annotation columns. Simulate serialization failure and assert
that no `.*.tmp` file remains.

Run:

```bash
python -m pytest -q test/is_boundary/test_annotation.py
```

Expected: all annotation tests pass.

- [ ] **Step 9: Commit the annotation stage**

```bash
git add scripts/is_boundary_io.py scripts/is_boundary_annotation.py scripts/annotate_is_change_boundaries.py test/is_boundary/test_annotation.py
git commit -m "feat: annotate IS loci with boundary events"
```

---

### Task 6: Boundary-aware Diagnostic Figures

**Files:**

- Modify: `scripts/plot_is_change_loci.py:27-390`
- Modify: `test/is_change/test_plot_is_change_loci.py`

**Interfaces:**

- Consumes: annotated primary loci plus the boundary-event map.
- Produces: existing figure manifest and PNGs with sample-specific boundary endpoint overlays and annotation titles.

- [ ] **Step 1: Write a failing overlay test**

Add a test that builds an annotated event with different A/B endpoints and
asserts the sample-specific axes receive different spans:

```python
def test_draw_boundary_support_uses_sample_specific_endpoint():
    figure, axes = plt.subplots(2, 1)
    event = pd.Series({
        "a_start": 100000, "a_end": 108000,
        "b_start": 180000, "b_end": 188000,
        "position_state": "shifted",
    })

    draw_boundary_support(axes[0], event, "a")
    draw_boundary_support(axes[1], event, "b")

    assert axes[0].patches[0].get_xy()[0][0] == 100000
    assert axes[1].patches[0].get_xy()[0][0] == 180000
    plt.close(figure)
```

- [ ] **Step 2: Run the focused test and confirm the missing-symbol failure**

Run:

```bash
python -m pytest -q test/is_change/test_plot_is_change_loci.py::test_draw_boundary_support_uses_sample_specific_endpoint
```

Expected: import fails for `draw_boundary_support`.

- [ ] **Step 3: Add optional event input and endpoint overlays**

Add `--boundary-events` with default `None` so disabling this feature preserves
the existing plotting path. When the argument is present, validate the event
schema and implement:

```python
BOUNDARY_COLORS = {"a": "#2563eb", "b": "#16a34a"}


def draw_boundary_support(axis, event, sample):
    start = event[f"{sample}_start"]
    end = event[f"{sample}_end"]
    if pd.isna(start) or pd.isna(end):
        return
    axis.axvspan(
        int(start), int(end), facecolor="none",
        edgecolor=BOUNDARY_COLORS[sample], linewidth=1.2,
    )
```

Draw A support on the A contact map, B support on the B contact map, and both
outlines on every IS axis. For complex events, draw the diagnostic component
with a gray hatched span. Add a figure title containing `position_state`,
`strength_state`, `is_consistency_state`, and `evidence_tier`.

- [ ] **Step 4: Make plotting consume annotated primary loci**

When `--boundary-events` is present, require the annotation columns and fail if
`confidence_tier` is not uniformly `primary`. When it is absent, keep the
existing unannotated schema and title. Retain current all-scale plotting and
shared contact normalization. Do not alter IS-change event red markers.

- [ ] **Step 5: Run plot and figure tests**

Run:

```bash
python -m pytest -q test/is_change/test_plot_is_change_loci.py
```

Expected: all plot tests pass and temporary figures are cleaned on failure.

- [ ] **Step 6: Commit the visualization update**

```bash
git add scripts/plot_is_change_loci.py test/is_change/test_plot_is_change_loci.py
git commit -m "feat: overlay TAD boundary events on IS plots"
```

---

### Task 7: Optional Constrained Set-level Permutation

**Files:**

- Create: `scripts/is_boundary_permutation.py`
- Create: `scripts/permute_is_boundary_association.py`
- Create: `test/is_boundary/test_permutation.py`

**Interfaces:**

- Consumes: one annotated locus tier, candidate regions, event map, mask, and
  the full `ISBoundarySettings` object so geometry and nested permutation
  settings cannot disagree.
- Produces: one tidy set-level permutation TSV; it does not alter locus annotations.

- [ ] **Step 1: Write failing constrained-sampling tests**

Create two candidate regions on different chromosomes and assert every sampled
anchor remains in its original region, has identical length, stays on the 4 kb
grid, and avoids masked bins.

```python
def test_randomized_loci_stay_in_candidate_region_and_avoid_mask(settings):
    loci = pd.DataFrame([
        {"locus_id": "l1", "region_id": "r1", "chrom": "A02", "event_start": 0, "event_end": 4000, "summit": 0, "is_change_type": "transition", "mask_status": "clean"},
        {"locus_id": "l2", "region_id": "r2", "chrom": "A03", "event_start": 100000, "event_end": 108000, "summit": np.nan, "is_change_type": "broad_shift", "mask_status": "clean"},
    ])
    regions = pd.DataFrame([
        {"region_id": "r1", "chrom": "A02", "region_start": 0, "region_end": 40000},
        {"region_id": "r2", "chrom": "A03", "region_start": 100000, "region_end": 160000},
    ])
    masked = pd.DataFrame([
        {"chrom": "A02", "start": 12000, "end": 16000, "reason": "IS_NA_A"},
        {"chrom": "A02", "start": 16000, "end": 20000, "reason": "IS_NA_B"},
    ])
    randomized = randomize_loci_once(loci, regions, masked, settings, np.random.default_rng(7))
    joined = randomized.merge(regions, on="region_id", suffixes=("", "_region"))
    assert (joined.chrom == joined.chrom_region).all()
    assert (joined.event_start >= joined.region_start).all()
    assert (joined.event_end <= joined.region_end).all()
    assert ((joined.event_start % 4000) == 0).all()
    overlaps = randomized.merge(masked, on="chrom")
    assert not ((overlaps.event_start < overlaps.end) & (overlaps.event_end > overlaps.start)).any()
```

- [ ] **Step 2: Run the focused test and confirm the missing-module failure**

Run:

```bash
python -m pytest -q test/is_boundary/test_permutation.py::test_randomized_loci_stay_in_candidate_region_and_avoid_mask
```

Expected: collection fails because `is_boundary_permutation` does not exist.

- [ ] **Step 3: Implement candidate-region-constrained randomization**

For each locus, enumerate every grid-aligned start that keeps the original
event length inside its own `region_id`. Remove starts whose candidate interval
overlaps the mask. Sample one valid start with replacement from that locus's
eligible list. Raise `IS_BOUNDARY_ANNOTATION_FAILED` with
`stage=permutation` and the locus ID if the eligible list is empty.

```python
def randomize_loci_once(loci, candidate_regions, masked, settings, rng):
    regions = candidate_regions.set_index("region_id", verify_integrity=True)
    randomized = []
    for locus in loci.itertuples(index=False):
        if locus.mask_status != "clean":
            raise ISBoundaryError(
                "IS_BOUNDARY_INPUT_SCHEMA_INVALID",
                locus_id=locus.locus_id,
                reason="permutation_requires_clean_locus",
            )
        region = regions.loc[locus.region_id]
        if str(region.chrom) != str(locus.chrom):
            raise ISBoundaryError(
                "IS_BOUNDARY_INPUT_SCHEMA_INVALID",
                locus_id=locus.locus_id,
                reason="region_chromosome_mismatch",
            )
        length = int(locus.event_end) - int(locus.event_start)
        grid = settings.base_bin_bp
        first = ((int(region.region_start) + grid - 1) // grid) * grid
        last = ((int(region.region_end) - length) // grid) * grid
        candidates = np.arange(first, last + grid, grid, dtype=int)
        local_mask = masked.loc[masked.chrom.astype(str) == str(locus.chrom)]
        eligible = [
            int(start)
            for start in candidates
            if not (
                (local_mask.start < start + length)
                & (local_mask.end > start)
            ).any()
        ]
        if not eligible:
            raise ISBoundaryError(
                "IS_BOUNDARY_ANNOTATION_FAILED",
                stage="permutation",
                locus_id=locus.locus_id,
            )
        record = locus._asdict()
        record["event_start"] = int(rng.choice(eligible))
        record["event_end"] = record["event_start"] + length
        if record["is_change_type"] != "broad_shift":
            record["summit"] = record["event_start"]
        randomized.append(record)
    return pd.DataFrame(randomized, columns=loci.columns)
```

Use the same endpoint-only association function from
`scripts/is_boundary_annotation.py`; do not duplicate shifted-event overlap
logic. Run primary and preliminary strata separately. Since accepted IS-change
loci are required to have `mask_status=clean`, excluding masked candidate
starts preserves their evaluability class as well as chromosome, event length,
geometry, and candidate-region membership.

- [ ] **Step 4: Write failing empirical-P and BH tests**

```python
def test_empirical_p_has_plus_one_correction_and_bh_is_monotone():
    null = np.array([0.1, 0.2, 0.3, 0.4])
    assert empirical_upper_p(0.35, null) == 0.4
    assert empirical_lower_p(0.15, null) == 0.4
    adjusted = bh_adjust(np.array([0.01, 0.04, 0.03]))
    assert np.all((0 <= adjusted) & (adjusted <= 1))
    assert adjusted[0] <= adjusted[2]
```

Implement BH adjustment without adding a statistics dependency:

```python
def bh_adjust(pvalues):
    values = np.asarray(pvalues, dtype=float)
    order = np.argsort(values, kind="stable")
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    output = np.empty_like(ranked)
    output[order] = np.clip(ranked, 0.0, 1.0)
    return output
```

- [ ] **Step 5: Implement pre-specified statistics and output**

Calculate exactly:

```text
any_boundary_fraction
rearranged_given_boundary_fraction
gained_in_B_fraction
lost_in_B_fraction
shifted_fraction
median_nearest_boundary_distance_bp
```

For `rearranged_given_boundary_fraction`, a draw with zero boundary-associated
loci is missing rather than zero. Drop missing null draws for that statistic;
if fewer than 95% of the configured permutations remain, publish missing
inferential fields and `detail=insufficient_valid_null_draws`.

For each statistic emit `observed`, `null_mean`, `null_q025`, `null_q975`,
`fold_enrichment`, `empirical_p`, `bh_qvalue`, `n_permutations`, and `seed`.
The exact output order is
`comparison,confidence_tier,statistic,observed,null_mean,null_q025,null_q975,fold_enrichment,empirical_p,bh_qvalue,n_permutations,seed,detail`.
Set `fold_enrichment` to missing when `null_mean == 0`.
Use:

```python
def empirical_upper_p(observed, null):
    return float((1 + np.count_nonzero(null >= observed)) / (len(null) + 1))


def empirical_lower_p(observed, null):
    return float((1 + np.count_nonzero(null <= observed)) / (len(null) + 1))
```

Use the upper tail for all fraction statistics and the lower tail for
`median_nearest_boundary_distance_bp`, because unusually small distances are
the association signal. Apply BH once across the six pre-specified valid
statistics within each confidence tier.

- [ ] **Step 6: Implement the optional CLI and deterministic-seed test**

Create CLI arguments:

```text
--loci
--candidate-regions
--events
--masked
--config
--output
--log
```

The CLI always runs when its output is explicitly targeted. The
`permutation.enabled` flag controls only whether Snakemake adds the permutation
output to the default target set. Two runs with the same seed and fixtures must
produce byte-identical TSV output.

- [ ] **Step 7: Run and commit the permutation module**

Run:

```bash
python -m pytest -q test/is_boundary/test_permutation.py
```

Expected: all permutation tests pass.

Commit:

```bash
git add scripts/is_boundary_permutation.py scripts/permute_is_boundary_association.py test/is_boundary/test_permutation.py
git commit -m "feat: add optional boundary association permutation"
```

---

### Task 8: Snakemake Wiring and Stable Output Targets

**Files:**

- Create: `rules/is_boundary_annotation.smk`
- Create: `test/is_boundary/test_workflow.py`
- Modify: `Snakefile:9-44,143-269`
- Modify: `rules/is_change.smk:85-126`
- Modify: `test/is_change/test_is_change_workflow.py`

**Interfaces:**

- Consumes: all CLIs created in Tasks 2, 3, 5, and 7.
- Produces: reusable per-sample unified boundaries, per-comparison event maps, annotated primary/preliminary outputs, QC, plots, and optional permutation targets.

- [ ] **Step 1: Write a failing static workflow-contract test**

```python
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def test_boundary_annotation_rules_and_all_outputs_are_wired():
    snakefile = (REPO / "Snakefile").read_text()
    rules = (REPO / "rules" / "is_boundary_annotation.smk").read_text()
    for name in (
        "rule unify_tad_boundaries:",
        "rule call_tad_boundary_events:",
        "rule annotate_is_change_boundaries:",
        "rule permute_is_boundary_association:",
    ):
        assert name in rules
    assert 'include: "rules/tad.smk"' in snakefile
    assert 'include: "rules/is_boundary_annotation.smk"' in snakefile
    assert snakefile.index('include: "rules/tad.smk"') < snakefile.index(
        'include: "rules/is_boundary_annotation.smk"'
    )
    assert "IS_BOUNDARY_ALL_OUTPUTS" in snakefile
    assert "PERMUTATION_ALL_OUTPUTS" in snakefile
```

- [ ] **Step 2: Run the static test and confirm the missing-rule-file failure**

Run:

```bash
python -m pytest -q test/is_boundary/test_workflow.py::test_boundary_annotation_rules_and_all_outputs_are_wired
```

Expected: fails because `rules/is_boundary_annotation.smk` is absent.

- [ ] **Step 3: Parse settings and add path helpers in Snakefile**

Import and parse after `TAD`:

```python
from is_boundary_config import parse_is_boundary_settings

IS_BOUNDARY = parse_is_boundary_settings(config, CWORLD_IS, TAD)
```

Add helpers using `.format()` because the repository's workflow test prohibits
Python f-strings in `Snakefile`:

```python
def unified_boundary_path(sample):
    return "{}/tad_boundaries/{}.multires.unified.boundaries.tsv".format(OUTDIR, sample)


def boundary_jitter_path(sample):
    return "{}/tad_boundaries/{}.boundary_jitter_offsets.tsv".format(OUTDIR, sample)


def boundary_view_path(sample):
    return "{}/tad_boundaries/{}.annotation_views.tsv".format(OUTDIR, sample)


def boundary_event_path(comparison):
    return "{}/tad_boundary_events/{}.boundary_events.tsv".format(OUTDIR, comparison)


def boundary_stability_path(comparison):
    return "{}/tad_boundary_events/{}.boundary_event_stability.tsv".format(OUTDIR, comparison)
```

- [ ] **Step 4: Add the per-sample unification rule**

In `rules/is_boundary_annotation.smk`, define inputs from
`tad_merged_boundary_path()`, ordered `tad_boundary_path()` lists, and
`get_sample_mcool()`. Use one thread and `../env/SNAP-HiC-Diff.yaml`. The
outputs are:

```text
results/tad_boundaries/{sample}.multires.unified.boundaries.tsv
results/tad_boundaries/{sample}.boundary_jitter_offsets.tsv
results/tad_boundaries/{sample}.annotation_views.tsv
```

Pass the exact CLI arguments from Task 2 and write stable log/benchmark files
under `results/logs/is_boundary/` and `results/benchmark/is_boundary/`.

- [ ] **Step 5: Add comparison event and annotation rules**

`call_tad_boundary_events` consumes two unified tables, two offset tables, and
two sensitivity-view tables through `get_pair_samples()`.
`annotate_is_change_boundaries` consumes the event map, primary, preliminary,
mask, and two IS files. Publish every output path from the approved design.

Use separate output names:

```text
results/tad_boundary_events/{comparison}.boundary_events.tsv
results/tad_boundary_events/{comparison}.boundary_event_stability.tsv
results/tad_boundary_events/{comparison}.matching_qc.tsv
results/tad_boundary_events/{comparison}.boundary_strength_by_scale.tsv
results/tad_boundary_events/{comparison}.boundary_strength_summary.tsv
results/tad_boundary_events/{comparison}.annotation_qc.tsv
results/tad_boundary_events/{comparison}.annotation_summary.tsv
results/multiscale/r4000/{comparison}/is_change/{comparison}_is_change_loci.boundary_annotated.tsv
results/multiscale/r4000/{comparison}/is_change/{comparison}_is_change_loci_preliminary.boundary_annotated.tsv
```

- [ ] **Step 6: Wire optional permutation without scheduling it by default**

Define `permute_is_boundary_association` unconditionally so it can be targeted.
Construct:

```python
PERMUTATION_ALL_OUTPUTS = (
    expand(
        OUTDIR + "/tad_boundary_events/{comparison}.permutation.tsv",
        comparison=get_window_compare_pairs(),
    )
    if IS_BOUNDARY.permutation.enabled
    else []
)
```

Do not put the permutation output in `rule all` while `enabled: false`.

- [ ] **Step 7: Make diagnostic plotting depend on annotation**

Modify `rule plot_is_change_loci` so `primary` points to the annotated primary
output and add `boundary_events` as an input when
`IS_BOUNDARY.enabled=true`; otherwise retain the original primary input and
omit the optional CLI argument. Update existing static workflow tests to
exercise both configurations.

- [ ] **Step 8: Add stable outputs to rule all**

Create `IS_BOUNDARY_ALL_OUTPUTS` from all samples and configured comparisons
when `IS_BOUNDARY.enabled` is true, otherwise set it to an empty list. Include
it and `PERMUTATION_ALL_OUTPUTS` in `rule all`. Keep the original raw/balanced
TAD outputs in `TAD_ALL_OUTPUTS` unchanged. Reject
`permutation.enabled=true` when `IS_BOUNDARY.enabled=false` during config
parsing.

- [ ] **Step 9: Write and run a targeted dry-run test**

The test executes:

```bash
snakemake -n -p --configfile config/config.yaml -- results/multiscale/r4000/H1398_vs_L134/is_change/H1398_vs_L134_is_change_loci.boundary_annotated.tsv
```

Expected: exit code 0; the plan contains `unify_tad_boundaries`,
`call_tad_boundary_events`, and `annotate_is_change_boundaries` for H1398 and
L134, but not `permute_is_boundary_association`.

Run:

```bash
python -m pytest -q test/is_boundary/test_workflow.py test/is_change/test_is_change_workflow.py
```

Expected: all workflow contract tests pass.

- [ ] **Step 10: Commit workflow integration**

```bash
git add Snakefile rules/is_boundary_annotation.smk rules/is_change.smk test/is_boundary/test_workflow.py test/is_change/test_is_change_workflow.py
git commit -m "feat: wire IS boundary annotation workflow"
```

---

### Task 9: Synthetic Benchmark, Documentation, and Final Verification

**Files:**

- Create: `test/is_boundary/test_synthetic_benchmark.py`
- Modify: `README.md:20-220`

**Interfaces:**

- Consumes: all public modules and workflow targets from Tasks 1-8.
- Produces: executable truth benchmarks, user-facing execution instructions, and final acceptance evidence.

- [ ] **Step 1: Write the noiseless synthetic truth benchmark**

Create ordered boundaries with at least 20 examples of each simple state.
Ensure unrelated boundary centers are separated by more than 200 kb. Inject
shifts of 40, 60, and 80 kb with `J=20 kb` and assert exact labels for all
noiseless examples.

```python
import numpy as np
import pandas as pd

from is_boundary_io import UNIFIED_BOUNDARY_COLUMNS
from is_boundary_matching import JitterEstimate, call_boundary_events


def make_unified(rows):
    records = []
    for boundary_id, center_bp in rows:
        records.append({
            "chrom": "A02",
            "start": center_bp - 4000,
            "end": center_bp + 4000,
            "unified_boundary_id": boundary_id,
            "center_bp": center_bp,
            "component_state": "simple",
            "signal_support": "both",
            "raw_boundary_ids": f"raw_{boundary_id}",
            "balanced_boundary_ids": f"balanced_{boundary_id}",
            "support_resolutions": "4000,8000",
            "n_supporting_resolutions": 2,
            "n_source_boundaries": 2,
        })
    return pd.DataFrame(records, columns=UNIFIED_BOUNDARY_COLUMNS)


def build_truth_fixture(n_per_class, conserved_offsets, shifts):
    rows_a = []
    rows_b = []
    truth = {}
    for index in range(n_per_class):
        block = 400000 + index * 1600000

        conserved_a = f"a_conserved_{index}"
        conserved_b = f"b_conserved_{index}"
        rows_a.append((conserved_a, block))
        rows_b.append((conserved_b, block + conserved_offsets[index % len(conserved_offsets)]))
        truth[(conserved_a, conserved_b)] = "conserved"

        shifted_a = f"a_shifted_{index}"
        shifted_b = f"b_shifted_{index}"
        rows_a.append((shifted_a, block + 400000))
        rows_b.append((shifted_b, block + 400000 + shifts[index % len(shifts)]))
        truth[(shifted_a, shifted_b)] = "shifted"

        lost_a = f"a_lost_{index}"
        rows_a.append((lost_a, block + 800000))
        truth[(lost_a, None)] = "lost_in_B"

        gained_b = f"b_gained_{index}"
        rows_b.append((gained_b, block + 1200000))
        truth[(None, gained_b)] = "gained_in_B"
    return make_unified(rows_a), make_unified(rows_b), truth


def event_key(row):
    boundary_a = None if pd.isna(row.boundary_a_id) else str(row.boundary_a_id)
    boundary_b = None if pd.isna(row.boundary_b_id) else str(row.boundary_b_id)
    return boundary_a, boundary_b


def test_noiseless_simple_events_are_exact(settings):
    boundaries_a, boundaries_b, truth = build_truth_fixture(
        n_per_class=20,
        conserved_offsets=(0,),
        shifts=(40000, 60000, 80000),
    )
    events = call_boundary_events(
        boundaries_a,
        boundaries_b,
        JitterEstimate(20000, "empirical_q95", 200, 18000.0),
        settings,
        "A_vs_B",
        "A",
        "B",
    )
    predicted = {
        event_key(row): row.position_state
        for row in events.itertuples(index=False)
    }
    assert predicted == truth
```

- [ ] **Step 2: Add the empirical-jitter macro-F1 benchmark**

Sample offsets from a fixed empirical-like vector on the 4 kb grid. Score
`conserved`, `gained_in_B`, `lost_in_B`, and `shifted` with
this dependency-free helper:

```python
def macro_f1(expected, observed, labels):
    scores = []
    for label in labels:
        true_positive = sum(e == label and o == label for e, o in zip(expected, observed))
        false_positive = sum(e != label and o == label for e, o in zip(expected, observed))
        false_negative = sum(e == label and o != label for e, o in zip(expected, observed))
        denominator = 2 * true_positive + false_positive + false_negative
        scores.append(0.0 if denominator == 0 else 2 * true_positive / denominator)
    return float(np.mean(scores))
```

Use the fixed technical-offset vector
`np.tile(np.array([0, 4000, 4000, 8000, 12000, 16000]), 40)`, split it between
the two samples, and pass it through `estimate_jitter()`. Build the truth
fixture with `conserved_offsets=(0, 4000, 8000, 12000, 16000)` and the same
40/60/80 kb shifts. Align observed and expected labels by `event_key()` and
require macro-F1 at least 0.90 for simple events meeting the specification's
spacing and shift-separation conditions. Do not add scikit-learn to the Conda
environment for this one metric.

Add a crowded-neighborhood regression with A centers at 100/108 kb and B
centers at 104/112 kb. Their 8 kb support intervals form a tied many-to-many
component, so assert one `complex` event rather than two arbitrary simple
matches. Add a second crowded but uniquely ordered case at A centers 300/340
kb and B centers 304/344 kb and assert two noncrossing `conserved` events.

Run:

```bash
python -m pytest -q test/is_boundary/test_synthetic_benchmark.py
```

Expected: both truth benchmarks pass.

- [ ] **Step 3: Document interpretation and targeted execution**

Add a README section named `IS-change and TAD-boundary annotation` containing:

- the two-axis position/strength interpretation;
- A-to-B gained/lost direction;
- raw/balanced provenance policy;
- 20 kb fallback jitter and 100 kb primary shift radius;
- primary versus preliminary output separation;
- the statement that cross-method agreement is not an independent replicate;
- all stable output paths;
- the targeted H1398-versus-L134 Snakemake command;
- the optional configuration switch for set-level permutation;
- links to the design, implementation plan, ADR-0001, and ADR-0002.

Document a deterministic visual-audit procedure: filter the figure manifest by
simple `position_state`, sort with seed `20260903`, and inspect up to 20 loci
per sufficiently populated class. Record `boundary_event_id`, reviewer,
`pass/fail`, and note in
`results/tad_boundary_events/{comparison}.visual_review.tsv`. This is a manual
QC artifact, not a Snakemake completion dependency and not a statistical test.

- [ ] **Step 4: Run the complete focused regression suite**

Run:

```bash
python -m pytest -q test/is_boundary test/is_change test/tad
```

Expected: all selected tests pass with zero failures.

- [ ] **Step 5: Run the production-target dry run**

Run:

```bash
CONDA_OFFLINE=true snakemake -n -p --use-conda --conda-prefix "$PWD/.snakemake/conda" --configfile config/config.yaml -- results/multiscale/r4000/H1398_vs_L134/is_change/H1398_vs_L134_is_change_loci.boundary_annotated.tsv results/multiscale/r4000/H1398_vs_L134/is_change/figures/manifest.tsv
```

Expected: exit code 0; the DAG includes required TAD, IS, unification, event,
annotation, and plotting jobs; the permutation rule is absent while disabled.

- [ ] **Step 6: Verify schemas, direction symmetry, and formatting**

Run:

```bash
python -m pytest -q test/is_boundary/test_matching.py -k "direction or schema"
git diff --check
```

Expected: selected tests pass and `git diff --check` prints no errors.

- [ ] **Step 7: Commit benchmark and user documentation**

```bash
git add README.md test/is_boundary/test_synthetic_benchmark.py
git commit -m "test: validate IS boundary annotation end to end"
```

- [ ] **Step 8: Review the complete branch before execution on production data**

Run:

```bash
git status --short
git log --oneline --decorate -10
```

Expected: no uncommitted changes created by these tasks, and the nine feature
commits appear in task order. Existing user-owned changes that predated plan
execution remain untouched and must be reported separately rather than added
to any feature commit.

## Completion Criteria

Implementation is complete only when all of the following are evidenced by
fresh command output:

- every `test/is_boundary` test passes;
- existing `test/is_change` and `test/tad` suites pass;
- the H1398-versus-L134 annotated target and diagnostic manifest dry-run with
  exit code 0;
- the disabled permutation rule is not scheduled by default;
- matching QC reports raw-only, balanced-only, every leave-one-resolution-out
  view, and 60/80/100 kb shift-radius concordance;
- strength QC reports categorical calls at robust-Z thresholds 1.5, 2.0, and
  2.5 without replacing the primary 2.0 call;
- the direction-swap invariant passes for every simple synthetic event;
- noiseless canonical cases are exact and empirical-jitter macro-F1 is at
  least 0.90 under the specification's separation conditions;
- primary and preliminary output schemas are both stable and separate;
- shifted-event association tests prove that the endpoint corridor is not
  treated as a boundary;
- ADR-0001 source-output separation remains intact;
- the deterministic visual-review manifest procedure is documented, with any
  completed manual review kept outside automated pass/fail criteria;
- `git diff --check` reports no formatting errors;
- final commits contain only files named by their corresponding tasks.
