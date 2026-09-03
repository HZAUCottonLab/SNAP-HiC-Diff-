from pathlib import Path
import sys

import pandas as pd
import pytest
import subprocess


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from is_boundary_annotation import (
    annotate_loci,
    associated_event_ids,
    evidence_tier,
    summarize_annotations,
)
from is_boundary_io import (
    ANNOTATION_COLUMNS,
    BOUNDARY_EVENT_COLUMNS,
    STRENGTH_LONG_COLUMNS,
    STRENGTH_SUMMARY_COLUMNS,
    atomic_write_tables,
)
from call_is_change_loci import PRIMARY_COLUMNS


def _event(**updates):
    row = {
        "comparison": "A_vs_B", "sample_a": "A", "sample_b": "B", "chrom": "chr1",
        "boundary_event_id": "e1", "position_state": "conserved",
        "position_substate": "conserved_exact", "boundary_a_id": "a1", "boundary_b_id": "b1",
        "a_start": 100000, "a_end": 104000, "b_start": 100000, "b_end": 104000,
        "event_span_start": 100000, "event_span_end": 104000, "shift_bp": 0,
        "match_stage": "fixed", "match_is_reciprocal": True, "second_best_margin_bp": 8000,
        "n_local_candidates": 1, "jitter_threshold_bp": 4000, "jitter_source": "empirical_q95",
        "max_shift_bp": 100000, "ambiguity_reason": "", "state_at_60000": "conserved",
        "state_at_80000": "conserved", "state_at_100000": "conserved",
        "position_state_stable": True, "signal_support_a": "both", "signal_support_b": "both",
        "resolution_support_a": "4000,8000", "resolution_support_b": "4000,8000",
    }
    row.update(updates)
    return row


def _locus(**updates):
    row = {
        "comparison": "A_vs_B", "region_id": "r1", "chrom": "chr1", "event_start": 100000,
        "event_end": 104000, "summit": 100000, "is_change_type": "transition",
        "confidence_tier": "primary", "support_scale_bands": "small,medium",
    }
    row.update(updates)
    return row


def _strength(**updates):
    row = {
        "comparison": "A_vs_B", "boundary_event_id": "e1", "position_state": "conserved",
        "strength_state": "weakened_in_B", "strength_support_bands": "small,medium",
        "is_consistency_state": "not_applicable", "small_band_score": -3.0,
        "medium_band_score": -3.0, "large_band_score": 0.0, "n_evaluable_bands": 3,
        "state_at_z1_5": "weakened_in_B", "state_at_z2_0": "weakened_in_B",
        "state_at_z2_5": "weakened_in_B", "strength_state_stable": True,
    }
    row.update(updates)
    return row


def _masked(rows=()):
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "reason"])


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


def test_annotation_keeps_one_row_per_locus_and_handles_masked_absence():
    loci = pd.DataFrame([
        _locus(region_id="none", summit=300000, event_start=300000, event_end=304000),
        _locus(region_id="masked", summit=300000, event_start=300000, event_end=304000),
        _locus(region_id="unique"),
        _locus(region_id="multiple"),
        _locus(region_id="complex"),
        _locus(region_id="gain_masked"),
    ])
    events = pd.DataFrame([
        _event(),
        _event(boundary_event_id="e2", a_start=100000, a_end=104000, b_start=100000, b_end=104000),
        _event(boundary_event_id="complex", position_state="complex", a_start=200000, a_end=204000, b_start=200000, b_end=204000),
        _event(boundary_event_id="gain", position_state="gained_in_B", a_start=float("nan"), a_end=float("nan"), b_start=400000, b_end=404000),
    ], columns=BOUNDARY_EVENT_COLUMNS)
    strength = pd.DataFrame([_strength(), _strength(boundary_event_id="e2")])

    clean = annotate_loci(loci.iloc[[0, 2, 3, 4]], events, strength, _masked(), type("S", (), {"base_bin_bp": 4000})())
    masked = annotate_loci(loci.iloc[[1]], events, strength, _masked([("chr1", 300000, 304000, "IS_NA_A")]), type("S", (), {"base_bin_bp": 4000})())
    gain_masked = annotate_loci(
        loci.iloc[[5]].assign(summit=400000, event_start=400000, event_end=404000), events, strength,
        _masked([("chr1", 400000, 404000, "IS_NA_A")]), type("S", (), {"base_bin_bp": 4000})(),
    )

    assert clean.columns.tolist() == list(loci.columns) + list(ANNOTATION_COLUMNS)
    assert len(clean) == 4
    assert clean.set_index("region_id").loc["none", "position_state"] == "no_boundary_association"
    assert clean.set_index("region_id").loc["none", "boundary_event_id"] == ""
    assert masked.iloc[0].position_state == "ambiguous"
    assert masked.iloc[0].annotation_reason == "mask_overlap:IS_NA_A"
    assert clean.set_index("region_id").loc["multiple", "position_state"] == "complex"
    assert clean.set_index("region_id").loc["multiple", "n_boundary_candidates"] == 2
    assert clean.set_index("region_id").loc["complex", "position_state"] == "complex"
    assert gain_masked.iloc[0].position_state == "ambiguous"


def test_tad_coverage_mask_blocks_zero_candidate_no_association():
    loci = pd.DataFrame([_locus(summit=300000, event_start=300000, event_end=304000)])
    tad_unevaluable = pd.DataFrame([
        {"sample": "B", "chrom": "chr1", "start": 300000, "end": 304000, "reason": "nonfinite_base_bins"},
    ])

    result = annotate_loci(
        loci, pd.DataFrame(columns=BOUNDARY_EVENT_COLUMNS), pd.DataFrame(), _masked(),
        type("S", (), {"base_bin_bp": 4000})(), tad_unevaluable=tad_unevaluable,
    ).iloc[0]

    assert result.position_state == "ambiguous"
    assert result.annotation_reason == "tad_coverage:B:nonfinite_base_bins"


def test_conserved_position_can_carry_a_weakened_strength_state():
    result = annotate_loci(
        pd.DataFrame([_locus()]), pd.DataFrame([_event()], columns=BOUNDARY_EVENT_COLUMNS),
        pd.DataFrame([_strength()]), _masked(), type("S", (), {"base_bin_bp": 4000})(),
    ).iloc[0]

    assert result.position_state == "conserved"
    assert result.strength_state == "weakened_in_B"
    assert result.evidence_tier == "high"


@pytest.mark.parametrize(
    ("event_updates", "locus_updates", "mask_reason"),
    [
        (
            {"position_state": "gained_in_B", "a_start": float("nan"), "a_end": float("nan"), "b_start": 400000, "b_end": 404000},
            {"summit": 400000, "event_start": 400000, "event_end": 404000}, "IS_NA_A",
        ),
        (
            {"position_state": "lost_in_B", "a_start": 400000, "a_end": 404000, "b_start": float("nan"), "b_end": float("nan")},
            {"summit": 400000, "event_start": 400000, "event_end": 404000}, "IS_NA_B",
        ),
    ],
)
def test_masked_gain_or_loss_is_ambiguous(event_updates, locus_updates, mask_reason):
    result = annotate_loci(
        pd.DataFrame([_locus(**locus_updates)]),
        pd.DataFrame([_event(**event_updates)], columns=BOUNDARY_EVENT_COLUMNS),
        pd.DataFrame([_strength(
            position_state=event_updates["position_state"], strength_state="not_applicable",
            is_consistency_state="unevaluable", strength_state_stable=False,
        )]),
        _masked([("chr1", 400000, 404000, mask_reason)]),
        type("S", (), {"base_bin_bp": 4000})(),
    ).iloc[0]

    assert result.position_state == "ambiguous"
    assert result.annotation_reason == f"mask_overlap:{mask_reason}"
    assert result.evidence_tier == "not_applicable"


@pytest.mark.parametrize(
    ("event_updates", "locus_updates"),
    [
        (
            {"position_state": "gained_in_B", "a_start": float("nan"), "a_end": float("nan"), "b_start": 400000, "b_end": 404000},
            {"summit": 400000, "event_start": 400000, "event_end": 404000},
        ),
        (
            {"position_state": "lost_in_B", "a_start": 400000, "a_end": 404000, "b_start": float("nan"), "b_end": float("nan")},
            {"summit": 400000, "event_start": 400000, "event_end": 404000},
        ),
    ],
)
def test_clean_unevaluable_gain_or_loss_keeps_position_with_low_evidence(event_updates, locus_updates):
    result = annotate_loci(
        pd.DataFrame([_locus(**locus_updates)]),
        pd.DataFrame([_event(**event_updates)], columns=BOUNDARY_EVENT_COLUMNS),
        pd.DataFrame([_strength(
            position_state=event_updates["position_state"], strength_state="not_applicable",
            strength_support_bands="", is_consistency_state="unevaluable",
            strength_state_stable=False,
        )]),
        _masked(), type("S", (), {"base_bin_bp": 4000})(),
    ).iloc[0]

    assert result.position_state == event_updates["position_state"]
    assert result.is_consistency_state == "unevaluable"
    assert result.evidence_tier == "low"


def test_annotation_canonicalizes_only_the_current_legacy_is_comparison(tmp_path):
    primary = pd.DataFrame([_locus(comparison="A_IS__vs__B_IS")], columns=PRIMARY_COLUMNS)
    preliminary = pd.DataFrame(columns=PRIMARY_COLUMNS)
    result = subprocess.run(
        _annotation_command(
            tmp_path, pd.DataFrame([_event()], columns=BOUNDARY_EVENT_COLUMNS), primary, preliminary
        ), text=True, capture_output=True, check=False,
    )

    assert result.returncode == 0, result.stderr
    annotated = pd.read_csv(tmp_path / "annotated.primary.tsv", sep="\t")
    assert annotated.loc[0, "comparison"] == "A_vs_B"


@pytest.mark.parametrize(
    ("updates", "expected"),
    [
        ({}, "high"),
        ({"is_consistency_state": "unevaluable"}, "moderate"),
        ({"signal_support_a": "raw_only", "resolution_support_a": "4000"}, "low"),
        ({"position_state_stable": False}, "moderate"),
        ({"jitter_source": "fallback_max_tad_resolution"}, "low"),
        ({"confidence_tier": "single_band_preliminary"}, "low"),
        ({"position_state": "complex"}, "not_applicable"),
        ({"position_state": "ambiguous"}, "not_applicable"),
    ],
)
def test_evidence_tier_outcomes(updates, expected):
    row = {
        "position_state": "conserved", "confidence_tier": "primary",
        "jitter_source": "empirical_q95", "signal_support_a": "both", "signal_support_b": "both",
        "resolution_support_a": "4000,8000", "resolution_support_b": "4000,8000",
        "is_consistency_state": "not_applicable", "strength_support_bands": "small,medium",
        "strength_state_stable": True, "position_state_stable": True,
    }
    row.update(updates)
    assert evidence_tier(pd.Series(row)) == expected


@pytest.mark.parametrize(
    ("position_state", "consistency", "support"),
    [
        ("gained_in_B", "direction_concordant", "small,medium"),
        ("lost_in_B", "direction_concordant", "small,medium"),
        ("shifted", "bidirectional_shift_support", "old_a:small,medium;new_b:small,medium"),
    ],
)
def test_direction_supported_rearrangements_can_reach_high_evidence(
    position_state, consistency, support
):
    row = pd.Series({
        "position_state": position_state, "confidence_tier": "primary",
        "jitter_source": "empirical_q95", "signal_support_a": "both", "signal_support_b": "both",
        "resolution_support_a": "4000,8000", "resolution_support_b": "4000,8000",
        "is_consistency_state": consistency, "strength_support_bands": support,
        "strength_state_stable": True, "position_state_stable": True,
    })

    assert evidence_tier(row) == "high"


def test_summary_uses_confidence_tier_denominators_and_explicit_state_order():
    primary = pd.DataFrame([
        {"confidence_tier": "primary", "position_state": "shifted", "strength_state": "not_applicable", "is_consistency_state": "direction_concordant", "evidence_tier": "high"},
        {"confidence_tier": "primary", "position_state": "conserved", "strength_state": "weakened_in_B", "is_consistency_state": "not_applicable", "evidence_tier": "high"},
    ])
    preliminary = pd.DataFrame([
        {"confidence_tier": "single_band_preliminary", "position_state": "conserved", "strength_state": "stable_strength", "is_consistency_state": "not_applicable", "evidence_tier": "low"},
    ])

    summary = summarize_annotations(primary, preliminary)

    assert summary.columns.tolist() == [
        "confidence_tier", "position_state", "strength_state", "is_consistency_state",
        "evidence_tier", "n_loci", "denominator", "fraction",
    ]
    assert summary.loc[summary.confidence_tier.eq("primary"), "denominator"].eq(2).all()
    assert summary.iloc[0].position_state == "conserved"


def test_annotation_qc_counts_each_mask_reason():
    from is_boundary_annotation import annotation_qc

    loci = pd.DataFrame([
        {"annotation_reason": "mask_overlap:IS_NA_A", "n_boundary_candidates": 0},
        {"annotation_reason": "mask_overlap:IS_NA_A,FLANK_COVERAGE", "n_boundary_candidates": 1},
        {"annotation_reason": "unique_boundary_event", "n_boundary_candidates": 1},
    ])
    events = pd.DataFrame({"jitter_source": ["empirical_q95"], "position_state_stable": [True]})
    strength_long = pd.DataFrame({"robust_z": [float("nan"), 1.0]})
    strength_summary = pd.DataFrame({"strength_state_stable": [True]})

    qc = annotation_qc(loci, loci.iloc[0:0], events, strength_long, strength_summary)

    reasons = qc[qc.metric.eq("masked_loci_reason")].set_index("detail")["value"].to_dict()
    assert reasons == {"FLANK_COVERAGE": 1, "IS_NA_A": 2}


def test_atomic_table_write_removes_temporary_files_on_serialization_failure(tmp_path, monkeypatch):
    original = pd.DataFrame.to_csv
    calls = 0

    def fail_second(self, *args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("simulated serialization failure")
        return original(self, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_csv", fail_second)
    with pytest.raises(OSError, match="simulated"):
        atomic_write_tables([
            (pd.DataFrame({"value": [1]}), tmp_path / "one.tsv", ("value",)),
            (pd.DataFrame({"value": [2]}), tmp_path / "two.tsv", ("value",)),
        ])
    assert not list(tmp_path.glob(".*.tmp"))


def _write_tracks(tmp_path):
    rows = [
        {"chrom": "chr1", "viewpoint": start, "IS_distance": scale, "IS": 0.0}
        for scale in (32000, 60000, 92000, 120000, 152000, 180000, 212000, 240000, 272000, 300000)
        for start in range(0, 408000, 4000)
    ]
    for name in ("a", "b"):
        pd.DataFrame(rows).to_csv(tmp_path / f"{name}.is.tsv", sep="\t", index=False)


def _annotation_command(tmp_path, events, primary, preliminary):
    events.to_csv(tmp_path / "events.tsv", sep="\t", index=False)
    primary.to_csv(tmp_path / "primary.tsv", sep="\t", index=False)
    preliminary.to_csv(tmp_path / "preliminary.tsv", sep="\t", index=False)
    (tmp_path / "masked.tsv").write_text("")
    _write_tracks(tmp_path)
    return [
        sys.executable, str(REPO / "scripts" / "annotate_is_change_boundaries.py"),
        "--comparison", "A_vs_B", "--events", str(tmp_path / "events.tsv"),
        "--primary", str(tmp_path / "primary.tsv"), "--preliminary", str(tmp_path / "preliminary.tsv"),
        "--masked", str(tmp_path / "masked.tsv"), "--is-a", str(tmp_path / "a.is.tsv"),
        "--is-b", str(tmp_path / "b.is.tsv"), "--config", str(REPO / "config" / "config.yaml"),
        "--out-primary", str(tmp_path / "annotated.primary.tsv"),
        "--out-preliminary", str(tmp_path / "annotated.preliminary.tsv"),
        "--out-strength", str(tmp_path / "strength.tsv"),
        "--out-strength-summary", str(tmp_path / "strength-summary.tsv"),
        "--out-qc", str(tmp_path / "qc.tsv"), "--out-summary", str(tmp_path / "summary.tsv"),
        "--log", str(tmp_path / "annotation.log"),
    ]


def test_annotation_cli_writes_valid_header_only_outputs(tmp_path):
    empty_events = pd.DataFrame(columns=BOUNDARY_EVENT_COLUMNS)
    empty_loci = pd.DataFrame(columns=PRIMARY_COLUMNS)
    result = subprocess.run(
        _annotation_command(tmp_path, empty_events, empty_loci, empty_loci),
        text=True, capture_output=True, check=False,
    )

    assert result.returncode == 0, result.stderr
    assert pd.read_csv(tmp_path / "annotated.primary.tsv", sep="\t").columns.tolist() == (
        list(PRIMARY_COLUMNS) + list(ANNOTATION_COLUMNS)
    )
    assert pd.read_csv(tmp_path / "strength.tsv", sep="\t").columns.tolist() == list(STRENGTH_LONG_COLUMNS)
    assert pd.read_csv(tmp_path / "strength-summary.tsv", sep="\t").columns.tolist() == list(STRENGTH_SUMMARY_COLUMNS)
    assert not list(tmp_path.glob(".*.tmp"))


def test_annotation_cli_preserves_input_columns_before_annotation_columns(tmp_path):
    primary = pd.DataFrame([_locus()], columns=PRIMARY_COLUMNS)
    preliminary = pd.DataFrame(columns=PRIMARY_COLUMNS)
    result = subprocess.run(
        _annotation_command(
            tmp_path, pd.DataFrame([_event()], columns=BOUNDARY_EVENT_COLUMNS), primary, preliminary
        ), text=True, capture_output=True, check=False,
    )

    assert result.returncode == 0, result.stderr
    observed = pd.read_csv(tmp_path / "annotated.primary.tsv", sep="\t")
    assert observed.columns.tolist() == list(primary.columns) + list(ANNOTATION_COLUMNS)
    assert observed.loc[0, "region_id"] == "r1"
    assert observed.loc[0, "boundary_event_id"] == "e1"
