from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest
import yaml


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from cworld_is_config import parse_cworld_is_settings
from is_boundary_config import parse_is_boundary_settings
from is_boundary_matching import (
    JitterEstimate,
    _candidate_data,
    align_ordered,
    call_boundary_events,
    compare_event_views,
    estimate_jitter,
)
from is_boundary_io import (
    BOUNDARY_EVENT_COLUMNS,
    BOUNDARY_STABILITY_COLUMNS,
    BOUNDARY_VIEW_COLUMNS,
    UNIFIED_BOUNDARY_COLUMNS,
)
from is_boundary_unification import JITTER_OFFSET_COLUMNS
from tad_config import parse_tad_settings


@pytest.fixture
def settings():
    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    return parse_is_boundary_settings(
        config, parse_cworld_is_settings(config), parse_tad_settings(config)
    )


def unified_frame(rows, signal_support="both", resolutions="4000,8000"):
    records = []
    for chrom, start, end, boundary_id in rows:
        records.append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "unified_boundary_id": boundary_id,
            "center_bp": (start + end) // 2,
            "component_state": "simple",
            "signal_support": signal_support,
            "raw_boundary_ids": f"raw_{boundary_id}",
            "balanced_boundary_ids": f"balanced_{boundary_id}",
            "support_resolutions": resolutions,
            "n_supporting_resolutions": len(resolutions.split(",")),
            "n_source_boundaries": 2,
        })
    return pd.DataFrame(records, columns=UNIFIED_BOUNDARY_COLUMNS)


def view_frame(sample, view, boundaries):
    return boundaries.assign(sample=sample, view=view).loc[:, BOUNDARY_VIEW_COLUMNS]


def offsets_frame(sample, values):
    return pd.DataFrame([
        {
            "sample": sample,
            "chrom": "A02",
            "source": "raw_balanced",
            "comparison": "raw_vs_balanced",
            "resolution_a": "multires",
            "resolution_b": "multires",
            "boundary_a_id": f"a{index}",
            "boundary_b_id": f"b{index}",
            "center_a": 0,
            "center_b": value,
            "offset_bp": value,
            "eligible": True,
        }
        for index, value in enumerate(values)
    ], columns=JITTER_OFFSET_COLUMNS)


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

    assert events.columns.tolist() == list(BOUNDARY_EVENT_COLUMNS)
    assert events["position_state"].tolist() == [
        "conserved", "conserved", "shifted", "lost_in_B", "gained_in_B"
    ]
    assert events["position_substate"].tolist()[:2] == [
        "conserved_exact", "conserved_within_jitter"
    ]
    shifted = events.loc[events.position_state == "shifted"].iloc[0]
    assert shifted.boundary_a_id == "a3"
    assert shifted.boundary_b_id == "b3"
    assert shifted.shift_bp == 72000
    matched_a = events[["a_start", "a_end"]].dropna()
    assert (matched_a["a_start"] <= matched_a["a_end"]).all()
    assert (events["b_start"].dropna() <= events["b_end"].dropna()).all()
    assert align_ordered(a.iloc[[2, 3]], b.iloc[[2, 3]], 100000) == [(0, 0)]


def test_tied_local_component_is_complex_and_large_jitter_preserves_conservation(settings):
    tied_a = unified_frame([("A02", 0, 8000, "a1"), ("A02", 8000, 16000, "a2")])
    tied_b = unified_frame([("A02", 4000, 12000, "b1"), ("A02", 12000, 20000, "b2")])
    tied = call_boundary_events(
        tied_a, tied_b, JitterEstimate(4000, "empirical_q95", 200, 4000.0),
        settings, "H1398_vs_L134", "H1398", "L134",
    )
    assert tied["position_state"].tolist() == ["complex"]
    assert tied.loc[0, "n_local_candidates"] == 4

    a = unified_frame([("A02", 0, 8000, "a")])
    b = unified_frame([("A02", 72000, 80000, "b")])
    conserved = call_boundary_events(
        a, b, JitterEstimate(100000, "empirical_q95", 200, 100000.0),
        settings, "H1398_vs_L134", "H1398", "L134",
    )
    assert conserved.loc[0, "position_state"] == "conserved"
    assert conserved.loc[0, "position_substate"] == "conserved_within_jitter"


def test_event_order_and_ids_do_not_depend_on_input_row_order(settings):
    a = unified_frame([("A02", 0, 8000, "a2"), ("A01", 0, 8000, "a1")])
    b = unified_frame([("A02", 0, 8000, "b2"), ("A01", 0, 8000, "b1")])
    jitter = JitterEstimate(4000, "empirical_q95", 200, 4000.0)

    original = call_boundary_events(a, b, jitter, settings, "H1398_vs_L134", "H1398", "L134")
    reordered = call_boundary_events(
        a.iloc[::-1], b.iloc[::-1], jitter, settings, "H1398_vs_L134", "H1398", "L134"
    )

    assert original[["chrom", "boundary_event_id"]].to_dict("records") == (
        reordered[["chrom", "boundary_event_id"]].to_dict("records")
    )


def test_stage2_uses_residual_candidates_after_fixed_match(settings):
    a = unified_frame([
        ("A02", 36000, 44000, "a_fixed"),
        ("A02", 56000, 64000, "a_shifted"),
    ])
    b = unified_frame([
        ("A02", 36000, 44000, "b_fixed"),
        ("A02", 128000, 136000, "b_shifted"),
    ])

    events = call_boundary_events(
        a, b, JitterEstimate(16000, "empirical_q95", 200, 16000.0), settings,
        "sample_a_vs_sample_b", "sample_a", "sample_b",
    )

    shifted = events.loc[events.boundary_a_id == "a_shifted"].iloc[0]
    assert shifted.position_state == "shifted"
    assert shifted.boundary_b_id == "b_shifted"


def test_ordered_alignment_enumerates_only_local_band_candidates(monkeypatch):
    a = unified_frame([
        ("A02", index * 1_000_000, index * 1_000_000 + 8000, f"a{index}")
        for index in range(100)
    ])
    b = unified_frame([
        ("A02", index * 1_000_000 + 72000, index * 1_000_000 + 80000, f"b{index}")
        for index in range(100)
    ])
    candidates = _candidate_data(a, b, 100000)
    assert sum(len(values) for values in candidates[0].values()) == len(a)

    calls = 0
    original_iterrows = pd.DataFrame.iterrows

    def counted_iterrows(frame, *args, **kwargs):
        nonlocal calls
        calls += 1
        return original_iterrows(frame, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "iterrows", counted_iterrows)
    _candidate_data(a, b, 100000)
    assert calls <= 1
    assert align_ordered(a, b, 100000) == [(index, index) for index in range(100)]


def test_view_stability_output_order_ignores_input_row_order(settings):
    a = unified_frame([("A02", 0, 8000, "a1"), ("A02", 40000, 48000, "a2")])
    b = unified_frame([("A02", 0, 8000, "b1"), ("A02", 40000, 48000, "b2")])
    events = call_boundary_events(
        a, b, JitterEstimate(4000, "empirical_q95", 200, 4000.0), settings,
        "sample_a_vs_sample_b", "sample_a", "sample_b",
    )
    views_a = pd.concat([
        view_frame("sample_a", "raw_only_input", a),
        view_frame("sample_a", "leave_r4000_out", a),
        view_frame("sample_a", "balanced_only_input", a),
    ], ignore_index=True)
    views_b = pd.concat([
        view_frame("sample_b", "raw_only_input", b),
        view_frame("sample_b", "leave_r4000_out", b),
        view_frame("sample_b", "balanced_only_input", b),
    ], ignore_index=True)

    first = compare_event_views(
        events.copy(), views_a, views_b, JitterEstimate(4000, "empirical_q95", 200, 4000.0),
        settings, "sample_a_vs_sample_b", "sample_a", "sample_b",
    )
    reordered = compare_event_views(
        events.iloc[::-1].copy(), views_a.iloc[::-1], views_b.iloc[::-1],
        JitterEstimate(4000, "empirical_q95", 200, 4000.0), settings,
        "sample_a_vs_sample_b", "sample_a", "sample_b",
    )

    assert first.to_dict("records") == reordered.to_dict("records")
    assert first.to_dict("records") == first.sort_values(
        ["comparison", "boundary_event_id", "view"], kind="stable"
    ).to_dict("records")


def test_swapping_direction_reverses_gain_loss_and_shift_sign(settings):
    a = unified_frame([("A02", 0, 8000, "a1"), ("A02", 100000, 108000, "a2")])
    b = unified_frame([("A02", 8000, 16000, "b1"), ("A02", 300000, 308000, "b2")])
    jitter = JitterEstimate(4000, "empirical_q95", 200, 4000.0)
    forward = call_boundary_events(a, b, jitter, settings, "H1398_vs_L134", "H1398", "L134")
    reverse = call_boundary_events(b, a, jitter, settings, "L134_vs_H1398", "L134", "H1398")

    assert forward["position_state"].tolist() == ["shifted", "lost_in_B", "gained_in_B"]
    assert set(reverse["position_state"]) == {"shifted", "lost_in_B", "gained_in_B"}
    assert forward.loc[0, "shift_bp"] == -reverse.loc[0, "shift_bp"]
    assert abs(forward.loc[0, "shift_bp"]) == abs(reverse.loc[0, "shift_bp"])
    assert (
        forward.loc[forward.position_state == "lost_in_B", "boundary_a_id"].item()
        == reverse.loc[reverse.position_state == "gained_in_B", "boundary_b_id"].item()
    )


def test_sensitivity_views_mark_only_agreeing_events_stable(settings):
    a = unified_frame([("A02", 0, 8000, "a1"), ("A02", 40000, 48000, "a2")])
    b = unified_frame([("A02", 0, 8000, "b1"), ("A02", 40000, 48000, "b2")])
    events = call_boundary_events(
        a, b, JitterEstimate(4000, "empirical_q95", 200, 4000.0), settings,
        "H1398_vs_L134", "H1398", "L134",
    )
    views_a = pd.concat([
        view_frame("H1398", "unified", a),
        view_frame("H1398", "raw_only_input", a),
        view_frame("H1398", "leave_r4000_out", a),
    ], ignore_index=True)
    views_b = pd.concat([
        view_frame("L134", "unified", b),
        view_frame("L134", "raw_only_input", b),
        view_frame("L134", "leave_r4000_out", b.iloc[[0]]),
    ], ignore_index=True)

    stability = compare_event_views(
        events, views_a, views_b, JitterEstimate(4000, "empirical_q95", 200, 4000.0),
        settings, "H1398_vs_L134", "H1398", "L134",
    )
    assert stability.columns.tolist() == list(BOUNDARY_STABILITY_COLUMNS)
    assert set(stability["view"]) == {"raw_only_input", "leave_r4000_out"}
    target = stability[stability.boundary_event_id == events.loc[1, "boundary_event_id"]]
    assert target.loc[target.view == "leave_r4000_out", "sensitivity_position_state"].item() == "lost_in_B"
    assert events["position_state_stable"].tolist() == [True, False]


def test_cli_writes_atomic_outputs_and_rejects_same_sample(tmp_path):
    boundaries = unified_frame([("A02", 0, 8000, "a1")])
    boundaries_a = tmp_path / "a.unified.tsv"
    boundaries_b = tmp_path / "b.unified.tsv"
    offsets_a = tmp_path / "a.offsets.tsv"
    offsets_b = tmp_path / "b.offsets.tsv"
    views_a = tmp_path / "a.views.tsv"
    views_b = tmp_path / "b.views.tsv"
    boundaries.to_csv(boundaries_a, sep="\t", index=False)
    boundaries.to_csv(boundaries_b, sep="\t", index=False)
    offsets_frame("H1398", range(0, 400000, 4000)).to_csv(offsets_a, sep="\t", index=False)
    offsets_frame("L134", range(0, 400000, 4000)).to_csv(offsets_b, sep="\t", index=False)
    view_frame("H1398", "unified", boundaries).to_csv(views_a, sep="\t", index=False)
    view_frame("L134", "unified", boundaries).to_csv(views_b, sep="\t", index=False)
    events = tmp_path / "events.tsv"
    stability = tmp_path / "stability.tsv"
    qc = tmp_path / "qc.tsv"
    command = [
        sys.executable, str(REPO / "scripts" / "call_tad_boundary_events.py"),
        "--comparison", "H1398_vs_L134", "--sample-a", "H1398", "--sample-b", "L134",
        "--boundaries-a", str(boundaries_a), "--boundaries-b", str(boundaries_b),
        "--offsets-a", str(offsets_a), "--offsets-b", str(offsets_b),
        "--views-a", str(views_a), "--views-b", str(views_b),
        "--config", str(REPO / "config" / "config.yaml"), "--output-events", str(events),
        "--output-stability", str(stability), "--output-qc", str(qc), "--log", str(tmp_path / "run.log"),
    ]
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    assert result.returncode == 0, result.stderr
    assert pd.read_csv(events, sep="\t").columns.tolist() == list(BOUNDARY_EVENT_COLUMNS)
    assert pd.read_csv(stability, sep="\t").columns.tolist() == list(BOUNDARY_STABILITY_COLUMNS)
    assert pd.read_csv(qc, sep="\t").columns.tolist() == ["metric", "value", "detail"]
    assert not list(tmp_path.glob(".*.tmp"))

    same_sample = command.copy()
    same_sample[same_sample.index("L134")] = "H1398"
    rejected = subprocess.run(same_sample, text=True, capture_output=True, check=False)
    assert rejected.returncode == 1
    assert "IS_BOUNDARY_REFERENCE_MISMATCH" in rejected.stderr

    off_grid = boundaries.copy()
    off_grid.loc[0, ["start", "end", "center_bp"]] = [2000, 10000, 6000]
    off_grid.to_csv(boundaries_b, sep="\t", index=False)
    rejected_grid = subprocess.run(command, text=True, capture_output=True, check=False)
    assert rejected_grid.returncode == 1
    assert "IS_BOUNDARY_REFERENCE_MISMATCH" in rejected_grid.stderr
