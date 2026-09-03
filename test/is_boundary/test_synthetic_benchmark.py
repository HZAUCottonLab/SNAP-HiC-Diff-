"""Deterministic truth benchmarks for IS-change boundary annotation."""

from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import pytest
import yaml


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from cworld_is_config import parse_cworld_is_settings
from is_boundary_annotation import annotate_loci, associated_event_ids
from is_boundary_config import parse_is_boundary_settings
from is_boundary_io import UNIFIED_BOUNDARY_COLUMNS
from is_boundary_matching import (
    JitterEstimate,
    align_ordered,
    call_boundary_events,
    estimate_jitter,
)
from tad_config import parse_tad_settings


@pytest.fixture
def settings():
    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    return parse_is_boundary_settings(
        config, parse_cworld_is_settings(config), parse_tad_settings(config)
    )


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
        rows_b.append((
            conserved_b,
            block + conserved_offsets[index % len(conserved_offsets)],
        ))
        truth[(conserved_a, conserved_b)] = "conserved"

        shifted_a = f"a_shifted_{index}"
        shifted_b = f"b_shifted_{index}"
        rows_a.append((shifted_a, block + 400000))
        rows_b.append((
            shifted_b,
            block + 400000 + shifts[index % len(shifts)],
        ))
        truth[(shifted_a, shifted_b)] = "shifted"

        lost_a = f"a_lost_{index}"
        rows_a.append((lost_a, block + 800000))
        truth[(lost_a, None)] = "lost_in_B"

        gained_b = f"b_gained_{index}"
        rows_b.append((gained_b, block + 1200000))
        truth[(None, gained_b)] = "gained_in_B"
    return make_unified(rows_a), make_unified(rows_b), truth


def event_key(row):
    return _identifier(row.boundary_a_id), _identifier(row.boundary_b_id)


def _identifier(value):
    return None if pd.isna(value) or value == "" else str(value)


def macro_f1(expected, observed, labels):
    scores = []
    for label in labels:
        true_positive = sum(e == label and o == label for e, o in zip(expected, observed))
        false_positive = sum(e != label and o == label for e, o in zip(expected, observed))
        false_negative = sum(e == label and o != label for e, o in zip(expected, observed))
        denominator = 2 * true_positive + false_positive + false_negative
        scores.append(0.0 if denominator == 0 else 2 * true_positive / denominator)
    return float(np.mean(scores))


def aligned_truth_labels(truth, predicted):
    assert set(predicted) == set(truth)
    return [truth[key] for key in truth], [predicted[key] for key in truth]


def _call(boundaries_a, boundaries_b, jitter, settings, comparison="A_vs_B"):
    return call_boundary_events(
        boundaries_a, boundaries_b, jitter, settings, comparison, "A", "B"
    )


def _predicted_states(events):
    return {event_key(row): row.position_state for row in events.itertuples(index=False)}


def test_noiseless_simple_events_are_exact_and_directionally_symmetric(settings):
    boundaries_a, boundaries_b, truth = build_truth_fixture(
        n_per_class=20,
        conserved_offsets=(0,),
        shifts=(40000, 60000, 80000),
    )
    events = _call(
        boundaries_a,
        boundaries_b,
        JitterEstimate(20000, "empirical_q95", 200, 18000.0),
        settings,
    )
    predicted = _predicted_states(events)

    assert predicted == truth
    assert set(predicted.values()) == {
        "conserved", "shifted", "lost_in_B", "gained_in_B"
    }

    reverse = call_boundary_events(
        boundaries_b,
        boundaries_a,
        JitterEstimate(20000, "empirical_q95", 200, 18000.0),
        settings,
        "B_vs_A",
        "B",
        "A",
    )
    reverse_by_key = {event_key(row): row for row in reverse.itertuples(index=False)}
    reverse_state = {
        "conserved": "conserved",
        "shifted": "shifted",
        "lost_in_B": "gained_in_B",
        "gained_in_B": "lost_in_B",
    }
    for event in events.itertuples(index=False):
        reversed_event = reverse_by_key[(
            _identifier(event.boundary_b_id),
            _identifier(event.boundary_a_id),
        )]
        assert reversed_event.position_state == reverse_state[event.position_state]
        if event.position_state == "shifted":
            assert event.shift_bp == -reversed_event.shift_bp


def test_empirical_jitter_truth_benchmark_has_high_macro_f1(settings):
    technical_offsets = np.tile(np.array([0, 4000, 4000, 8000, 12000, 16000]), 40)
    offsets_a = pd.DataFrame({"offset_bp": technical_offsets[:120], "eligible": True})
    offsets_b = pd.DataFrame({"offset_bp": technical_offsets[120:], "eligible": True})
    jitter = estimate_jitter(offsets_a, offsets_b, settings)
    boundaries_a, boundaries_b, truth = build_truth_fixture(
        n_per_class=40,
        conserved_offsets=(0, 4000, 8000, 12000, 16000),
        shifts=(40000, 60000, 80000),
    )
    predicted = _predicted_states(_call(boundaries_a, boundaries_b, jitter, settings))
    expected, observed = aligned_truth_labels(truth, predicted)

    assert jitter.source == "empirical_q95"
    assert jitter.threshold_bp == 16000
    assert macro_f1(
        expected,
        observed,
        ("conserved", "gained_in_B", "lost_in_B", "shifted"),
    ) >= 0.90


def test_truth_alignment_rejects_extra_predicted_event():
    truth = {("a", "b"): "conserved"}
    predicted = {**truth, (None, "extra"): "gained_in_B"}

    with pytest.raises(AssertionError):
        aligned_truth_labels(truth, predicted)


def test_crowded_components_are_complex_only_when_tied(settings):
    tied = _call(
        make_unified([("a1", 100000), ("a2", 108000)]),
        make_unified([("b1", 104000), ("b2", 112000)]),
        JitterEstimate(4000, "empirical_q95", 200, 4000.0),
        settings,
    )
    ordered = _call(
        make_unified([("a1", 300000), ("a2", 340000)]),
        make_unified([("b1", 304000), ("b2", 344000)]),
        JitterEstimate(4000, "empirical_q95", 200, 4000.0),
        settings,
    )

    assert tied.position_state.tolist() == ["complex"]
    assert tied.loc[0, "n_local_candidates"] == 4
    assert ordered.position_state.tolist() == ["conserved", "conserved"]


def test_shifted_endpoint_only_annotation_respects_mask_and_is_tier(settings):
    events = _call(
        make_unified([("a", 100000)]),
        make_unified([("b", 180000)]),
        JitterEstimate(20000, "empirical_q95", 200, 18000.0),
        settings,
    )
    shifted = events.iloc[0]
    between = pd.Series({
        "is_change_type": "transition", "summit": 140000,
        "event_start": 140000, "event_end": 144000,
    })
    endpoint = between.copy()
    endpoint["summit"] = 100000
    endpoint["event_start"] = 100000
    endpoint["event_end"] = 104000
    loci = pd.DataFrame([
        {
            "comparison": "A_vs_B", "region_id": "corridor", "chrom": "A02",
            **between.to_dict(), "confidence_tier": "primary",
            "support_scale_bands": "small,medium",
        },
        {
            "comparison": "A_vs_B", "region_id": "endpoint", "chrom": "A02",
            **endpoint.to_dict(), "confidence_tier": "single_band_preliminary",
            "support_scale_bands": "small",
        },
    ])
    masked = pd.DataFrame([
        {"chrom": "A02", "start": 140000, "end": 144000, "reason": "IS_NA_A"}
    ])
    annotated = annotate_loci(
        loci,
        events,
        pd.DataFrame(),
        masked,
        type("Settings", (), {"base_bin_bp": 4000})(),
    ).set_index("region_id")

    assert shifted.position_state == "shifted"
    assert associated_event_ids(between, events, 4000) == []
    assert associated_event_ids(endpoint, events, 4000) == [shifted.boundary_event_id]
    assert annotated.loc["corridor", "position_state"] == "ambiguous"
    assert annotated.loc["corridor", "annotation_reason"] == "mask_overlap:IS_NA_A"
    assert annotated.loc["endpoint", "position_state"] == "shifted"
    assert annotated.loc["endpoint", "evidence_tier"] == "low"


def test_ordered_alignment_scales_with_local_candidates():
    boundaries_a = make_unified([
        (f"a{index}", index * 400000 + 400000) for index in range(1000)
    ])
    boundaries_b = make_unified([
        (f"b{index}", index * 400000 + 472000) for index in range(1000)
    ])

    started = time.perf_counter()
    pairs = align_ordered(boundaries_a, boundaries_b, 100000)

    assert pairs == [(index, index) for index in range(1000)]
    assert time.perf_counter() - started < 2.0
