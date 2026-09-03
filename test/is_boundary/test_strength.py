from pathlib import Path
import sys

import numpy as np
import pandas as pd
import yaml


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from cworld_is_config import parse_cworld_is_settings
from is_boundary_config import parse_is_boundary_settings
from is_change_config import parse_is_change_settings
from is_boundary_io import STRENGTH_LONG_COLUMNS, STRENGTH_SUMMARY_COLUMNS
from is_boundary_strength import (
    ROLE_BY_STATE,
    compute_strength_long,
    local_valley_strength,
    summarize_strength,
)
from tad_config import parse_tad_settings


def settings():
    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    cworld = parse_cworld_is_settings(config)
    return (
        parse_is_boundary_settings(
            config, cworld, parse_tad_settings(config)
        ),
        parse_is_change_settings(config, cworld),
    )


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


def test_local_valley_strength_requires_each_flank():
    track = pd.DataFrame({
        "start": np.arange(9) * 4000,
        "IS": [np.nan, np.nan, np.nan, -2.0, 0.4, 0.7, 1.0, 0.9, 0.8],
    })

    observation = local_valley_strength(
        track, 12000, 20000, 32000, 0.50, min_flank_valid_fraction=0.80
    )

    assert np.isnan(observation.strength)
    assert observation.detection == "insufficient_coverage"


def event_frame(rows):
    return pd.DataFrame(rows, columns=[
        "comparison", "boundary_event_id", "chrom", "position_state",
        "a_start", "a_end", "b_start", "b_end",
    ])


def strength_rows(events, is_settings, effects):
    rows = []
    band_by_scale = {
        scale: band
        for band, scales in is_settings.scale_bands.items()
        for scale in scales
    }
    for event in events.itertuples(index=False):
        roles = ROLE_BY_STATE[event.position_state]
        for role in roles:
            for scale in is_settings.scales_bp:
                delta = effects(event.boundary_event_id, role, band_by_scale[scale], scale)
                rows.append({
                    "comparison": event.comparison,
                    "boundary_event_id": event.boundary_event_id,
                    "chrom": event.chrom,
                    "position_state": event.position_state,
                    "coordinate_role": role,
                    "IS_distance": scale,
                    "scale_band": band_by_scale[scale],
                    "a_start": 100000,
                    "a_end": 104000,
                    "b_start": 100000,
                    "b_end": 104000,
                    "strength_a": 2.0,
                    "strength_b": 2.0 + delta,
                    "delta_strength": delta,
                    "robust_z": np.nan,
                    "valley_a_start": 100000,
                    "valley_b_start": 100000,
                    "detection_a": "valley",
                    "detection_b": "valley",
                    "valid_fraction_a": 1.0,
                    "valid_fraction_b": 1.0,
                    "reference_n": np.nan,
                    "reference_median": np.nan,
                    "reference_mad": np.nan,
                })
    return pd.DataFrame(rows, columns=STRENGTH_LONG_COLUMNS)


def reference_events(extra):
    rows = [
        ("A_vs_B", f"reference_{index}", "chr1", "conserved", 0, 4000, 0, 4000)
        for index in range(60)
    ]
    rows.extend(extra)
    return event_frame(rows)


def test_summary_uses_robust_z_two_band_calls_and_threshold_sensitivity():
    boundary_settings, is_settings = settings()
    events = reference_events([
        ("A_vs_B", "event_target", "chr1", "conserved", 0, 4000, 0, 4000),
        ("A_vs_B", "event_weakened", "chr1", "conserved", 0, 4000, 0, 4000),
        ("A_vs_B", "event_stable", "chr1", "conserved", 0, 4000, 0, 4000),
    ])

    def effects(event_id, _role, band, _scale):
        if event_id.startswith("reference_"):
            return (-1.0, 0.0, 1.0)[int(event_id.rsplit("_", 1)[1]) % 3]
        if event_id == "event_target" and band in {"small", "medium"}:
            return 5.0
        if event_id == "event_weakened" and band in {"small", "medium"}:
            return -5.0
        return 0.0

    summary = summarize_strength(
        strength_rows(events, is_settings, effects), events, boundary_settings, is_settings
    )
    target = summary.loc[summary.boundary_event_id == "event_target"].iloc[0]
    weakened = summary.loc[summary.boundary_event_id == "event_weakened"].iloc[0]
    stable = summary.loc[summary.boundary_event_id == "event_stable"].iloc[0]

    assert summary.columns.tolist() == list(STRENGTH_SUMMARY_COLUMNS)
    assert target.strength_state == "strengthened_in_B"
    assert target.strength_support_bands == "small,medium"
    assert weakened.strength_state == "weakened_in_B"
    assert weakened.strength_support_bands == "small,medium"
    assert stable.strength_state == "stable_strength"
    assert target.state_at_z1_5 == "strengthened_in_B"
    assert target.state_at_z2_0 == target.strength_state
    assert target.state_at_z2_5 == "strengthened_in_B"
    assert bool(target.strength_state_stable)


def test_zero_mad_makes_affected_scales_unevaluable():
    boundary_settings, is_settings = settings()
    events = reference_events([
        ("A_vs_B", "event_target", "chr1", "conserved", 0, 4000, 0, 4000),
    ])

    def effects(event_id, _role, band, _scale):
        if band == "small":
            return 5.0 if event_id == "event_target" else 0.0
        if event_id == "event_target" and band == "medium":
            return 5.0
        return (-1.0, 0.0, 1.0)[int(event_id.rsplit("_", 1)[1]) % 3] if event_id.startswith("reference_") else 0.0

    summary = summarize_strength(
        strength_rows(events, is_settings, effects), events, boundary_settings, is_settings
    )
    target = summary.loc[summary.boundary_event_id == "event_target"].iloc[0]

    assert np.isnan(target.small_band_score)
    assert target.strength_state == "discordant_or_unevaluable"


def test_rearrangements_use_directional_is_consistency_states():
    boundary_settings, is_settings = settings()
    events = reference_events([
        ("A_vs_B", "gain", "chr1", "gained_in_B", np.nan, np.nan, 0, 4000),
        ("A_vs_B", "loss", "chr1", "lost_in_B", 0, 4000, np.nan, np.nan),
        ("A_vs_B", "shift", "chr1", "shifted", 0, 4000, 40000, 44000),
        ("A_vs_B", "shift_opposite", "chr1", "shifted", 0, 4000, 40000, 44000),
    ])

    def effects(event_id, role, _band, _scale):
        if event_id.startswith("reference_"):
            return (-1.0, 0.0, 1.0)[int(event_id.rsplit("_", 1)[1]) % 3]
        if event_id == "gain":
            return 5.0
        if event_id == "loss":
            return -5.0
        if event_id == "shift":
            return -5.0 if role == "old_a" else 5.0
        if event_id == "shift_opposite":
            return 5.0
        return 0.0

    summary = summarize_strength(
        strength_rows(events, is_settings, effects), events, boundary_settings, is_settings
    ).set_index("boundary_event_id")

    assert summary.loc["gain", "strength_state"] == "not_applicable"
    assert summary.loc["gain", "is_consistency_state"] == "direction_concordant"
    assert summary.loc["loss", "is_consistency_state"] == "direction_concordant"
    assert summary.loc["shift", "is_consistency_state"] == "bidirectional_shift_support"
    assert summary.loc["shift_opposite", "is_consistency_state"] == "discordant"


def test_compute_strength_long_assigns_each_coordinate_role():
    boundary_settings, is_settings = settings()
    events = event_frame([
        ("A_vs_B", "conserved", "chr1", "conserved", 400000, 404000, 404000, 408000),
        ("A_vs_B", "gained", "chr1", "gained_in_B", np.nan, np.nan, 600000, 604000),
        ("A_vs_B", "lost", "chr1", "lost_in_B", 800000, 804000, np.nan, np.nan),
        ("A_vs_B", "shifted", "chr1", "shifted", 1000000, 1004000, 1200000, 1204000),
        ("A_vs_B", "complex", "chr1", "complex", 1400000, 1404000, 1404000, 1408000),
    ])
    starts = np.arange(0, 1800000, 4000)

    def track(valleys):
        rows = []
        for scale in is_settings.scales_bp:
            values = np.zeros(len(starts))
            values[np.isin(starts, valleys)] = -2.0
            rows.append(pd.DataFrame({
                "chrom": "chr1", "start": starts, "IS_distance": scale, "IS": values,
            }))
        return pd.concat(rows, ignore_index=True)

    track_a = track([400000, 600000, 800000, 1000000, 1200000, 1400000, 1404000])
    track_b = track([404000, 600000, 800000, 1000000, 1200000, 1400000, 1404000])
    strength_long = compute_strength_long(
        events, track_a, track_b, boundary_settings, is_settings
    )

    assert strength_long.columns.tolist() == list(STRENGTH_LONG_COLUMNS)
    roles = strength_long.groupby("boundary_event_id").coordinate_role.agg(list).to_dict()
    assert {event_id: list(dict.fromkeys(value)) for event_id, value in roles.items()} == {
        "complex": ["diagnostic"],
        "conserved": ["shared"],
        "gained": ["new_b"],
        "lost": ["old_a"],
        "shifted": ["old_a", "new_b"],
    }
    assert len(strength_long) == 6 * len(is_settings.scales_bp)
    gained = strength_long.loc[strength_long.boundary_event_id == "gained"].iloc[0]
    lost = strength_long.loc[strength_long.boundary_event_id == "lost"].iloc[0]
    assert (gained.a_start, gained.b_start) == (600000, 600000)
    assert (lost.a_start, lost.b_start) == (800000, 800000)
