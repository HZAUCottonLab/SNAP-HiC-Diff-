from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
import pandas as pd
import pytest


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from is_change_detection import (
    build_valid_mask,
    classify_scale_events,
    cluster_multiscale_events,
    score_scale,
)


@pytest.fixture
def settings():
    return SimpleNamespace(
        resolution_bp=4000,
        local_radius_bins=1,
        min_local_valid_fraction=2 / 3,
        flank_fraction_of_scale=0.25,
        min_flank_bins=2,
        min_flank_valid_fraction=0.80,
        max_missing_run_bins=1,
        high_percentile=0.95,
        broad_percentile=0.90,
        spatial_cluster_tolerance_bp=12000,
        min_primary_scale_bands=2,
        scale_bands={
            "small": (32000, 60000, 92000),
            "medium": (120000, 152000, 180000),
            "large": (212000, 240000, 272000, 300000),
        },
    )


def score_track_from_z(z, resolution_bp, scale_bp, settings):
    """Minimal scored geometry fixture for classification-only tests."""
    z = np.asarray(z, dtype=float)
    coordinates = np.arange(len(z)) * resolution_bp
    nonzero = np.abs(z) > 0
    amplitude_percentile = np.where(nonzero, 0.95, 0.10)
    transition_percentile = np.full(len(z), 0.10)
    pulse_percentile = np.full(len(z), 0.10)
    pulse_supported = np.zeros(len(z), dtype=bool)

    changes = np.flatnonzero(np.diff(z) != 0) + 1
    if len(changes) == 1 and z[changes[0] - 1] * z[changes[0]] < 0:
        transition_percentile[changes[0]] = 0.99

    starts = np.flatnonzero(np.r_[nonzero[0], nonzero[1:] & ~nonzero[:-1]])
    ends = np.flatnonzero(np.r_[nonzero[:-1] & ~nonzero[1:], nonzero[-1]])
    for start, end in zip(starts, ends):
        if end - start + 1 >= 2 and np.all(np.sign(z[start : end + 1]) == np.sign(z[start])):
            pulse_percentile[start : end + 1] = 0.99
            pulse_supported[start : end + 1] = True

    track = pd.DataFrame(
        {
            "chrom": "chr1",
            "coordinate": coordinates,
            "scale_bp": scale_bp,
            "z": z,
            "amplitude_percentile": amplitude_percentile,
            "transition_percentile": transition_percentile,
            "pulse_percentile": pulse_percentile,
            "pulse_direction_supported": pulse_supported,
        }
    )
    track.attrs["settings"] = settings
    return track


def event(chrom="chr1", scale_bp=32000, summit=100000, event_type="transition"):
    return {
        "chrom": chrom,
        "scale_bp": scale_bp,
        "summit": summit,
        "event_start": summit,
        "event_end": summit,
        "event_type": event_type,
        "amplitude_percentile": 0.96,
        "transition_percentile": 0.99 if event_type == "transition" else 0.10,
        "pulse_percentile": 0.99 if event_type == "focal_pulse" else 0.10,
    }


def test_zero_is_is_valid_but_nan_is_is_missing():
    a = np.array([0.0, np.nan, -1.0, -1.0])
    b = np.array([0.0, 0.2, 1.0, 1.0])

    valid = build_valid_mask(a, b, np.ones(4, bool), np.ones(4, bool))

    assert valid.tolist() == [True, False, True, True]


def test_transition_is_localized_at_state_change(settings):
    z = np.array([-2.0, -2.0, -2.0, -2.0, 2.0, 2.0, 2.0, 2.0])

    events = classify_scale_events(score_track_from_z(z, 4000, 32000, settings))

    assert events.loc[events["event_type"] == "transition", "summit"].tolist() == [16000]


def test_pulse_is_not_a_transition(settings):
    z = np.array([0.0, 0.0, 0.0, 3.0, 3.0, 0.0, 0.0, 0.0])

    events = classify_scale_events(score_track_from_z(z, 4000, 32000, settings))

    assert set(events["event_type"]) == {"focal_pulse"}


def test_singleton_spike_is_not_a_pulse(settings):
    z = np.array([0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0])

    events = classify_scale_events(score_track_from_z(z, 4000, 32000, settings))

    assert events.empty


def test_broad_amplitude_run_is_reported_as_an_interval(settings):
    track = pd.DataFrame(
        {
            "chrom": "chr1",
            "coordinate": np.arange(8) * 4000,
            "scale_bp": 32000,
            "amplitude_percentile": [0.10, 0.10, 0.92, 0.92, 0.92, 0.92, 0.10, 0.10],
            "transition_percentile": 0.10,
            "pulse_percentile": 0.10,
            "pulse_direction_supported": False,
        }
    )
    track.attrs["settings"] = settings

    events = classify_scale_events(track)

    assert events[["event_type", "event_start", "event_end"]].to_dict("records") == [
        {"event_type": "broad_shift", "event_start": 8000, "event_end": 24000}
    ]

def test_discontinuous_broad_runs_are_not_compacted(settings):
    track = pd.DataFrame(
        {
            "chrom": "chr1",
            "coordinate": np.arange(12) * 4000,
            "scale_bp": 32000,
            "amplitude_percentile": [
                0.92, 0.92, 0.92, 0.92, 0.10, 0.10,
                0.92, 0.92, 0.92, 0.92, 0.10, 0.10,
            ],
            "transition_percentile": 0.10,
            "pulse_percentile": 0.10,
            "pulse_direction_supported": False,
        }
    )
    track.attrs["settings"] = settings

    events = classify_scale_events(track)

    assert events[["event_start", "event_end"]].to_dict("records") == [
        {"event_start": 0, "event_end": 16000},
        {"event_start": 24000, "event_end": 40000},
    ]

def test_eighty_percent_valid_flank_remains_evaluable(settings):
    n_bins = 16
    balanced = np.ones(n_bins, dtype=bool)
    balanced[5] = False
    track = score_scale(
        "chr1",
        np.arange(n_bins) * settings.resolution_bp,
        np.zeros(n_bins),
        np.arange(n_bins, dtype=float),
        balanced,
        np.ones(n_bins, dtype=bool),
        80000,
        settings,
    )

    row = track.loc[7]
    assert row["transition_evaluable"]
    assert row["pulse_evaluable"]


def test_two_bin_invalid_flank_run_only_disables_flank_evidence(settings):
    n_bins = 16
    balanced = np.ones(n_bins, dtype=bool)
    balanced[[4, 5]] = False
    track = score_scale(
        "chr1",
        np.arange(n_bins) * settings.resolution_bp,
        np.zeros(n_bins),
        np.arange(n_bins, dtype=float),
        balanced,
        np.ones(n_bins, dtype=bool),
        80000,
        settings,
    )

    row = track.loc[7]
    assert row["amplitude_evaluable"]
    assert not row["transition_evaluable"]
    assert not row["pulse_evaluable"]


def test_scales_from_one_band_are_not_independent_support(settings):
    events = pd.DataFrame(
        [
            event(scale_bp=32000, summit=100000),
            event(scale_bp=60000, summit=104000),
        ]
    )

    clustered = cluster_multiscale_events(events, settings)

    assert clustered.loc[0, "support_scale_bands"] == "small"
    assert clustered.loc[0, "confidence_tier"] == "single_band_preliminary"


def test_two_distinct_bands_make_a_primary_call_within_twelve_kb(settings):
    events = pd.DataFrame(
        [
            event(scale_bp=32000, summit=100000),
            event(scale_bp=120000, summit=112000),
        ]
    )

    clustered = cluster_multiscale_events(events, settings)

    assert len(clustered) == 1
    assert clustered.loc[0, "support_scale_bands"] == "small,medium"
    assert clustered.loc[0, "confidence_tier"] == "primary"


def test_distant_point_events_remain_separate_loci(settings):
    events = pd.DataFrame(
        [
            event(scale_bp=32000, summit=100000),
            event(scale_bp=120000, summit=116000),
        ]
    )

    clustered = cluster_multiscale_events(events, settings)

    assert clustered["summit"].tolist() == [100000, 116000]


def test_compatible_transition_and_pulse_cluster_as_mixed(settings):
    events = pd.DataFrame(
        [
            event(scale_bp=32000, summit=100000, event_type="transition"),
            event(scale_bp=120000, summit=104000, event_type="focal_pulse"),
        ]
    )

    clustered = cluster_multiscale_events(events, settings)

    assert clustered.loc[0, "event_type"] == "mixed"
