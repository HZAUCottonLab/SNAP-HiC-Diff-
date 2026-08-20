from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
import pandas as pd
import pytest


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

import is_change_detection
from is_change_detection import (
    build_valid_mask,
    classify_scale_events,
    cluster_multiscale_events,
    refine_point_clusters,
    recover_residual_events,
    residual_extrema,
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

    changes = np.flatnonzero(np.diff(z) != 0) + 1
    if len(changes) == 1 and z[changes[0] - 1] * z[changes[0]] < 0:
        transition_percentile[changes[0]] = 0.99

    starts = np.flatnonzero(np.r_[nonzero[0], nonzero[1:] & ~nonzero[:-1]])
    ends = np.flatnonzero(np.r_[nonzero[:-1] & ~nonzero[1:], nonzero[-1]])
    for start, end in zip(starts, ends):
        if end - start + 1 >= 2 and np.all(np.sign(z[start : end + 1]) == np.sign(z[start])):
            pulse_percentile[start : end + 1] = 0.99

    track = pd.DataFrame(
        {
            "chrom": "chr1",
            "coordinate": coordinates,
            "scale_bp": scale_bp,
            "z": z,
            "amplitude_percentile": amplitude_percentile,
            "transition_percentile": transition_percentile,
            "pulse_percentile": pulse_percentile,
        }
    )
    track.attrs["settings"] = settings
    return track


def aligned_score_tracks(settings, z_by_scale, valid_by_scale=None, amplitude=None):
    """Build aligned score_scale-like frames for configured IS scales."""
    valid_by_scale = valid_by_scale or {}
    amplitude = amplitude or {}
    n_bins = len(next(iter(z_by_scale.values())))
    coordinates = np.arange(n_bins) * settings.resolution_bp
    return {
        scale_bp: pd.DataFrame(
            {
                "chrom": "chr1",
                "coordinate": coordinates,
                "scale_bp": scale_bp,
                "valid": valid_by_scale.get(scale_bp, np.ones(n_bins, dtype=bool)),
                "z": z_by_scale[scale_bp],
                "amplitude_percentile": amplitude.get(
                    scale_bp, np.full(n_bins, settings.high_percentile)
                ),
            }
        )
        for scale_bp in settings.scales_bp
    }


def event(chrom="chr1", scale_bp=32000, summit=100000, event_type="transition"):
    return {
        "chrom": chrom,
        "scale_bp": scale_bp,
        "summit": summit,
        "event_start": summit,
        "event_end": summit + 4000,
        "event_type": event_type,
        "amplitude_percentile": 0.96,
        "transition_percentile": 0.99 if event_type == "transition" else 0.10,
        "pulse_percentile": 0.99 if event_type == "focal_pulse" else 0.10,
    }



def point_cluster(summit, support_bands="small,medium", support_scales="32000,120000"):
    return {
        "chrom": "chr1",
        "event_type": "transition",
        "event_start": summit,
        "event_end": summit + 4000,
        "summit": summit,
        "confidence_tier": "primary",
        "support_scale_bands": support_bands,
        "support_scales": support_scales,
        "n_supporting_scales": len(support_scales.split(",")),
        "max_amplitude_percentile": 0.96,
        "max_transition_percentile": 0.99,
        "max_pulse_percentile": 0.95,
    }


def consensus_and_extrema(score_tracks, settings):
    consensus = is_change_detection.build_consensus_track(score_tracks, settings)
    return consensus, is_change_detection.local_consensus_extrema(consensus, settings)


def test_zero_is_is_valid_but_nan_is_is_missing():
    a = np.array([0.0, np.nan, -1.0, -1.0])
    b = np.array([0.0, 0.2, 1.0, 1.0])

    valid = build_valid_mask(a, b, np.ones(4, bool), np.ones(4, bool))

    assert valid.tolist() == [True, False, True, True]


def test_consensus_median_uses_signed_z_and_common_valid_segments(settings):
    # median([-3, -1, 1]) == -1; a false/NaN row creates a segment break.
    settings.scales_bp = (32000, 120000, 212000)
    tracks = aligned_score_tracks(
        settings,
        {
            32000: [-3.0, -3.0, -3.0],
            120000: [-1.0, np.nan, -1.0],
            212000: [1.0, 1.0, 1.0],
        },
        {120000: [True, False, True]},
    )

    consensus = is_change_detection.build_consensus_track(tracks, settings)

    assert consensus.columns.tolist() == [
        "chrom",
        "coordinate",
        "consensus_z",
        "valid",
        "segment_id",
    ]
    assert consensus["valid"].tolist() == [True, False, True]
    assert np.allclose(consensus["consensus_z"], [-1.0, np.nan, -1.0], equal_nan=True)
    assert consensus["segment_id"].iloc[[0, 2]].tolist() == [0, 1]
    assert pd.isna(consensus["segment_id"].iloc[1])


def test_consensus_rejects_adjacent_chromosome_transition(settings):
    settings.scales_bp = (32000, 120000, 212000)
    tracks = aligned_score_tracks(
        settings,
        {
            32000: [-1.0, 3.0, -1.0],
            120000: [-1.0, 3.0, -1.0],
            212000: [-1.0, 3.0, -1.0],
        },
    )
    for track in tracks.values():
        track.loc[1:, "chrom"] = "chr2"

    with pytest.raises(ValueError, match="one chromosome"):
        is_change_detection.build_consensus_track(tracks, settings)


def test_plateau_aware_extrema_emit_one_peak_and_valley_with_prominence(settings):
    # A valley, two-bin peak plateau, and valley expose limiting-flank prominence.
    consensus = pd.DataFrame(
        {
            "chrom": "chr1",
            "coordinate": np.arange(7) * settings.resolution_bp,
            "consensus_z": [-1.0, 3.0, -2.0, 3.0, 3.0, -1.0, 3.0],
            "valid": True,
            "segment_id": 0,
        }
    )

    extrema = is_change_detection.local_consensus_extrema(consensus, settings)

    peak = extrema.loc[extrema["plateau_start"] == 12000].iloc[0]
    valley = extrema.loc[extrema["coordinate"] == 8000].iloc[0]
    assert peak[["direction", "plateau_start", "plateau_end"]].to_dict() == {
        "direction": "peak",
        "plateau_start": 12000,
        "plateau_end": 16000,
    }
    assert peak[["left_flank_coordinate", "right_flank_coordinate"]].to_dict() == {
        "left_flank_coordinate": 8000,
        "right_flank_coordinate": 20000,
    }
    assert peak["prominence"] == 4.0
    assert valley["direction"] == "valley"
    assert valley["prominence"] == 5.0
    assert valley[["left_flank_coordinate", "right_flank_coordinate"]].to_dict() == {
        "left_flank_coordinate": 4000,
        "right_flank_coordinate": 12000,
    }


def test_segment_edges_and_invalid_gaps_have_no_recovery_prominence(settings):
    # Neither condition can supply both opposite-polarity flanks.
    consensus = pd.DataFrame(
        {
            "chrom": "chr1",
            "coordinate": np.arange(7) * settings.resolution_bp,
            "consensus_z": [-1.0, 3.0, -1.0, np.nan, -1.0, 3.0, -1.0],
            "valid": [True, True, True, False, True, True, True],
            "segment_id": [0, 0, 0, np.nan, 1, 1, 1],
        }
    )

    extrema = is_change_detection.local_consensus_extrema(consensus, settings)

    assert 12000 not in extrema["coordinate"].tolist()
    assert extrema["coordinate"].tolist() == [4000, 20000]
    assert extrema["prominence"].isna().all()


def test_directional_support_requires_distinct_bands_and_high_amplitude_when_requested(settings):
    settings.scales_bp = (32000, 60000, 120000)
    tracks = aligned_score_tracks(
        settings,
        {
            32000: [1.0, 1.0],
            60000: [1.0, 1.0],
            120000: [-1.0, -1.0],
        },
    )

    assert is_change_detection.directional_support(tracks, 0, "peak", settings) == (
        ["small"],
        [32000, 60000],
    )
    tracks[120000].loc[:, "z"] = 1.0
    assert is_change_detection.directional_support(tracks, 0, "peak", settings) == (
        ["small", "medium"],
        [32000, 60000, 120000],
    )
    for track in tracks.values():
        track.loc[:, "amplitude_percentile"] = 0.10
    assert is_change_detection.directional_support(
        tracks, 0, "peak", settings, require_high_amplitude=True
    ) == ([], [])
    assert is_change_detection.directional_support(tracks, 4000, "peak", settings) == (
        ["small", "medium"],
        [32000, 60000, 120000],
    )
    assert is_change_detection.directional_support(tracks, 9999, "peak", settings) == (
        [],
        [],
    )


def test_transition_is_localized_at_state_change(settings):
    z = np.array([-2.0, -2.0, -2.0, -2.0, 2.0, 2.0, 2.0, 2.0])

    events = classify_scale_events(score_track_from_z(z, 4000, 32000, settings))

    assert events.loc[events["event_type"] == "transition", "summit"].tolist() == [16000]


def test_pulse_is_not_a_transition(settings):
    z = np.array([0.0, 0.0, 0.0, 3.0, 3.0, 0.0, 0.0, 0.0])

    events = classify_scale_events(score_track_from_z(z, 4000, 32000, settings))

    assert set(events["event_type"]) == {"focal_pulse"}


def test_score_scale_calls_same_direction_peak_and_trough_pulses(settings):
    settings.local_radius_bins = 0
    settings.high_percentile = 0.80
    settings.broad_percentile = 1.10
    delta = np.array([-1., 1., -1., -1., 1., 10., 10., -1., 1., -1., -10., -10., 1., -1., 1.])
    track = score_scale(
        "chr1",
        np.arange(len(delta)) * settings.resolution_bp,
        np.zeros(len(delta)),
        delta,
        np.ones(len(delta), dtype=bool),
        np.ones(len(delta), dtype=bool),
        32000,
        settings,
    )

    events = classify_scale_events(track)
    pulses = events.loc[events["event_type"] == "focal_pulse"]
    assert {20000, 40000}.issubset(set(pulses["summit"]))
    assert (pulses["event_end"] - pulses["event_start"] == settings.resolution_bp).all()


def test_singleton_spike_is_not_a_pulse(settings):
    settings.local_radius_bins = 0
    settings.high_percentile = 0.80
    settings.broad_percentile = 1.10
    delta = np.array([-1., 1., -1., -1., 1., 10., -1., 1., -1., 1., -1., 1., -1., 1., -1.])
    track = score_scale(
        "chr1",
        np.arange(len(delta)) * settings.resolution_bp,
        np.zeros(len(delta)),
        delta,
        np.ones(len(delta), dtype=bool),
        np.ones(len(delta), dtype=bool),
        32000,
        settings,
    )

    events = classify_scale_events(track)

    assert events.empty
    assert not track.loc[5, "pulse_evaluable"]


def test_two_valid_bins_in_a_three_bin_local_window_are_evaluable(settings):
    balanced = np.ones(16, dtype=bool)
    balanced[6] = False
    track = score_scale(
        "chr1",
        np.arange(16) * settings.resolution_bp,
        np.zeros(16),
        np.arange(16, dtype=float),
        balanced,
        np.ones(16, dtype=bool),
        80000,
        settings,
    )

    assert track.loc[7, "amplitude_evaluable"]

def test_broad_amplitude_run_is_reported_as_an_interval(settings):
    track = pd.DataFrame(
        {
            "chrom": "chr1",
            "coordinate": np.arange(8) * 4000,
            "scale_bp": 32000,
            "amplitude_percentile": [0.10, 0.10, 0.92, 0.92, 0.92, 0.92, 0.10, 0.10],
            "transition_percentile": 0.10,
            "pulse_percentile": 0.10,
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
        }
    )
    track.attrs["settings"] = settings

    events = classify_scale_events(track)

    assert events[["event_start", "event_end"]].to_dict("records") == [
        {"event_start": 0, "event_end": 16000},
        {"event_start": 24000, "event_end": 40000},
    ]

def test_long_broad_run_produces_one_interval(settings):
    n_bins = 1000
    track = pd.DataFrame(
        {
            "chrom": "chr1",
            "coordinate": np.arange(n_bins) * settings.resolution_bp,
            "scale_bp": 32000,
            "amplitude_percentile": np.full(n_bins, 0.92),
            "transition_percentile": 0.10,
            "pulse_percentile": 0.10,
        }
    )
    track.attrs["settings"] = settings

    events = classify_scale_events(track)

    assert events[["event_type", "event_start", "event_end"]].to_dict("records") == [
        {"event_type": "broad_shift", "event_start": 0, "event_end": n_bins * settings.resolution_bp}
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
    assert clustered.loc[0, "event_start"] == 100000
    assert clustered.loc[0, "event_end"] == 116000
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


def test_adjacent_broad_intervals_remain_separate(settings):
    events = pd.DataFrame(
        [
            {"chrom": "chr1", "scale_bp": 32000, "event_type": "broad_shift", "summit": np.nan, "event_start": 0, "event_end": 16000},
            {"chrom": "chr1", "scale_bp": 120000, "event_type": "broad_shift", "summit": np.nan, "event_start": 16000, "event_end": 32000},
        ]
    )

    clustered = cluster_multiscale_events(events, settings)

    assert clustered[["event_start", "event_end"]].to_dict("records") == [
        {"event_start": 0, "event_end": 16000},
        {"event_start": 16000, "event_end": 32000},
    ]


def test_compatible_transition_and_pulse_cluster_as_mixed(settings):
    events = pd.DataFrame(
        [
            event(scale_bp=32000, summit=100000, event_type="transition"),
            event(scale_bp=120000, summit=104000, event_type="focal_pulse"),
        ]
    )

    clustered = cluster_multiscale_events(events, settings)

    assert clustered.loc[0, "event_type"] == "mixed"



def test_refinement_moves_a_positive_event_to_nearest_supported_peak(settings):
    settings.scales_bp = (32000, 120000, 212000)
    tracks = aligned_score_tracks(
        settings,
        {scale: [-1.0, 3.0, -1.0, 1.0, 2.0, -1.0, 4.0, -1.0] for scale in settings.scales_bp},
    )
    consensus, extrema = consensus_and_extrema(tracks, settings)

    refined = refine_point_clusters(
        pd.DataFrame([point_cluster(12000)]), consensus, extrema, tracks, settings
    )

    row = refined.iloc[0]
    assert row[["summit", "event_start", "event_end"]].tolist() == [16000, 16000, 20000]
    assert row[["discovery_mode", "original_summit", "refinement_status", "refinement_direction"]].tolist() == [
        "refined", 12000, "refined", "peak"
    ]
    assert row["refinement_support_scale_bands"] == "small,medium,large"


def test_refinement_moves_a_negative_event_to_nearest_supported_valley(settings):
    settings.scales_bp = (32000, 120000, 212000)
    tracks = aligned_score_tracks(
        settings,
        {scale: [1.0, -3.0, 1.0, -1.0, -2.0, 1.0, -4.0, 1.0] for scale in settings.scales_bp},
    )
    consensus, extrema = consensus_and_extrema(tracks, settings)

    refined = refine_point_clusters(
        pd.DataFrame([point_cluster(12000)]), consensus, extrema, tracks, settings
    )

    row = refined.iloc[0]
    assert row[["summit", "event_start", "event_end"]].tolist() == [16000, 16000, 20000]
    assert len(refined) == 1


def test_refinement_chooses_nearest_plateau_member_and_stops_at_invalid_gap(settings):
    settings.scales_bp = (32000, 120000, 212000)
    valid = [True, True, True, True, True, True, False, True, True]
    tracks = aligned_score_tracks(
        settings,
        {scale: [-1.0, -1.0, 1.0, 3.0, 3.0, -1.0, np.nan, 1.0, -1.0] for scale in settings.scales_bp},
        {scale: valid for scale in settings.scales_bp},
    )
    consensus, extrema = consensus_and_extrema(tracks, settings)

    refined = refine_point_clusters(
        pd.DataFrame([point_cluster(8000), point_cluster(28000)]),
        consensus,
        extrema,
        tracks,
        settings,
    ).set_index("original_summit")

    assert refined.loc[8000, "summit"] == 12000
    assert refined.loc[28000, ["summit", "event_start", "event_end"]].tolist() == [28000, 28000, 32000]
    assert refined.loc[28000, "refinement_status"] == "unchanged"


def test_refinement_preserves_ambiguous_or_single_band_event_and_broad_shift(settings):
    settings.scales_bp = (32000, 60000, 120000)
    tracks = aligned_score_tracks(
        settings,
        {
            32000: [0.0, -1.0, 1.0, 3.0, -1.0],
            60000: [0.0, -1.0, 1.0, 3.0, -1.0],
            120000: [0.0, -1.0, -1.0, -1.0, -1.0],
        },
    )
    consensus, extrema = consensus_and_extrema(tracks, settings)
    broad = point_cluster(16000)
    broad.update({"event_type": "broad_shift", "summit": np.nan, "event_end": 24000})

    refined = refine_point_clusters(
        pd.DataFrame([point_cluster(0), point_cluster(8000), broad]),
        consensus,
        extrema,
        tracks,
        settings,
    )
    ambiguous = refined.loc[refined["original_summit"] == 0].iloc[0]
    single_band = refined.loc[refined["original_summit"] == 8000].iloc[0]
    observed_broad = refined.loc[refined["event_type"] == "broad_shift"].iloc[0]

    assert ambiguous["refinement_status"] == "ambiguous"
    assert pd.isna(ambiguous["refinement_direction"])
    assert ambiguous["refinement_support_scale_bands"] == ""
    assert single_band[["summit", "event_start", "event_end"]].tolist() == [8000, 8000, 12000]
    assert single_band[["refinement_status", "refinement_direction", "refinement_support_scale_bands"]].tolist() == [
        "unchanged", "peak", "small"
    ]
    assert observed_broad[["event_start", "event_end", "discovery_mode", "refinement_status"]].tolist() == [
        16000, 24000, "original", "unchanged"
    ]
    assert pd.isna(observed_broad["original_summit"])


def test_refinement_collision_keeps_stronger_multiband_event(settings):
    settings.scales_bp = (32000, 120000, 212000)
    tracks = aligned_score_tracks(
        settings,
        {scale: [-1.0, 1.0, 3.0, 3.0, 1.0, -1.0] for scale in settings.scales_bp},
    )
    consensus, extrema = consensus_and_extrema(tracks, settings)
    weaker = point_cluster(4000, "small", "32000")
    weaker["max_transition_percentile"] = 0.90

    refined = refine_point_clusters(
        pd.DataFrame([weaker, point_cluster(8000)]), consensus, extrema, tracks, settings
    )

    assert len(refined) == 1
    row = refined.iloc[0]
    assert row[["summit", "original_summit", "support_scale_bands"]].tolist() == [
        8000, 8000, "small,medium"
    ]


def recovery_tracks(settings, z_by_scale, amplitude_by_scale=None, evidence=0.95):
    tracks = aligned_score_tracks(
        settings, z_by_scale, amplitude=amplitude_by_scale
    )
    for track in tracks.values():
        track["transition_percentile"] = evidence
        track["pulse_percentile"] = evidence
    return tracks


def test_recovery_adds_supported_peak_and_valley_at_global_prominence_cutoff(settings):
    settings.scales_bp = (32000, 120000, 212000)
    peak_tracks = recovery_tracks(
        settings, {scale: [0.0, -5.0, 5.0, -5.0, 0.0] for scale in settings.scales_bp}
    )
    valley_tracks = recovery_tracks(
        settings, {scale: [0.0, 5.0, -5.0, 5.0, 0.0] for scale in settings.scales_bp}
    )
    peak_consensus, _ = consensus_and_extrema(peak_tracks, settings)
    valley_consensus, _ = consensus_and_extrema(valley_tracks, settings)

    recovered = recover_residual_events(
        {
            "peak": residual_extrema(peak_consensus, set(), settings),
            "valley": residual_extrema(valley_consensus, set(), settings),
        },
        {"peak": peak_tracks, "valley": valley_tracks},
        settings,
    )

    peak = recovered["peak"].iloc[0]
    valley = recovered["valley"].iloc[0]
    assert peak[["event_type", "event_start", "event_end", "summit"]].tolist() == [
        "recovered_peak", 8000, 12000, 8000
    ]
    assert valley[["event_type", "event_start", "event_end", "summit"]].tolist() == [
        "recovered_valley", 8000, 12000, 8000
    ]
    assert peak["confidence_tier"] == valley["confidence_tier"] == "primary"
    assert peak["refinement_support_scale_bands"] == "small,medium,large"
    assert pd.isna(peak["original_summit"])
    assert peak["refinement_status"] == "recovered"


def test_recovery_uses_one_comparison_wide_prominence_cutoff(settings):
    settings.scales_bp = (32000, 120000, 212000)
    tracks = recovery_tracks(
        settings, {scale: [1.0, 1.0, 1.0] for scale in settings.scales_bp}
    )

    def residuals(prominences):
        return pd.DataFrame(
            [
                {
                    "chrom": "chr1",
                    "coordinate": index * settings.resolution_bp,
                    "plateau_start": index * settings.resolution_bp,
                    "plateau_end": index * settings.resolution_bp,
                    "direction": "peak",
                    "segment_id": 0,
                    "prominence": prominence,
                    "left_flank_coordinate": np.nan,
                    "right_flank_coordinate": np.nan,
                }
                for index, prominence in enumerate(prominences)
            ]
        )

    weak, strong = residuals([1.0, 2.0, 3.0]), residuals([10.0, 20.0, 30.0])
    recovered = recover_residual_events(
        {"weak": weak, "strong": strong},
        {"weak": tracks, "strong": tracks},
        settings,
    )

    assert weak["prominence"].quantile(settings.high_percentile) < 3.0
    assert pd.concat((weak, strong))["prominence"].quantile(
        settings.high_percentile
    ) > 3.0
    assert recovered["weak"].empty
    assert recovered["strong"]["summit"].tolist() == [8000]


def test_recovery_rejects_single_band_low_amplitude_and_low_prominence_extrema(settings):
    settings.scales_bp = (32000, 60000, 120000)
    high = [0.0, -5.0, 5.0, -5.0, 0.0]
    low = [0.0, -2.0, 2.0, -2.0, 0.0]
    single_band = recovery_tracks(
        settings,
        {32000: high, 60000: high, 120000: [-0.1 * value for value in high]},
    )
    low_amplitude = recovery_tracks(
        settings,
        {scale: high for scale in settings.scales_bp},
        {scale: np.full(len(high), settings.high_percentile - 0.01) for scale in settings.scales_bp},
    )
    low_prominence = recovery_tracks(
        settings, {scale: low for scale in settings.scales_bp}
    )
    consensuses = {
        name: consensus_and_extrema(tracks, settings)[0]
        for name, tracks in {
            "single_band": single_band,
            "low_amplitude": low_amplitude,
            "low_prominence": low_prominence,
        }.items()
    }

    recovered = recover_residual_events(
        {name: residual_extrema(consensus, set(), settings) for name, consensus in consensuses.items()},
        {
            "single_band": single_band,
            "low_amplitude": low_amplitude,
            "low_prominence": low_prominence,
        },
        settings,
    )

    assert all(events.empty for events in recovered.values())


def test_recovery_excludes_initial_final_occupied_plateaus_and_invalid_bins(settings):
    settings.scales_bp = (32000, 120000, 212000)
    consensus = pd.DataFrame(
        {
            "chrom": "chr1",
            "coordinate": np.arange(12) * settings.resolution_bp,
            "consensus_z": [0.0, 5.0, -5.0, 4.0, 4.0, -5.0, 5.0, 0.0, np.nan, 0.0, 4.0, -4.0],
            "valid": [True] * 8 + [False] + [True] * 3,
            "segment_id": [0.0] * 8 + [np.nan] + [1.0] * 3,
        }
    )

    residuals = residual_extrema(consensus, {12000, 16000}, settings)

    assert not (
        residuals["plateau_start"].le(12000)
        & residuals["plateau_end"].ge(12000)
    ).any()
    assert not (
        residuals["plateau_start"].le(16000)
        & residuals["plateau_end"].ge(16000)
    ).any()
    assert (residuals["plateau_end"] < 32000).all()
    assert np.isfinite(residuals["prominence"]).all()
    assert (residuals["prominence"] >= 0).all()


def test_recovery_never_uses_context_or_segment_edge_extrema(settings):
    settings.scales_bp = (32000, 120000, 212000)
    consensus = pd.DataFrame(
        {
            "chrom": "chr1",
            "coordinate": np.arange(9) * settings.resolution_bp,
            "consensus_z": [5.0, -3.0, 3.0, -3.0, 5.0, np.nan, 5.0, -3.0, 5.0],
            "valid": [True] * 5 + [False] + [True] * 3,
            "segment_id": [0.0] * 5 + [np.nan] + [1.0] * 3,
        }
    )

    residuals = residual_extrema(consensus, set(), settings)

    assert residuals[["plateau_start", "plateau_end"]].to_dict("records") == [
        {"plateau_start": 8000, "plateau_end": 8000}
    ]


def test_recovery_collapses_duplicate_plateau_rows_by_evidence(settings):
    settings.scales_bp = (32000, 120000, 212000)
    z = [0.0, -5.0, 5.0, -5.0, 0.0]
    weaker = recovery_tracks(
        settings, {scale: z for scale in settings.scales_bp}, evidence=0.95
    )
    stronger = recovery_tracks(
        settings, {scale: z for scale in settings.scales_bp}, evidence=0.99
    )
    consensus, _ = consensus_and_extrema(weaker, settings)
    residuals = residual_extrema(consensus, set(), settings)

    recovered = recover_residual_events(
        {"weaker": residuals, "stronger": residuals.copy()},
        {"weaker": weaker, "stronger": stronger},
        settings,
    )

    assert recovered["weaker"].empty
    assert recovered["stronger"]["max_transition_percentile"].tolist() == [0.99]
