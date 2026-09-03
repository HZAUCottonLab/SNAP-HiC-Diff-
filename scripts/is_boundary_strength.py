"""Measure multiscale IS-valley strength at TAD-boundary event coordinates."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_prominences

from is_boundary_io import STRENGTH_LONG_COLUMNS, STRENGTH_SUMMARY_COLUMNS


ROLE_BY_STATE = {
    "conserved": ("shared",),
    "gained_in_B": ("new_b",),
    "lost_in_B": ("old_a",),
    "shifted": ("old_a", "new_b"),
    "complex": ("diagnostic",),
    "ambiguous": ("diagnostic",),
}
_BANDS = ("small", "medium", "large")


@dataclass(frozen=True)
class StrengthObservation:
    strength: float
    valley_start: float
    detection: str
    valid_fraction: float


def _finite_runs(mask):
    indices = np.flatnonzero(mask)
    if not len(indices):
        return ()
    return tuple(part for part in np.split(indices, np.flatnonzero(np.diff(indices) > 1) + 1) if len(part))


def _valid_fraction(values):
    return float(np.isfinite(values).mean()) if len(values) else 0.0


def _max_invalid_run(valid):
    run = maximum = 0
    for observed in valid:
        if observed:
            run = 0
        else:
            run += 1
            maximum = max(maximum, run)
    return maximum


def _evaluable(values, fraction, max_missing_run_bins):
    valid = np.isfinite(values)
    return (
        len(valid) > 0
        and valid.mean() >= fraction
        and _max_invalid_run(valid) <= max_missing_run_bins
    )


def local_valley_strength(
    track, start, end, scale_bp, min_valid_fraction, min_flank_valid_fraction=None,
    resolution_bp=4000, max_missing_run_bins=1,
):
    """Return the strongest plateau-aware IS valley inside one boundary interval."""
    radius = int(scale_bp)
    flank_threshold = min_valid_fraction if min_flank_valid_fraction is None else min_flank_valid_fraction
    if track.empty:
        return _missing_observation()
    lower = max(0, int(start) - radius)
    first = ((lower + resolution_bp - 1) // resolution_bp) * resolution_bp
    upper = min(int(end) + radius, int(track.start.max()) + resolution_bp)
    starts = np.arange(first, upper, resolution_bp, dtype=int)
    values = (
        track.set_index("start").IS.reindex(starts).to_numpy(dtype=float)
        if len(starts) else np.array([], dtype=float)
    )
    finite = np.isfinite(values)
    valid_fraction = _valid_fraction(values)
    central = values[(starts >= start) & (starts < end)]
    left = values[(starts >= start - radius) & (starts < start)]
    right = values[(starts >= end) & (starts < end + radius)]
    if (
        not _evaluable(values, min_valid_fraction, max_missing_run_bins)
        or not _evaluable(central, min_valid_fraction, max_missing_run_bins)
        or not _evaluable(left, flank_threshold, max_missing_run_bins)
        or not _evaluable(right, flank_threshold, max_missing_run_bins)
    ):
        return StrengthObservation(float("nan"), float("nan"), "insufficient_coverage", valid_fraction)

    candidates = []
    for run in _finite_runs(finite):
        run_values = values[run]
        peaks, _ = find_peaks(-run_values, plateau_size=True)
        for peak, prominence in zip(peaks, peak_prominences(-run_values, peaks)[0]):
            index = run[peak]
            if start <= starts[index] < end:
                candidates.append((float(prominence), int(starts[index])))
    if not candidates:
        return StrengthObservation(0.0, float("nan"), "no_valley", valid_fraction)
    strength, valley_start = max(candidates, key=lambda item: (item[0], -item[1]))
    return StrengthObservation(strength, valley_start, "valley", valid_fraction)


def _missing_observation():
    return StrengthObservation(float("nan"), float("nan"), "insufficient_coverage", 0.0)


def _intervals(event, role):
    a = (event.a_start, event.a_end) if pd.notna(event.a_start) and pd.notna(event.a_end) else None
    b = (event.b_start, event.b_end) if pd.notna(event.b_start) and pd.notna(event.b_end) else None
    if role == "new_b":
        return b, b
    if role == "old_a":
        return a, a
    return a, b


def _track_index(track):
    if isinstance(track, dict):
        return track
    return {
        (chrom, int(scale)): group[["start", "IS"]].sort_values("start", kind="stable")
        for (chrom, scale), group in track.groupby(["chrom", "IS_distance"], sort=False, observed=True)
    }


def _observation(track, interval, radius, is_settings):
    if interval is None:
        return _missing_observation()
    return local_valley_strength(
        track, int(interval[0]), int(interval[1]), radius,
        is_settings.min_local_valid_fraction, is_settings.min_flank_valid_fraction,
        is_settings.resolution_bp, is_settings.max_missing_run_bins,
    )


def _band_by_scale(is_settings):
    return {
        int(scale): band
        for band, scales in is_settings.scale_bands.items()
        for scale in scales
    }


def _standardize_strength_long(strength_long, boundary_settings):
    if strength_long.empty:
        return strength_long.loc[:, STRENGTH_LONG_COLUMNS].copy()
    result = strength_long.copy()
    result[["robust_z", "reference_n", "reference_median", "reference_mad"]] = np.nan
    reference = result.loc[
        (result.position_state == "conserved")
        & (result.coordinate_role == "shared")
        & np.isfinite(result.delta_strength),
        ["IS_distance", "delta_strength"],
    ]
    for scale, indices in result.groupby("IS_distance", sort=False).groups.items():
        values = reference.loc[reference.IS_distance == scale, "delta_strength"].to_numpy(dtype=float)
        count = len(values)
        median = float(np.median(values)) if count else float("nan")
        mad = float(np.median(np.abs(values - median))) if count else float("nan")
        result.loc[indices, ["reference_n", "reference_median", "reference_mad"]] = count, median, mad
        if count < boundary_settings.strength.min_reference_boundaries or not np.isfinite(mad) or mad == 0:
            continue
        finite = np.isfinite(result.loc[indices, "delta_strength"])
        target = np.asarray(indices)[finite]
        result.loc[target, "robust_z"] = (
            result.loc[target, "delta_strength"] - median
        ) / (1.4826 * mad)
    return result.loc[:, STRENGTH_LONG_COLUMNS]


def compute_strength_long(events, track_a, track_b, boundary_settings, is_settings):
    """Evaluate both balanced IS tracks for every event role and configured scale."""
    tracks_a, tracks_b = _track_index(track_a), _track_index(track_b)
    band_by_scale = _band_by_scale(is_settings)
    rows = []
    for event in events.itertuples(index=False):
        for role in ROLE_BY_STATE[event.position_state]:
            interval_a, interval_b = _intervals(event, role)
            for scale in is_settings.scales_bp:
                radius = int(round(int(scale) * boundary_settings.strength.local_radius_scale_factor))
                observed_a = _observation(
                    tracks_a.get((event.chrom, int(scale)), pd.DataFrame(columns=["start", "IS"])),
                    interval_a, radius, is_settings,
                )
                observed_b = _observation(
                    tracks_b.get((event.chrom, int(scale)), pd.DataFrame(columns=["start", "IS"])),
                    interval_b, radius, is_settings,
                )
                rows.append({
                    "comparison": event.comparison,
                    "boundary_event_id": event.boundary_event_id,
                    "chrom": event.chrom,
                    "position_state": event.position_state,
                    "coordinate_role": role,
                    "IS_distance": int(scale),
                    "scale_band": band_by_scale[int(scale)],
                    "a_start": interval_a[0] if interval_a else np.nan,
                    "a_end": interval_a[1] if interval_a else np.nan,
                    "b_start": interval_b[0] if interval_b else np.nan,
                    "b_end": interval_b[1] if interval_b else np.nan,
                    "strength_a": observed_a.strength,
                    "strength_b": observed_b.strength,
                    "delta_strength": (
                        observed_b.strength - observed_a.strength
                        if np.isfinite(observed_a.strength) and np.isfinite(observed_b.strength)
                        else np.nan
                    ),
                    "robust_z": np.nan,
                    "valley_a_start": observed_a.valley_start,
                    "valley_b_start": observed_b.valley_start,
                    "detection_a": observed_a.detection,
                    "detection_b": observed_b.detection,
                    "valid_fraction_a": observed_a.valid_fraction,
                    "valid_fraction_b": observed_b.valid_fraction,
                    "reference_n": np.nan,
                    "reference_median": np.nan,
                    "reference_mad": np.nan,
                })
    return _standardize_strength_long(
        pd.DataFrame(rows, columns=STRENGTH_LONG_COLUMNS), boundary_settings
    )


def _band_scores(rows, role, boundary_settings):
    scores = {}
    for band in _BANDS:
        values = rows.loc[
            (rows.coordinate_role == role) & (rows.scale_band == band),
            ["IS_distance", "robust_z"],
        ].groupby("IS_distance", sort=False).robust_z.median().dropna()
        scores[band] = (
            float(values.median())
            if len(values) >= boundary_settings.strength.min_valid_scales_per_band
            else float("nan")
        )
    return scores


def _strength_state(scores, threshold, stable_threshold):
    valid = {band: value for band, value in scores.items() if np.isfinite(value)}
    positive = [band for band, value in valid.items() if value >= threshold]
    negative = [band for band, value in valid.items() if value <= -threshold]
    if len(positive) >= 2:
        return "strengthened_in_B", positive
    if len(negative) >= 2:
        return "weakened_in_B", negative
    if len(valid) >= 2 and all(abs(value) < stable_threshold for value in valid.values()):
        return "stable_strength", []
    return "discordant_or_unevaluable", []


def _direction_state(scores, sign, threshold):
    valid = [value for value in scores.values() if np.isfinite(value)]
    if len(valid) < 2:
        return "unevaluable"
    expected = sum(value * sign >= threshold for value in valid) >= 2
    opposite = sum(value * sign <= -threshold for value in valid) >= 2
    if expected:
        return "supported"
    if opposite:
        return "opposite"
    return "absent"


def _is_consistency(position_state, scores_by_role, threshold):
    if position_state in {"conserved", "complex", "ambiguous"}:
        return "not_applicable"
    if position_state == "gained_in_B":
        return "direction_concordant" if _direction_state(scores_by_role["new_b"], 1, threshold) == "supported" else (
            "unevaluable" if _direction_state(scores_by_role["new_b"], 1, threshold) == "unevaluable" else "discordant"
        )
    if position_state == "lost_in_B":
        state = _direction_state(scores_by_role["old_a"], -1, threshold)
        return "direction_concordant" if state == "supported" else "unevaluable" if state == "unevaluable" else "discordant"
    old = _direction_state(scores_by_role["old_a"], -1, threshold)
    new = _direction_state(scores_by_role["new_b"], 1, threshold)
    if "unevaluable" in {old, new}:
        return "unevaluable"
    if "opposite" in {old, new}:
        return "discordant"
    if old == new == "supported":
        return "bidirectional_shift_support"
    if "supported" in {old, new}:
        return "partial_shift_support"
    return "discordant"


def summarize_strength(strength_long, events, boundary_settings, is_settings):
    """Classify strength at conserved boundaries and directionality at rearrangements."""
    strength_long = _standardize_strength_long(strength_long, boundary_settings)
    rows = []
    for event in events.itertuples(index=False):
        event_rows = strength_long[strength_long.boundary_event_id == event.boundary_event_id]
        scores_by_role = {
            role: _band_scores(event_rows, role, boundary_settings)
            for role in ROLE_BY_STATE[event.position_state]
        }
        scores = scores_by_role.get("shared") or scores_by_role.get("new_b") or scores_by_role.get("old_a") or {band: np.nan for band in _BANDS}
        if event.position_state == "conserved":
            states = {
                threshold: _strength_state(
                    scores, threshold, boundary_settings.strength.stable_abs_z_threshold
                )[0]
                for threshold in boundary_settings.strength.sensitivity_thresholds
            }
            primary, support = _strength_state(
                scores,
                boundary_settings.strength.robust_z_threshold,
                boundary_settings.strength.stable_abs_z_threshold,
            )
            consistency = "not_applicable"
        else:
            primary, support, consistency = "not_applicable", [], _is_consistency(
                event.position_state, scores_by_role,
                boundary_settings.strength.robust_z_threshold,
            )
            states = {threshold: "not_applicable" for threshold in boundary_settings.strength.sensitivity_thresholds}
        state_names = {f"state_at_z{threshold:.1f}": state for threshold, state in states.items()}
        rows.append({
            "comparison": event.comparison,
            "boundary_event_id": event.boundary_event_id,
            "position_state": event.position_state,
            "strength_state": primary,
            "strength_support_bands": ",".join(support),
            "is_consistency_state": consistency,
            "small_band_score": scores["small"],
            "medium_band_score": scores["medium"],
            "large_band_score": scores["large"],
            "n_evaluable_bands": sum(np.isfinite(value) for value in scores.values()),
            "state_at_z1_5": state_names.get("state_at_z1.5", "not_applicable"),
            "state_at_z2_0": state_names.get("state_at_z2.0", "not_applicable"),
            "state_at_z2_5": state_names.get("state_at_z2.5", "not_applicable"),
            "strength_state_stable": len(set(states.values())) == 1,
        })
    return pd.DataFrame(rows, columns=STRENGTH_SUMMARY_COLUMNS)
