"""Describe multiscale insulation-score changes without imputing missing bins."""

from math import ceil

import numpy as np
import pandas as pd


_EVENT_COLUMNS = [
    "chrom",
    "scale_bp",
    "event_type",
    "event_start",
    "event_end",
    "summit",
    "amplitude_percentile",
    "transition_percentile",
    "pulse_percentile",
]
_CLUSTER_COLUMNS = [
    "chrom",
    "event_type",
    "event_start",
    "event_end",
    "summit",
    "confidence_tier",
    "support_scale_bands",
    "support_scales",
    "n_supporting_scales",
    "max_amplitude_percentile",
    "max_transition_percentile",
    "max_pulse_percentile",
]


def build_valid_mask(is_a, is_b, balanced_a, balanced_b):
    """Return bins with finite IS values and valid balanced bins in both samples."""
    return (
        np.isfinite(is_a)
        & np.isfinite(is_b)
        & np.asarray(balanced_a, dtype=bool)
        & np.asarray(balanced_b, dtype=bool)
    )


def robust_zscore(delta, valid):
    """Robustly standardize valid deltas, leaving all unavailable positions missing."""
    out = np.full(np.asarray(delta).shape, np.nan, dtype=float)
    values = np.asarray(delta, dtype=float)[valid]
    if not len(values):
        return out
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    out[valid] = 0.0 if mad == 0 else (values - median) / (1.4826 * mad)
    return out


def flank_bins(scale_bp, settings):
    """Return the fixed genomic flank width for one IS scale."""
    return max(
        settings.min_flank_bins,
        round(scale_bp * settings.flank_fraction_of_scale / settings.resolution_bp),
    )


def _max_invalid_run(valid):
    run = maximum = 0
    for observed in valid:
        if observed:
            run = 0
        else:
            run += 1
            maximum = max(maximum, run)
    return maximum


def _flank_is_evaluable(valid, settings):
    return (
        len(valid) > 0
        and valid.mean() >= settings.min_flank_valid_fraction
        and _max_invalid_run(valid) <= settings.max_missing_run_bins
    )


def _percentiles(values, evaluable):
    out = np.full(len(values), np.nan, dtype=float)
    indices = np.flatnonzero(evaluable & np.isfinite(values))
    if len(indices):
        out[indices] = (
            pd.Series(values[indices]).rank(method="average", pct=True).to_numpy()
        )
    return out


def score_scale(
    chrom,
    coordinates,
    is_a,
    is_b,
    balanced_a,
    balanced_b,
    scale_bp,
    settings,
):
    """Score one fixed-coordinate chromosome/scale IS series.

    IS values are already log2-normalized; this function only differences and
    robustly standardizes them within the supplied chromosome and scale.
    """
    coordinates = np.asarray(coordinates)
    is_a = np.asarray(is_a, dtype=float)
    is_b = np.asarray(is_b, dtype=float)
    balanced_a = np.asarray(balanced_a, dtype=bool)
    balanced_b = np.asarray(balanced_b, dtype=bool)
    n_bins = len(coordinates)
    if any(
        len(values) != n_bins
        for values in (is_a, is_b, balanced_a, balanced_b)
    ):
        raise ValueError("coordinates, IS values, and balanced masks must have equal length")

    valid = build_valid_mask(is_a, is_b, balanced_a, balanced_b)
    delta = is_b - is_a
    z = robust_zscore(delta, valid)
    amplitude = np.full(n_bins, np.nan)
    transition = np.full(n_bins, np.nan)
    pulse = np.full(n_bins, np.nan)
    amplitude_evaluable = np.zeros(n_bins, dtype=bool)
    transition_evaluable = np.zeros(n_bins, dtype=bool)
    pulse_evaluable = np.zeros(n_bins, dtype=bool)
    radius = settings.local_radius_bins
    flank_width = flank_bins(scale_bp, settings)

    for center in np.flatnonzero(valid):
        local_start = max(0, center - radius)
        local_end = min(n_bins, center + radius + 1)
        local_valid = valid[local_start:local_end]
        if local_valid.mean() >= settings.min_local_valid_fraction:
            amplitude[center] = np.median(
                np.abs(z[local_start:local_end][local_valid])
            )
            amplitude_evaluable[center] = True

        left_valid = valid[center - flank_width:center]
        right_valid = valid[center + 1:center + flank_width + 1]
        if (
            len(left_valid) != flank_width
            or len(right_valid) != flank_width
            or not _flank_is_evaluable(left_valid, settings)
            or not _flank_is_evaluable(right_valid, settings)
        ):
            continue

        left_z = z[center - flank_width:center][left_valid]
        right_z = z[center + 1:center + flank_width + 1][right_valid]
        joint_flanks = np.concatenate((left_z, right_z))
        transition[center] = abs(np.median(right_z) - np.median(left_z))
        pulse[center] = z[center] - np.median(joint_flanks)
        transition_evaluable[center] = True
        pulse_evaluable[center] = True

    frame = pd.DataFrame(
        {
            "chrom": chrom,
            "coordinate": coordinates,
            "scale_bp": scale_bp,
            "valid": valid,
            "delta": delta,
            "z": z,
            "amplitude": amplitude,
            "transition": transition,
            "pulse": pulse,
            "amplitude_evaluable": amplitude_evaluable,
            "transition_evaluable": transition_evaluable,
            "pulse_evaluable": pulse_evaluable,
        }
    )
    frame["amplitude_percentile"] = _percentiles(amplitude, amplitude_evaluable)
    frame["transition_percentile"] = _percentiles(
        transition, transition_evaluable
    )
    frame["pulse_percentile"] = _percentiles(pulse, pulse_evaluable)
    frame.attrs["settings"] = settings
    return frame


def label_event(row, settings):
    """Return a descriptive curve-geometry label or None."""
    high_a = row.amplitude_percentile >= settings.high_percentile
    high_c = row.transition_percentile >= settings.high_percentile
    high_p = row.pulse_percentile >= settings.high_percentile
    if high_a and high_c and high_p:
        return "mixed"
    if high_c and row.amplitude_percentile >= settings.broad_percentile:
        return "transition"
    if high_a and high_p:
        return "focal_pulse"
    return None


def _event_row(row, event_type):
    return {
        "chrom": row["chrom"],
        "scale_bp": row["scale_bp"],
        "event_type": event_type,
        "event_start": row["coordinate"],
        "event_end": row["coordinate"],
        "summit": row["coordinate"],
        "amplitude_percentile": row["amplitude_percentile"],
        "transition_percentile": row["transition_percentile"],
        "pulse_percentile": row["pulse_percentile"],
    }


def classify_scale_events(score_track, settings=None):
    """Classify non-opaque point and broad IS-change evidence for one scale."""
    settings = settings or score_track.attrs.get("settings")
    if settings is None:
        raise ValueError("settings are required to classify a score track")
    required = {
        "chrom",
        "coordinate",
        "scale_bp",
        "amplitude_percentile",
        "transition_percentile",
        "pulse_percentile",
    }
    missing = required.difference(score_track.columns)
    if missing:
        raise ValueError(f"score track missing columns: {', '.join(sorted(missing))}")

    track = score_track.reset_index(drop=True)
    point_labels = np.full(len(track), None, dtype=object)
    events = []
    for index, row in track.iterrows():
        event_type = label_event(row, settings)
        point_labels[index] = event_type
        if event_type:
            events.append(_event_row(row, event_type))

    broad_mask = (
        track["amplitude_percentile"].ge(settings.broad_percentile).to_numpy()
        & np.isfinite(track["amplitude_percentile"].to_numpy())
        & (point_labels == None)
    )
    min_broad_bins = max(
        3, ceil(track["scale_bp"].iloc[0] / (2 * settings.resolution_bp))
    )
    broad_indices = np.flatnonzero(broad_mask)
    run_start = 0
    while run_start < len(broad_indices):
        run_end = run_start
        while (
            run_end + 1 < len(broad_indices)
            and broad_indices[run_end + 1] == broad_indices[run_end] + 1
            and track.loc[broad_indices[run_end + 1], "coordinate"]
            - track.loc[broad_indices[run_end], "coordinate"]
            == settings.resolution_bp
        ):
            run_end += 1
        run = broad_indices[run_start:run_end + 1]
        if len(run) < min_broad_bins:
            run_start = run_end + 1
            continue
        rows = track.iloc[run]
        events.append(
            {
                "chrom": rows["chrom"].iloc[0],
                "scale_bp": rows["scale_bp"].iloc[0],
                "event_type": "broad_shift",
                "event_start": rows["coordinate"].iloc[0],
                "event_end": rows["coordinate"].iloc[-1] + settings.resolution_bp,
                "summit": np.nan,
                "amplitude_percentile": rows["amplitude_percentile"].max(),
                "transition_percentile": rows["transition_percentile"].max(),
                "pulse_percentile": rows["pulse_percentile"].max(),
            }
        )
        run_start = run_end + 1
    return pd.DataFrame(events, columns=_EVENT_COLUMNS).sort_values(
        ["event_start", "event_type"], kind="stable"
    ).reset_index(drop=True)


def _scale_band(scale_bp, settings):
    for band, scales in settings.scale_bands.items():
        if scale_bp in scales:
            return band
    raise ValueError(f"scale {scale_bp} is not assigned to a scale band")


def _finite_max(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(values.max()) if len(values) else np.nan


def _support_summary(rows, settings):
    scales = sorted({int(scale) for scale in rows["scale_bp"]})
    bands = [
        band
        for band in settings.scale_bands
        if any(scale in settings.scale_bands[band] for scale in scales)
    ]
    return bands, scales


def _point_summit(rows, event_type):
    candidates = []
    for coordinate, at_coordinate in rows.groupby("summit", sort=False):
        band_count = at_coordinate["_band"].nunique()
        if event_type == "transition":
            values = at_coordinate["transition_percentile"].to_numpy(float)
        elif event_type == "focal_pulse":
            values = at_coordinate["pulse_percentile"].to_numpy(float)
        else:
            values = np.minimum.reduce(
                [
                    at_coordinate[column].to_numpy(float)
                    for column in (
                        "amplitude_percentile",
                        "transition_percentile",
                        "pulse_percentile",
                    )
                ]
            )
        values = values[np.isfinite(values)]
        evidence = float(np.median(values)) if len(values) else -np.inf
        candidates.append((band_count, evidence, coordinate))
    return sorted(candidates, key=lambda value: (-value[0], -value[1], value[2]))[0][2]


def _cluster_point_rows(rows, settings):
    rows = rows.copy()
    rows["_band"] = rows["scale_bp"].map(lambda scale: _scale_band(scale, settings))
    types = set(rows["event_type"])
    event_type = (
        "mixed"
        if "mixed" in types or {"transition", "focal_pulse"}.issubset(types)
        else next(iter(types))
    )
    bands, scales = _support_summary(rows, settings)
    return {
        "chrom": rows["chrom"].iloc[0],
        "event_type": event_type,
        "event_start": rows["event_start"].min(),
        "event_end": rows["event_end"].max(),
        "summit": _point_summit(rows, event_type),
        "confidence_tier": (
            "primary"
            if len(bands) >= settings.min_primary_scale_bands
            else "single_band_preliminary"
        ),
        "support_scale_bands": ",".join(bands),
        "support_scales": ",".join(map(str, scales)),
        "n_supporting_scales": len(scales),
        "max_amplitude_percentile": _finite_max(rows["amplitude_percentile"]),
        "max_transition_percentile": _finite_max(rows["transition_percentile"]),
        "max_pulse_percentile": _finite_max(rows["pulse_percentile"]),
    }


def _cluster_broad_rows(rows, settings):
    rows = rows.copy()
    rows["_band"] = rows["scale_bp"].map(lambda scale: _scale_band(scale, settings))
    bands, scales = _support_summary(rows, settings)
    return {
        "chrom": rows["chrom"].iloc[0],
        "event_type": "broad_shift",
        "event_start": rows["event_start"].min(),
        "event_end": rows["event_end"].max(),
        "summit": np.nan,
        "confidence_tier": (
            "primary"
            if len(bands) >= settings.min_primary_scale_bands
            else "single_band_preliminary"
        ),
        "support_scale_bands": ",".join(bands),
        "support_scales": ",".join(map(str, scales)),
        "n_supporting_scales": len(scales),
        "max_amplitude_percentile": _finite_max(rows["amplitude_percentile"]),
        "max_transition_percentile": _finite_max(rows["transition_percentile"]),
        "max_pulse_percentile": _finite_max(rows["pulse_percentile"]),
    }


def cluster_multiscale_events(events, settings):
    """Cluster point events by tolerance and broad intervals by overlap only."""
    if events.empty:
        return pd.DataFrame(columns=_CLUSTER_COLUMNS)

    events = events.copy()
    if "event_start" not in events:
        events["event_start"] = events["summit"]
    if "event_end" not in events:
        events["event_end"] = events["summit"]
    for column in (
        "amplitude_percentile",
        "transition_percentile",
        "pulse_percentile",
    ):
        if column not in events:
            events[column] = np.nan

    clusters = []
    point_events = events.loc[events["event_type"] != "broad_shift"].dropna(
        subset=["summit"]
    )
    for _, chromosome_events in point_events.groupby("chrom", sort=False):
        ordered = chromosome_events.sort_values("summit", kind="stable")
        group = []
        previous_summit = None
        for index, row in ordered.iterrows():
            if (
                previous_summit is not None
                and row["summit"] - previous_summit
                > settings.spatial_cluster_tolerance_bp
            ):
                clusters.append(_cluster_point_rows(ordered.loc[group], settings))
                group = []
            group.append(index)
            previous_summit = row["summit"]
        if group:
            clusters.append(_cluster_point_rows(ordered.loc[group], settings))

    broad_events = events.loc[events["event_type"] == "broad_shift"]
    for _, chromosome_events in broad_events.groupby("chrom", sort=False):
        ordered = chromosome_events.sort_values("event_start", kind="stable")
        group = []
        current_end = None
        for index, row in ordered.iterrows():
            if current_end is not None and row["event_start"] >= current_end:
                clusters.append(_cluster_broad_rows(ordered.loc[group], settings))
                group = []
                current_end = None
            group.append(index)
            current_end = max(current_end, row["event_end"]) if current_end is not None else row["event_end"]
        if group:
            clusters.append(_cluster_broad_rows(ordered.loc[group], settings))

    return pd.DataFrame(clusters, columns=_CLUSTER_COLUMNS).sort_values(
        ["chrom", "event_start", "event_type"], kind="stable"
    ).reset_index(drop=True)

