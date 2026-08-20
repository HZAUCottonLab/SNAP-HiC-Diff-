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

CONSENSUS_COLUMNS = ["chrom", "coordinate", "consensus_z", "valid", "segment_id"]
EXTREMUM_COLUMNS = [
    "chrom",
    "coordinate",
    "plateau_start",
    "plateau_end",
    "direction",
    "segment_id",
    "prominence",
    "left_flank_coordinate",
    "right_flank_coordinate",
]


def build_consensus_track(score_tracks: dict[int, pd.DataFrame], settings) -> pd.DataFrame:
    """Combine aligned score tracks into a common-valid signed-z consensus."""
    scales = tuple(settings.scales_bp)
    missing_scales = [scale for scale in scales if scale not in score_tracks]
    if missing_scales:
        raise ValueError(f"score tracks missing configured scales: {missing_scales}")

    tracks = [score_tracks[scale].reset_index(drop=True) for scale in scales]
    required = {"chrom", "coordinate", "valid", "z"}
    if any(required.difference(track.columns) for track in tracks):
        raise ValueError("score tracks must contain chrom, coordinate, valid, and z")
    reference = tracks[0]
    if reference["chrom"].nunique(dropna=False) != 1:
        raise ValueError("score tracks must describe one chromosome")
    if any(
        not reference["chrom"].equals(track["chrom"])
        or not reference["coordinate"].equals(track["coordinate"])
        for track in tracks[1:]
    ):
        raise ValueError("score tracks must share chromosome and coordinate series")

    coordinates = reference["coordinate"].to_numpy()
    if np.any(coordinates % settings.resolution_bp):
        raise ValueError("score track coordinates must be aligned to resolution_bp")
    if len(np.unique(coordinates)) != len(coordinates):
        raise ValueError("score track coordinates must be unique")
    if np.any(coordinates[1:] <= coordinates[:-1]):
        raise ValueError("score track coordinates must be strictly increasing")

    z_values = [track["z"].to_numpy(dtype=float) for track in tracks]
    valid = np.logical_and.reduce(
        [
            track["valid"].fillna(False).to_numpy(dtype=bool) & np.isfinite(z)
            for track, z in zip(tracks, z_values)
        ]
    )
    consensus_z = np.full(len(reference), np.nan, dtype=float)
    if valid.any():
        consensus_z[valid] = np.median(np.asarray(z_values)[:, valid], axis=0)

    segment_id = np.full(len(reference), np.nan)
    current_segment = -1
    for index in np.flatnonzero(valid):
        if (
            index == 0
            or not valid[index - 1]
            or coordinates[index] - coordinates[index - 1] != settings.resolution_bp
        ):
            current_segment += 1
        segment_id[index] = current_segment

    return pd.DataFrame(
        {
            "chrom": reference["chrom"],
            "coordinate": coordinates,
            "consensus_z": consensus_z,
            "valid": valid,
            "segment_id": segment_id,
        },
        columns=CONSENSUS_COLUMNS,
    )


def local_consensus_extrema(consensus: pd.DataFrame, settings) -> pd.DataFrame:
    """Extract one strict peak or valley for each equal-valued plateau."""
    required = set(CONSENSUS_COLUMNS)
    if required.difference(consensus.columns):
        raise ValueError("consensus track is missing required columns")

    track = consensus.reset_index(drop=True)
    coordinates = track["coordinate"].to_numpy()
    segments = track["segment_id"].to_numpy()
    z = track["consensus_z"].to_numpy(dtype=float)
    valid = (
        track["valid"].fillna(False).to_numpy(dtype=bool)
        & np.isfinite(z)
        & pd.notna(track["segment_id"]).to_numpy()
    )
    records = []
    start = 0
    while start < len(track):
        if not valid[start]:
            start += 1
            continue
        end = start
        while (
            end + 1 < len(track)
            and valid[end + 1]
            and segments[end + 1] == segments[start]
            and coordinates[end + 1] - coordinates[end] == settings.resolution_bp
            and np.isclose(z[end + 1], z[end], rtol=0, atol=1e-12)
        ):
            end += 1
        left = start - 1
        right = end + 1
        if (
            left >= 0
            and right < len(track)
            and valid[left]
            and valid[right]
            and segments[left] == segments[start] == segments[right]
            and coordinates[start] - coordinates[left] == settings.resolution_bp
            and coordinates[right] - coordinates[end] == settings.resolution_bp
        ):
            direction = None
            if z[start] > z[left] and z[start] > z[right]:
                direction = "peak"
            elif z[start] < z[left] and z[start] < z[right]:
                direction = "valley"
            if direction:
                records.append(
                    {
                        "chrom": track.at[start, "chrom"],
                        "coordinate": coordinates[start],
                        "plateau_start": coordinates[start],
                        "plateau_end": coordinates[end],
                        "direction": direction,
                        "segment_id": segments[start],
                        "prominence": np.nan,
                        "left_flank_coordinate": np.nan,
                        "right_flank_coordinate": np.nan,
                        "_z": z[start],
                    }
                )
        start = end + 1

    left_flanks = {}
    for record in records:
        opposite = "valley" if record["direction"] == "peak" else "peak"
        left = left_flanks.get((record["segment_id"], opposite))
        record["_left"] = left
        if left is not None:
            record["left_flank_coordinate"] = left["coordinate"]
        left_flanks[(record["segment_id"], record["direction"])] = record

    right_flanks = {}
    for record in reversed(records):
        opposite = "valley" if record["direction"] == "peak" else "peak"
        right = right_flanks.get((record["segment_id"], opposite))
        if right is not None:
            record["right_flank_coordinate"] = right["coordinate"]
        left = record["_left"]
        if left is not None and right is not None:
            drops = (
                (record["_z"] - left["_z"], record["_z"] - right["_z"])
                if record["direction"] == "peak"
                else (left["_z"] - record["_z"], right["_z"] - record["_z"])
            )
            record["prominence"] = min(drops)
        right_flanks[(record["segment_id"], record["direction"])] = record

    return pd.DataFrame(records, columns=EXTREMUM_COLUMNS)


def directional_support(
    score_tracks: dict[int, pd.DataFrame],
    coordinate: int,
    direction: str,
    settings,
    require_high_amplitude: bool = False,
) -> tuple[list[str], list[int]]:
    """Return configured scale bands and scales supporting one signed extremum."""
    if direction not in {"peak", "valley"}:
        raise ValueError("direction must be peak or valley")
    scales = tuple(settings.scales_bp)
    missing_scales = [scale for scale in scales if scale not in score_tracks]
    if missing_scales:
        raise ValueError(f"score tracks missing configured scales: {missing_scales}")

    supporting_scales = []
    for scale in scales:
        rows = score_tracks[scale].loc[
            score_tracks[scale]["coordinate"] == coordinate
        ]
        if rows.empty:
            return [], []
        row = rows.iloc[0]
        z = row["z"]
        if (
            not pd.notna(row["valid"])
            or not bool(row["valid"])
            or not np.isfinite(z)
            or (direction == "peak" and z <= 0)
            or (direction == "valley" and z >= 0)
        ):
            continue
        if require_high_amplitude and (
            "amplitude_percentile" not in row
            or not np.isfinite(row["amplitude_percentile"])
            or row["amplitude_percentile"] < settings.high_percentile
        ):
            continue
        supporting_scales.append(scale)

    support_bands = [
        band
        for band, band_scales in settings.scale_bands.items()
        if any(scale in band_scales for scale in supporting_scales)
    ]
    return support_bands, sorted(supporting_scales)


_REFINEMENT_COLUMNS = [
    "discovery_mode",
    "original_summit",
    "refinement_status",
    "refinement_direction",
    "refinement_support_scale_bands",
]


def select_nearest_extremum(
    original_summit, direction, extrema, consensus, resolution_bp
):
    """Select the closest member of either nearest same-direction plateau."""
    original = consensus.loc[consensus["coordinate"] == original_summit]
    if original.empty or not bool(original.iloc[0]["valid"]):
        return None
    segment_id = original.iloc[0]["segment_id"]
    if pd.isna(segment_id):
        return None

    candidates = extrema.loc[
        (extrema["direction"] == direction)
        & (extrema["segment_id"] == segment_id)
    ]
    if candidates.empty:
        return None
    containing = candidates.loc[
        candidates["plateau_start"].le(original_summit)
        & candidates["plateau_end"].ge(original_summit)
    ]
    if containing.empty:
        left = candidates.loc[candidates["plateau_end"] < original_summit]
        right = candidates.loc[candidates["plateau_start"] > original_summit]
        candidates = pd.concat(
            (
                left.loc[[left["plateau_end"].idxmax()]] if not left.empty else left,
                right.loc[[right["plateau_start"].idxmin()]] if not right.empty else right,
            )
        )
    else:
        candidates = containing

    choices = []
    for _, extremum in candidates.iterrows():
        members = np.arange(
            extremum["plateau_start"],
            extremum["plateau_end"] + resolution_bp,
            resolution_bp,
        )
        coordinate = min(members, key=lambda member: (abs(member - original_summit), member))
        z = consensus.loc[consensus["coordinate"] == coordinate, "consensus_z"].iloc[0]
        choices.append((coordinate, z))
    return min(
        choices,
        key=lambda choice: (abs(choice[0] - original_summit), -abs(choice[1]), choice[0]),
    )[0]


def refine_point_clusters(clusters, consensus, extrema, score_tracks, settings):
    """Move point clusters to supported consensus extrema and record provenance."""
    columns = list(clusters.columns) + [
        column for column in _REFINEMENT_COLUMNS if column not in clusters
    ]
    if clusters.empty:
        return clusters.reindex(columns=columns)

    records = []
    for _, row in clusters.iterrows():
        record = row.to_dict()
        if row["event_type"] == "broad_shift":
            record.update(
                discovery_mode="original",
                original_summit=np.nan,
                refinement_status="unchanged",
                refinement_direction=np.nan,
                refinement_support_scale_bands="",
            )
            records.append(record)
            continue

        original_summit = int(row["summit"])
        record.update(
            discovery_mode="original",
            original_summit=original_summit,
            refinement_status="ambiguous",
            refinement_direction=np.nan,
            refinement_support_scale_bands="",
        )
        original = consensus.loc[consensus["coordinate"] == original_summit]
        if (
            original.empty
            or not bool(original.iloc[0]["valid"])
            or not np.isfinite(original.iloc[0]["consensus_z"])
            or original.iloc[0]["consensus_z"] == 0
        ):
            records.append(record)
            continue

        direction = "peak" if original.iloc[0]["consensus_z"] > 0 else "valley"
        candidate = select_nearest_extremum(
            original_summit, direction, extrema, consensus, settings.resolution_bp
        )
        final_summit = original_summit
        if candidate is not None:
            bands, _ = directional_support(score_tracks, candidate, direction, settings)
            if len(bands) >= settings.min_primary_scale_bands:
                final_summit = int(candidate)
        bands, _ = directional_support(score_tracks, final_summit, direction, settings)
        record.update(
            refinement_status="refined" if final_summit != original_summit else "unchanged",
            refinement_direction=direction,
            refinement_support_scale_bands=",".join(bands),
        )
        if final_summit != original_summit:
            record.update(
                discovery_mode="refined",
                summit=final_summit,
                event_start=final_summit,
                event_end=final_summit + settings.resolution_bp,
            )
        records.append(record)

    refined = pd.DataFrame(records, columns=columns)
    collisions = {}
    for index, row in refined.iterrows():
        if row["event_type"] == "broad_shift" or pd.isna(row["refinement_direction"]):
            continue
        plateau = extrema.loc[
            (extrema["direction"] == row["refinement_direction"])
            & extrema["plateau_start"].le(row["summit"])
            & extrema["plateau_end"].ge(row["summit"])
        ]
        if not plateau.empty:
            extremum = plateau.iloc[0]
            key = (
                row["chrom"],
                row["refinement_direction"],
                extremum["segment_id"],
                extremum["plateau_start"],
                extremum["plateau_end"],
            )
            collisions.setdefault(key, []).append(index)

    keep = set(refined.index)
    for indices in collisions.values():
        if len(indices) < 2 or not (refined.loc[indices, "refinement_status"] == "refined").any():
            continue
        winner = min(indices, key=lambda index: _evidence_key(refined.loc[index]))
        keep.difference_update(set(indices) - {winner})
    return refined.loc[sorted(keep)].reset_index(drop=True)


def _evidence_key(row):
    band_count = len(set(filter(None, str(row["support_scale_bands"]).split(","))))
    scores = [
        row["max_transition_percentile"],
        row["max_pulse_percentile"],
        row["max_amplitude_percentile"],
    ]
    return (
        -band_count,
        *(-float(score) if np.isfinite(score) else np.inf for score in scores),
        row["original_summit"] if pd.notna(row["original_summit"]) else row["summit"],
    )


def residual_extrema(consensus: pd.DataFrame, occupied_coordinates: set[int], settings) -> pd.DataFrame:
    """Return finite-prominence consensus extrema not occupied by point calls."""
    residuals = local_consensus_extrema(consensus, settings)
    prominence = pd.to_numeric(residuals["prominence"], errors="coerce")
    eligible = np.isfinite(prominence) & prominence.ge(0)
    for coordinate in occupied_coordinates:
        eligible &= ~(
            residuals["plateau_start"].le(coordinate)
            & residuals["plateau_end"].ge(coordinate)
        )
    residuals = residuals.loc[eligible].copy()
    for index, extremum in residuals.iterrows():
        members = consensus.loc[
            consensus["coordinate"].between(
                extremum["plateau_start"], extremum["plateau_end"]
            )
        ]
        coordinate = min(
            members.itertuples(),
            key=lambda member: (-abs(member.consensus_z), member.coordinate),
        ).coordinate
        residuals.at[index, "coordinate"] = coordinate
    return residuals.reset_index(drop=True)


def recover_residual_events(
    residuals_by_region, score_tracks_by_region, settings
) -> dict[object, pd.DataFrame]:
    """Recover high-prominence, high-amplitude extrema across candidate cores."""
    recovered = {region_id: [] for region_id in residuals_by_region}
    pool = []
    for region_id, residuals in residuals_by_region.items():
        if residuals.empty:
            continue
        for residual in residuals.to_dict("records"):
            residual["_region_id"] = region_id
            pool.append(residual)
    if not pool:
        return {
            region_id: pd.DataFrame(columns=_CLUSTER_COLUMNS + _REFINEMENT_COLUMNS)
            for region_id in recovered
        }

    prominence_cutoff = pd.DataFrame(pool)["prominence"].quantile(
        settings.high_percentile
    )
    candidates = []
    for residual in pool:
        if residual["prominence"] < prominence_cutoff:
            continue
        region_id = residual["_region_id"]
        if region_id not in score_tracks_by_region:
            continue
        coordinate = int(residual["coordinate"])
        direction = residual["direction"]
        bands, scales = directional_support(
            score_tracks_by_region[region_id],
            coordinate,
            direction,
            settings,
            require_high_amplitude=True,
        )
        if len(bands) < settings.min_primary_scale_bands:
            continue
        evidence = {}
        for column in (
            "amplitude_percentile",
            "transition_percentile",
            "pulse_percentile",
        ):
            evidence[f"max_{column}"] = _finite_max(
                [
                    track.loc[track["coordinate"] == coordinate, column].iloc[0]
                    if column in track and not track.loc[track["coordinate"] == coordinate].empty
                    else np.nan
                    for scale, track in score_tracks_by_region[region_id].items()
                    if scale in scales
                ]
            )
        candidates.append(
            {
                "chrom": residual["chrom"],
                "event_type": "recovered_peak" if direction == "peak" else "recovered_valley",
                "event_start": coordinate,
                "event_end": coordinate + settings.resolution_bp,
                "summit": coordinate,
                "confidence_tier": "primary",
                "support_scale_bands": ",".join(bands),
                "support_scales": ",".join(map(str, scales)),
                "n_supporting_scales": len(scales),
                "max_amplitude_percentile": evidence["max_amplitude_percentile"],
                "max_transition_percentile": evidence["max_transition_percentile"],
                "max_pulse_percentile": evidence["max_pulse_percentile"],
                "discovery_mode": "recovered",
                "original_summit": np.nan,
                "refinement_status": "recovered",
                "refinement_direction": direction,
                "refinement_support_scale_bands": ",".join(bands),
                "_region_id": region_id,
                "_plateau_key": (
                    residual["chrom"],
                    direction,
                    residual["segment_id"],
                    residual["plateau_start"],
                    residual["plateau_end"],
                ),
            }
        )

    winners = {}
    for candidate in candidates:
        key = candidate["_plateau_key"]
        if key not in winners or (
            _evidence_key(candidate), str(candidate["_region_id"])
        ) < (
            _evidence_key(winners[key]), str(winners[key]["_region_id"])
        ):
            winners[key] = candidate
    for candidate in winners.values():
        recovered[candidate["_region_id"]].append(candidate)
    return {
        region_id: pd.DataFrame(records, columns=_CLUSTER_COLUMNS + _REFINEMENT_COLUMNS)
        for region_id, records in recovered.items()
    }


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
        transition_evaluable[center] = True
        if any(
            0 <= neighbor < n_bins and valid[neighbor] and z[center] * z[neighbor] > 0
            for neighbor in (center - 1, center + 1)
        ):
            pulse[center] = abs(z[center] - np.median(joint_flanks))
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


def _event_row(row, event_type, settings):
    return {
        "chrom": row["chrom"],
        "scale_bp": row["scale_bp"],
        "event_type": event_type,
        "event_start": row["coordinate"],
        "event_end": row["coordinate"] + settings.resolution_bp,
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
            events.append(_event_row(row, event_type, settings))

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

