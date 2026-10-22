#!/usr/bin/env python3
"""Call multiscale IS-change loci within already selected candidate regions."""

import argparse
import os
from pathlib import Path
import tempfile

import cooler
import numpy as np
import pandas as pd
import yaml

from cworld_is_config import parse_cworld_is_settings
from is_change_config import parse_is_change_settings
from is_change_detection import (
    build_consensus_track,
    classify_scale_events,
    cluster_multiscale_events,
    local_consensus_extrema,
    refine_point_clusters,
    recover_residual_events,
    residual_extrema,
    score_scale,
)


REFINEMENT_COLUMNS = [
    "discovery_mode", "original_summit", "refinement_status",
    "refinement_direction", "refinement_support_scale_bands",
    "insulation_direction",
]
# lower IS means a stronger boundary, so a valley in B is a strengthened
# boundary and a peak in B is a weakened one.
INSULATION_DIRECTION = {"valley": "stronger_in_B", "peak": "weaker_in_B"}
PRIMARY_COLUMNS = [
    "comparison", "region_id", "chrom", "event_start", "event_end", "summit",
    "is_change_type", "confidence_tier", "support_scale_bands", "support_scales",
    "max_amplitude_percentile", "max_transition_percentile", "max_pulse_percentile",
    "broad_shift_length_bp", "mask_status", "coverage_state_summary", "rank_within_region",
] + REFINEMENT_COLUMNS
MASK_COLUMNS = ["chrom", "start", "end", "reason", "scale_bp"]
TRACK_VALUE_COLUMNS = ["start", "IS"]


def index_is_track(track):
    return {
        (chrom, int(scale)): group[TRACK_VALUE_COLUMNS].copy().reset_index(drop=True)
        for (chrom, scale), group in track.groupby(
            ["chrom", "IS_distance"], sort=False, observed=True
        )
    }


def track_values(track_index, chrom, scale):
    return track_index.get(
        (chrom, scale), pd.DataFrame(columns=TRACK_VALUE_COLUMNS)
    )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-regions", required=True)
    parser.add_argument("--is-a", required=True)
    parser.add_argument("--is-b", required=True)
    parser.add_argument("--cool-a", required=True, nargs="+")
    parser.add_argument("--cool-b", required=True, nargs="+")
    parser.add_argument("--config", required=True)
    parser.add_argument("--out-primary", required=True)
    parser.add_argument("--out-preliminary", required=True)
    parser.add_argument("--out-masked", required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def log(message, path):
    print(message)
    log_path = Path(path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as handle:
        handle.write(f"{message}\n")


def _integer_column(frame, column, label):
    values = pd.to_numeric(frame[column], errors="raise")
    if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
        raise ValueError(f"{label} {column} must contain finite integers")
    return values.astype(int)


def read_regions(path, resolution_bp):
    regions = pd.read_csv(path, sep="\t")
    required = {"region_id", "chrom", "region_start", "region_end"}
    missing = required.difference(regions.columns)
    if missing:
        raise ValueError(f"candidate regions missing columns: {', '.join(sorted(missing))}")
    regions = regions[["region_id", "chrom", "region_start", "region_end"]].copy()
    regions["region_start"] = _integer_column(regions, "region_start", "candidate regions")
    regions["region_end"] = _integer_column(regions, "region_end", "candidate regions")
    if regions["region_id"].duplicated().any():
        raise ValueError("candidate region_id values must be unique")
    if (regions["region_start"] < 0).any() or (regions["region_end"] <= regions["region_start"]).any():
        raise ValueError("candidate regions must have non-negative, non-empty intervals")
    if ((regions["region_start"] % resolution_bp != 0) | (regions["region_end"] % resolution_bp != 0)).any():
        raise ValueError(f"candidate regions must align to bin size {resolution_bp}")
    return regions.sort_values(["chrom", "region_start", "region_end"], kind="stable")


def read_is_track(path, settings):
    track = pd.read_csv(path, sep="\t")
    required = {"chrom", "viewpoint", "IS_distance", "IS"}
    missing = required.difference(track.columns)
    if missing:
        raise ValueError(f"IS track {path} missing columns: {', '.join(sorted(missing))}")
    track = track[["chrom", "viewpoint", "IS_distance", "IS"]].copy()
    track["start"] = _integer_column(track, "viewpoint", f"IS track {path}")
    track["IS_distance"] = _integer_column(track, "IS_distance", f"IS track {path}")
    if (track["start"] < 0).any() or (track["start"] % settings.resolution_bp != 0).any():
        raise ValueError(f"IS track {path} viewpoints must align to bin size {settings.resolution_bp}")
    unknown = set(track["IS_distance"]).difference(settings.scales_bp)
    if unknown:
        raise ValueError(f"IS track {path} has unconfigured scales: {', '.join(map(str, sorted(unknown)))}")
    missing_scales = set(settings.scales_bp).difference(track["IS_distance"])
    if missing_scales:
        raise ValueError(f"IS track {path} missing configured scales: {', '.join(map(str, sorted(missing_scales)))}")
    if track.duplicated(["chrom", "start", "IS_distance"]).any():
        raise ValueError(f"IS track {path} has duplicate chrom/viewpoint/scale rows")
    track["IS"] = pd.to_numeric(track["IS"], errors="coerce")
    return track[["chrom", "start", "IS_distance", "IS"]]


def read_balanced_bins(paths, resolution_bp):
    combined = []
    for path in paths:
        bins = cooler.Cooler(str(path)).bins()[:]
        required = {"chrom", "start", "end", "weight"}
        missing = required.difference(bins.columns)
        if missing:
            raise ValueError(f"cooler {path} missing bins columns: {', '.join(sorted(missing))}")
        bins = bins[["chrom", "start", "end", "weight"]].copy()
        bins["start"] = _integer_column(bins, "start", f"cooler {path}")
        bins["end"] = _integer_column(bins, "end", f"cooler {path}")
        widths = bins["end"] - bins["start"]
        terminal_ends = bins.groupby("chrom", observed=True)["end"].transform("max")
        terminal_partial = (widths > 0) & (widths < resolution_bp) & (bins["end"] == terminal_ends)
        valid_widths = (widths == resolution_bp) | terminal_partial
        if not valid_widths.all() or (bins["start"] % resolution_bp != 0).any():
            raise ValueError(f"cooler {path} bins must be {resolution_bp} bp or terminal partial bins")
        if bins.duplicated(["chrom", "start"]).any():
            raise ValueError(f"cooler {path} has duplicate bins")
        combined.append(bins[["chrom", "start"]].assign(available=np.isfinite(bins["weight"])))
    if not combined:
        raise ValueError("at least one cooler is required")
    combined = pd.concat(combined, ignore_index=True)
    if combined.duplicated(["chrom", "start"]).any():
        raise ValueError("cooler inputs have duplicate chrom/start bins")
    return combined.sort_values(["chrom", "start"], kind="stable").reset_index(drop=True)


def context_frame(chrom, core_start, core_end, settings, track_a, track_b, bins_a, bins_b):
    start = max(0, ((core_start - settings.context_bp) // settings.resolution_bp) * settings.resolution_bp)
    end = ((core_end + settings.context_bp + settings.resolution_bp - 1) // settings.resolution_bp) * settings.resolution_bp
    grid = pd.DataFrame({"chrom": chrom, "start": np.arange(start, end, settings.resolution_bp)})
    base = grid.merge(bins_a.rename(columns={"available": "balanced_a"}), on=["chrom", "start"], how="left")
    base = base.merge(bins_b.rename(columns={"available": "balanced_b"}), on=["chrom", "start"], how="left")
    base[["balanced_a", "balanced_b"]] = base[["balanced_a", "balanced_b"]].fillna(False).astype(bool)
    return base


def _event_rows_for_region(region, settings, track_a, track_b, bins_a, bins_b):
    chrom, core_start, core_end = region[["chrom", "region_start", "region_end"]]
    base = context_frame(chrom, core_start, core_end, settings, track_a, track_b, bins_a, bins_b)
    events, masks, core_score_tracks = [], [], {}
    for scale in settings.scales_bp:
        frame = base.copy()
        for track_index, column in ((track_a, "is_a"), (track_b, "is_b")):
            values = track_values(track_index, chrom, scale)
            frame = frame.merge(
                values.rename(columns={"IS": column}), on="start", how="left"
            )
        frame["is_a"] = frame["is_a"].astype(float)
        frame["is_b"] = frame["is_b"].astype(float)
        scored = score_scale(chrom, frame["start"], frame["is_a"], frame["is_b"], frame["balanced_a"], frame["balanced_b"], scale, settings)
        core = frame[frame["start"].between(core_start, core_end - settings.resolution_bp)].reset_index(drop=True)
        core_scores = scored[scored["coordinate"].between(core_start, core_end - settings.resolution_bp)].reset_index(drop=True)
        core_score_tracks[scale] = core_scores
        for index, row in core.iterrows():
            reasons = []
            if not row["balanced_a"] and not row["balanced_b"]:
                reasons.append("BOTH_EMPTY")
            elif not row["balanced_a"]:
                reasons.append("A_EMPTY")
            elif not row["balanced_b"]:
                reasons.append("B_EMPTY")
            if not np.isfinite(row["is_a"]):
                reasons.append("IS_NA_A")
            if not np.isfinite(row["is_b"]):
                reasons.append("IS_NA_B")
            score = core_scores.iloc[index]
            if not reasons and (not score["transition_evaluable"] or not score["pulse_evaluable"]):
                reasons.append("FLANK_COVERAGE")
            masks.extend(
                (chrom, int(row["start"]), int(row["start"] + settings.resolution_bp), reason, int(scale))
                for reason in reasons
            )
        scale_events = classify_scale_events(scored, settings)
        points = scale_events[scale_events["event_type"] != "broad_shift"]
        points = points[points["summit"].between(core_start, core_end - settings.resolution_bp)]
        broad = scale_events[scale_events["event_type"] == "broad_shift"].copy()
        broad = broad[(broad["event_end"] > core_start) & (broad["event_start"] < core_end)]
        broad["event_start"] = broad["event_start"].clip(lower=core_start)
        broad["event_end"] = broad["event_end"].clip(upper=core_end)
        events.append(pd.concat((points, broad), ignore_index=True))
    events = pd.concat(events, ignore_index=True) if events else pd.DataFrame()
    clusters = cluster_multiscale_events(events, settings)
    return clusters, masks, core_score_tracks


def _coverage_summary(region, mask_rows):
    reasons = sorted({row[3] for row in mask_rows})
    return "clean" if not reasons else ",".join(reasons)


def output_rows(regions, settings, track_a, track_b, bins_a, bins_b, comparison):
    primary, preliminary, masks, analyses = [], [], [], []
    for _, region in regions.iterrows():
        clusters, region_masks, core_score_tracks = _event_rows_for_region(
            region, settings, track_a, track_b, bins_a, bins_b
        )
        masks.extend(region_masks)
        analyses.append((region, clusters, region_masks, core_score_tracks))

    residuals_by_region, score_tracks_by_region, refined_analyses = {}, {}, []
    for region, clusters, region_masks, core_score_tracks in analyses:
        region_id = region["region_id"]
        initial_summits = set(
            clusters.loc[
                clusters["event_type"].ne("broad_shift"), "summit"
            ].dropna().astype(int)
        )
        if core_score_tracks:
            consensus = build_consensus_track(core_score_tracks, settings)
            extrema = local_consensus_extrema(consensus, settings)
            clusters = refine_point_clusters(
                clusters, consensus, extrema, core_score_tracks, settings
            )
            final_summits = set(
                clusters.loc[
                    clusters["event_type"].ne("broad_shift"), "summit"
                ].dropna().astype(int)
            )
            residuals_by_region[region_id] = residual_extrema(
                consensus, initial_summits | final_summits, settings
            )
            score_tracks_by_region[region_id] = core_score_tracks
        else:
            clusters = clusters.copy()
            broad = clusters["event_type"].eq("broad_shift")
            clusters["discovery_mode"] = "original"
            clusters["original_summit"] = clusters["summit"].where(~broad, np.nan)
            clusters["refinement_status"] = "unchanged"
            clusters["refinement_direction"] = np.nan
            clusters["refinement_support_scale_bands"] = ""
            residuals_by_region[region_id] = pd.DataFrame()
        refined_analyses.append((region, clusters, region_masks))

    recovered_by_region = recover_residual_events(
        residuals_by_region, score_tracks_by_region, settings
    )
    for region, clusters, region_masks in refined_analyses:
        recovered = recovered_by_region[region["region_id"]]
        if not recovered.empty:
            clusters = pd.concat((clusters, recovered), ignore_index=True)
        summary = _coverage_summary(region, region_masks)
        for _, cluster in clusters.iterrows():
            row = {
                "comparison": comparison, "region_id": region["region_id"], "chrom": cluster["chrom"],
                "event_start": int(cluster["event_start"]), "event_end": int(cluster["event_end"]),
                "summit": cluster["summit"], "is_change_type": cluster["event_type"],
                "confidence_tier": cluster["confidence_tier"], "support_scale_bands": cluster["support_scale_bands"],
                "support_scales": cluster["support_scales"],
                "max_amplitude_percentile": cluster["max_amplitude_percentile"],
                "max_transition_percentile": cluster["max_transition_percentile"],
                "max_pulse_percentile": cluster["max_pulse_percentile"],
                "broad_shift_length_bp": int(cluster["event_end"] - cluster["event_start"]) if cluster["event_type"] == "broad_shift" else 0,
                "mask_status": "clean", "coverage_state_summary": summary, "rank_within_region": 0,
                "discovery_mode": cluster["discovery_mode"], "original_summit": cluster["original_summit"],
                "refinement_status": cluster["refinement_status"],
                "refinement_direction": cluster["refinement_direction"],
                "refinement_support_scale_bands": cluster["refinement_support_scale_bands"],
                "insulation_direction": INSULATION_DIRECTION.get(
                    cluster["refinement_direction"], ""
                ),
            }
            (primary if row["confidence_tier"] == "primary" else preliminary).append(row)
    def finish(rows):
        frame = pd.DataFrame(rows, columns=PRIMARY_COLUMNS)
        frame["_support_band_count"] = frame["support_scale_bands"].fillna("").map(
            lambda bands: len(set(filter(None, str(bands).split(","))))
        )
        frame = frame.sort_values(
            [
                "region_id", "_support_band_count", "max_transition_percentile",
                "max_pulse_percentile", "max_amplitude_percentile", "chrom",
                "event_start", "event_end", "is_change_type",
            ],
            ascending=[True, False, False, False, False, True, True, True, True],
            kind="stable",
            na_position="last",
        ).reset_index(drop=True)
        if not frame.empty:
            frame["rank_within_region"] = frame.groupby("region_id", sort=False).cumcount() + 1
        return frame.drop(columns="_support_band_count")
    masked = pd.DataFrame(sorted(set(masks)), columns=MASK_COLUMNS)
    return finish(primary), finish(preliminary), masked

def atomic_write(primary, preliminary, masked, paths):
    tables = ((primary, paths[0], True), (preliminary, paths[1], True), (masked, paths[2], False))
    temporary = []
    try:
        for table, target, header in tables:
            target = Path(target)
            target.parent.mkdir(parents=True, exist_ok=True)
            descriptor, name = tempfile.mkstemp(dir=target.parent, prefix=f".{target.name}.", suffix=".tmp")
            os.close(descriptor)
            temporary.append((name, target))
            table.to_csv(name, sep="\t", index=False, header=header)
        for name, target in temporary:
            os.replace(name, target)
    except BaseException:
        for name, _ in temporary:
            Path(name).unlink(missing_ok=True)
        raise


def main():
    args = parse_args()
    with open(args.config) as handle:
        config = yaml.safe_load(handle)
    settings = parse_is_change_settings(config, parse_cworld_is_settings(config))
    regions = read_regions(args.candidate_regions, settings.resolution_bp)
    track_a = index_is_track(read_is_track(args.is_a, settings))
    track_b = index_is_track(read_is_track(args.is_b, settings))
    bins_a = read_balanced_bins(args.cool_a, settings.resolution_bp)
    bins_b = read_balanced_bins(args.cool_b, settings.resolution_bp)
    comparison = f"{Path(args.is_a).stem}__vs__{Path(args.is_b).stem}"
    primary, preliminary, masked = output_rows(regions, settings, track_a, track_b, bins_a, bins_b, comparison)
    if list(primary.columns) != PRIMARY_COLUMNS or list(preliminary.columns) != PRIMARY_COLUMNS:
        raise AssertionError("invalid IS-change output schema")
    if list(masked.columns) != MASK_COLUMNS:
        raise AssertionError("invalid mask output schema")
    atomic_write(primary, preliminary, masked, (args.out_primary, args.out_preliminary, args.out_masked))
    log(f"Wrote {len(primary)} primary, {len(preliminary)} preliminary, and {len(masked)} masked rows", args.log)


if __name__ == "__main__":
    main()
