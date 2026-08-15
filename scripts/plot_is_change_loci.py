#!/usr/bin/env python3
"""Render diagnostic multiscale IS plots for the highest-ranked primary loci."""

import argparse
import os
from pathlib import Path
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from cworld_is_config import parse_cworld_is_settings
from is_change_config import parse_is_change_settings


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", required=True)
    parser.add_argument("--candidate-regions", required=True)
    parser.add_argument("--is-a", required=True)
    parser.add_argument("--is-b", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--manifest", required=True)
    return parser.parse_args()


def read_track(path, settings):
    track = pd.read_csv(path, sep="\t")
    required = {"chrom", "viewpoint", "IS_distance", "IS"}
    missing = required.difference(track.columns)
    if missing:
        raise ValueError(f"IS track {path} missing columns: {', '.join(sorted(missing))}")
    track = track[["chrom", "viewpoint", "IS_distance", "IS"]].copy()
    track["start"] = pd.to_numeric(track["viewpoint"], errors="raise").astype(int)
    track["IS_distance"] = pd.to_numeric(track["IS_distance"], errors="raise").astype(int)
    if (track["start"] % settings.resolution_bp != 0).any():
        raise ValueError(f"IS track {path} viewpoints must align to bin size {settings.resolution_bp}")
    if set(settings.scales_bp).difference(track["IS_distance"]):
        raise ValueError(f"IS track {path} is missing configured scales")
    if track.duplicated(["chrom", "start", "IS_distance"]).any():
        raise ValueError(f"IS track {path} has duplicate chrom/viewpoint/scale rows")
    track["IS"] = pd.to_numeric(track["IS"], errors="coerce")
    return track[["chrom", "start", "IS_distance", "IS"]]


def plot_region(region, events, track_a, track_b, settings, output):
    start = max(0, int(region.region_start) - settings.context_bp)
    end = int(region.region_end) + settings.context_bp
    figure, axes = plt.subplots(len(settings.scales_bp), 1, figsize=(13, 2.1 * len(settings.scales_bp)), sharex=True)
    for axis, scale in zip(np.atleast_1d(axes), settings.scales_bp):
        a = track_a[(track_a.chrom == region.chrom) & (track_a.IS_distance == scale) & track_a.start.between(start, end - settings.resolution_bp)]
        b = track_b[(track_b.chrom == region.chrom) & (track_b.IS_distance == scale) & track_b.start.between(start, end - settings.resolution_bp)]
        frame = pd.DataFrame({"start": np.arange(start, end, settings.resolution_bp)})
        frame = frame.merge(a[["start", "IS"]].rename(columns={"IS": "a"}), on="start", how="left")
        frame = frame.merge(b[["start", "IS"]].rename(columns={"IS": "b"}), on="start", how="left")
        unavailable = ~(np.isfinite(frame.a) & np.isfinite(frame.b))
        for point in frame.loc[unavailable, "start"]:
            axis.axvspan(point, point + settings.resolution_bp, color="0.85", zorder=0)
        axis.axvspan(region.region_start, region.region_end, color="#dbeafe", zorder=0)
        axis.plot(frame.start, frame.a, color="#2563eb", lw=0.8, label="A")
        axis.plot(frame.start, frame.b, color="#16a34a", lw=0.8, label="B")
        axis.plot(frame.start, frame.b - frame.a, color="#6b21a8", lw=0.8, label="B-A")
        for event in events.itertuples():
            if event.is_change_type == "broad_shift":
                axis.axvspan(event.event_start, event.event_end, color="#dc2626", alpha=0.20)
            elif pd.notna(event.summit):
                axis.axvline(event.summit, color="#dc2626", lw=1)
        axis.set_ylabel(f"{scale // 1000} kb")
        axis.grid(axis="y", alpha=0.2)
    axes[0].legend(loc="upper right", ncol=3, fontsize=8)
    axes[-1].set_xlabel(f"{region.chrom} coordinate (bp)")
    figure.tight_layout()
    figure.savefig(output, dpi=140)
    plt.close(figure)


def atomic_manifest(frame, path):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(dir=target.parent, prefix=f".{target.name}.", suffix=".tmp")
    os.close(descriptor)
    try:
        frame.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, target)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def main():
    args = parse_args()
    with open(args.config) as handle:
        config = yaml.safe_load(handle)
    settings = parse_is_change_settings(config, parse_cworld_is_settings(config))
    primary = pd.read_csv(args.primary, sep="\t")
    regions = pd.read_csv(args.candidate_regions, sep="\t")
    required = {"region_id", "chrom", "event_start", "event_end", "summit", "is_change_type"}
    missing = required.difference(primary.columns)
    if missing:
        raise ValueError(f"primary loci missing columns: {', '.join(sorted(missing))}")
    required_regions = {"region_id", "chrom", "region_start", "region_end"}
    if required_regions.difference(regions.columns):
        raise ValueError("candidate regions missing required coordinate columns")
    track_a, track_b = read_track(args.is_a, settings), read_track(args.is_b, settings)
    sort_columns = [column for column in ("rank_within_region", "max_amplitude_percentile", "event_start") if column in primary]
    top = primary.sort_values(
        sort_columns, ascending=[True, False, True][:len(sort_columns)], kind="stable"
    ).head(settings.diagnostic_plot_top_n)
    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for region_id in top.region_id.drop_duplicates():
        region = regions.loc[regions.region_id == region_id]
        if len(region) != 1:
            raise ValueError(f"primary locus references unknown or duplicate region_id {region_id}")
        region = region.iloc[0]
        events = primary[primary.region_id == region_id]
        comparison = str(events.iloc[0].get("comparison", "comparison"))
        output = output_dir / f"{comparison}_region_{region_id}.png"
        plot_region(region, events, track_a, track_b, settings, output)
        manifest.append({"region_id": region_id, "chrom": region.chrom, "plot_path": str(output.resolve())})
    atomic_manifest(pd.DataFrame(manifest, columns=["region_id", "chrom", "plot_path"]), args.manifest)


if __name__ == "__main__":
    main()
