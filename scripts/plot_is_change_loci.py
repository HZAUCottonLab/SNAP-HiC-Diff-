#!/usr/bin/env python3
"""Render diagnostic multiscale IS plots for the highest-ranked primary loci."""

import argparse
import os
from pathlib import Path
import tempfile
from typing import Mapping, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from matplotlib.image import AxesImage
import numpy as np
import pandas as pd
import yaml
import cooler

from cworld_is_config import parse_cworld_is_settings
from call_is_change_loci import read_balanced_bins
from is_change_config import ISChangeSettings, parse_is_change_settings


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", required=True)
    parser.add_argument("--candidate-regions", required=True)
    parser.add_argument("--is-a", required=True)
    parser.add_argument("--is-b", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--cool-a", required=True, nargs="+")
    parser.add_argument("--cool-b", required=True, nargs="+")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--manifest", required=True)
    return parser.parse_args()


def _integer_column(frame, column, label):
    values = pd.to_numeric(frame[column], errors="raise")
    if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
        raise ValueError(f"{label} {column} must contain finite integers")
    return values.astype(int)


def read_track(path, settings):
    track = pd.read_csv(path, sep="\t")
    required = {"chrom", "viewpoint", "IS_distance", "IS"}
    missing = required.difference(track.columns)
    if missing:
        raise ValueError("IS track {} missing columns: {}".format(path, ", ".join(sorted(missing))))
    track = track[["chrom", "viewpoint", "IS_distance", "IS"]].copy()
    track["start"] = _integer_column(track, "viewpoint", f"IS track {path}")
    track["IS_distance"] = _integer_column(track, "IS_distance", f"IS track {path}")
    if (track["start"] < 0).any() or (track["start"] % settings.resolution_bp != 0).any():
        raise ValueError(f"IS track {path} viewpoints must be non-negative and align to bin size {settings.resolution_bp}")
    unknown = set(track["IS_distance"]).difference(settings.scales_bp)
    if unknown:
        raise ValueError("IS track {} has unconfigured scales: {}".format(path, ", ".join(map(str, sorted(unknown)))))
    missing_scales = set(settings.scales_bp).difference(track["IS_distance"])
    if missing_scales:
        raise ValueError("IS track {} missing configured scales: {}".format(path, ", ".join(map(str, sorted(missing_scales)))))
    if track.duplicated(["chrom", "start", "IS_distance"]).any():
        raise ValueError(f"IS track {path} has duplicate chrom/viewpoint/scale rows")
    track["IS"] = pd.to_numeric(track["IS"], errors="coerce")
    return track[["chrom", "start", "IS_distance", "IS"]]

CONTACT_VMIN_PERCENTILE = 2.0
CONTACT_VMAX_PERCENTILE = 98.0
CONTACT_CMAP = "YlOrRd"
HEATMAP_HEIGHT_RATIO = 2.8
EVENT_COLOR = "red"
EVENT_ALPHA = 0.22
MAX_CONTACT_MATRIX_BINS = 2_000


def plot_interval(region: pd.Series, settings: ISChangeSettings) -> tuple[int, int]:
    """Return the common, resolution-aligned genomic extent for one region."""
    resolution = settings.resolution_bp
    context = settings.context_bp
    start = max(0, ((int(region.region_start) - context) // resolution) * resolution)
    end = ((int(region.region_end) + context + resolution - 1) // resolution) * resolution
    return start, end


def common_plot_interval(
    region: pd.Series,
    settings: ISChangeSettings,
    coolers_a: Mapping[str, cooler.Cooler],
    coolers_b: Mapping[str, cooler.Cooler],
    sample_a: str,
    sample_b: str,
) -> tuple[int, int]:
    """Clip one shared, resolution-aligned interval to both chromosome extents."""
    start, end = plot_interval(region, settings)
    chrom = region.chrom
    chrom_ends = []
    for coolers_by_chrom, sample_label in ((coolers_a, sample_a), (coolers_b, sample_b)):
        if chrom not in coolers_by_chrom:
            raise ValueError(f"{sample_label}: no iced cooler is available for {chrom}")
        chrom_ends.append(int(coolers_by_chrom[chrom].chromsizes[chrom]))
    common_start = 0
    common_end = (min(chrom_ends) // settings.resolution_bp) * settings.resolution_bp
    start = max(start, common_start)
    end = min(end, common_end)
    if end - start < settings.resolution_bp:
        raise ValueError(f"{chrom}: no common full plotting bin remains after chromosome-bound clipping")
    return start, end


def open_coolers_by_chrom(paths: Sequence[str], sample_label: str) -> dict[str, cooler.Cooler]:
    """Open each supplied iced cooler once and index it by the chromosomes it contains."""
    by_chrom: dict[str, cooler.Cooler] = {}
    for path in paths:
        if not path.endswith("_iced.cool"):
            raise ValueError(f"{sample_label}: requires _iced.cool cooler paths, got {path}")
        handle = cooler.Cooler(path)
        for chrom in handle.chromnames:
            if chrom in by_chrom:
                raise ValueError(f"{sample_label}: duplicate cooler for chromosome {chrom}")
            by_chrom[chrom] = handle
    if not by_chrom:
        raise ValueError(f"{sample_label}: no chromosomes found in supplied coolers")
    return by_chrom


def load_log_contact_matrix(
    coolers_by_chrom: Mapping[str, cooler.Cooler],
    chrom: str,
    start: int,
    end: int,
    sample_label: str,
) -> np.ma.MaskedArray:
    """Fetch one balanced cis matrix, retaining unavailable values as a mask."""
    if chrom not in coolers_by_chrom:
        raise ValueError(f"{sample_label}: no iced cooler is available for {chrom}")
    handle = coolers_by_chrom[chrom]
    if start < 0 or end <= start or end > int(handle.chromsizes[chrom]):
        raise ValueError(f"{sample_label}: invalid plotting interval {chrom}:{start}-{end}")
    binsize = int(handle.binsize)
    bin_count = (end + binsize - 1) // binsize - start // binsize
    if bin_count > MAX_CONTACT_MATRIX_BINS:
        raise ValueError(
            f"{sample_label}: plotting interval {chrom}:{start}-{end} exceeds dense plotting limit "
            f"of {MAX_CONTACT_MATRIX_BINS} bins"
        )
    raw = np.asarray(handle.matrix(balance=True).fetch(f"{chrom}:{start}-{end}"))
    finite = np.isfinite(raw)
    if np.any(raw[finite] < 0):
        raise ValueError(f"{sample_label}: balanced contacts contain finite negative values in {chrom}:{start}-{end}")
    if not finite.any():
        raise ValueError(f"{sample_label}: no finite balanced contacts in {chrom}:{start}-{end}")
    transformed = np.zeros(raw.shape, dtype=float)
    transformed[finite] = np.log1p(raw[finite])
    return np.ma.masked_array(transformed, mask=~finite)


def shared_contact_norm(matrix_a: np.ma.MaskedArray, matrix_b: np.ma.MaskedArray) -> Normalize:
    """Create the single linear colour normalisation used by both heatmaps."""
    values = np.concatenate((matrix_a.compressed(), matrix_b.compressed()))
    if values.size == 0:
        raise ValueError("No finite balanced contacts are available for either sample")
    vmin, vmax = np.percentile(values, (CONTACT_VMIN_PERCENTILE, CONTACT_VMAX_PERCENTILE))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin < 0.0:
        raise ValueError("Shared Hi-C colour scale has no finite non-negative range")
    return Normalize(vmin=float(vmin), vmax=float(vmin + 1.0 if vmax <= vmin else vmax), clip=True)


def scale_plot_frame(region, track_a, track_b, bins_a, bins_b, scale, settings, interval=None):
    start, end = plot_interval(region, settings) if interval is None else interval
    a = track_a[(track_a.chrom == region.chrom) & (track_a.IS_distance == scale) & track_a.start.between(start, end - settings.resolution_bp)]
    b = track_b[(track_b.chrom == region.chrom) & (track_b.IS_distance == scale) & track_b.start.between(start, end - settings.resolution_bp)]
    frame = pd.DataFrame({"chrom": region.chrom, "start": np.arange(start, end, settings.resolution_bp)})
    frame = frame.merge(a[["start", "IS"]].rename(columns={"IS": "a"}), on="start", how="left")
    frame = frame.merge(b[["start", "IS"]].rename(columns={"IS": "b"}), on="start", how="left")
    frame = frame.merge(
        bins_a.rename(columns={"available": "balanced_a"}), on=["chrom", "start"], how="left"
    )
    frame = frame.merge(
        bins_b.rename(columns={"available": "balanced_b"}), on=["chrom", "start"], how="left"
    )
    frame[["balanced_a", "balanced_b"]] = frame[["balanced_a", "balanced_b"]].fillna(False).astype(bool)
    frame["unavailable"] = ~(
        np.isfinite(frame["a"])
        & np.isfinite(frame["b"])
        & frame["balanced_a"]
        & frame["balanced_b"]
    )
    return frame


def draw_contact_heatmap(
    axis: Axes,
    matrix: np.ma.MaskedArray,
    start: int,
    end: int,
    sample_label: str,
    norm: Normalize,
    events: pd.DataFrame,
) -> AxesImage:
    cmap = plt.get_cmap(CONTACT_CMAP).copy()
    cmap.set_bad("0.85")
    image = axis.imshow(
        matrix,
        extent=(start, end, end, start),
        origin="upper",
        aspect="auto",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
    )
    axis.set_ylabel(f"{sample_label}\ncontact coordinate (bp)")
    for event in events.itertuples(index=False):
        if event.is_change_type == "broad_shift":
            axis.axvspan(event.event_start, event.event_end, color=EVENT_COLOR, alpha=EVENT_ALPHA, linewidth=0)
            axis.axhspan(event.event_start, event.event_end, color=EVENT_COLOR, alpha=EVENT_ALPHA, linewidth=0)
        elif pd.notna(event.summit):
            axis.axvline(event.summit, color=EVENT_COLOR, linewidth=1.0)
            axis.axhline(event.summit, color=EVENT_COLOR, linewidth=1.0)
    axis.set_xlim(start, end)
    axis.set_ylim(end, start)
    return image


def build_region_figure(
    region: pd.Series,
    events: pd.DataFrame,
    track_a: pd.DataFrame,
    track_b: pd.DataFrame,
    bins_a: pd.DataFrame,
    bins_b: pd.DataFrame,
    coolers_a: Mapping[str, cooler.Cooler],
    coolers_b: Mapping[str, cooler.Cooler],
    sample_a: str,
    sample_b: str,
    settings: ISChangeSettings,
) -> tuple[Figure, tuple[Axes, Axes], tuple[Axes, ...]]:
    start, end = common_plot_interval(
        region, settings, coolers_a, coolers_b, sample_a, sample_b
    )
    matrix_a = load_log_contact_matrix(coolers_a, region.chrom, start, end, sample_a)
    matrix_b = load_log_contact_matrix(coolers_b, region.chrom, start, end, sample_b)
    norm = shared_contact_norm(matrix_a, matrix_b)
    figure = plt.figure(
        figsize=(13, 2.1 * len(settings.scales_bp) + 2 * HEATMAP_HEIGHT_RATIO),
        constrained_layout=True,
    )
    grid = figure.add_gridspec(
        nrows=len(settings.scales_bp) + 2,
        ncols=1,
        height_ratios=[HEATMAP_HEIGHT_RATIO, HEATMAP_HEIGHT_RATIO, *([1.0] * len(settings.scales_bp))],
    )
    heatmap_a = figure.add_subplot(grid[0, 0])
    heatmap_b = figure.add_subplot(grid[1, 0], sharex=heatmap_a)
    is_axes = tuple(
        figure.add_subplot(grid[index + 2, 0], sharex=heatmap_a)
        for index in range(len(settings.scales_bp))
    )
    image = draw_contact_heatmap(heatmap_a, matrix_a, start, end, sample_a, norm, events)
    draw_contact_heatmap(heatmap_b, matrix_b, start, end, sample_b, norm, events)
    figure.colorbar(image, ax=(heatmap_a, heatmap_b), label="log1p(balanced contact; shared P2-P98)")
    for axis, scale in zip(is_axes, settings.scales_bp):
        frame = scale_plot_frame(
            region, track_a, track_b, bins_a, bins_b, scale, settings, interval=(start, end)
        )
        for point in frame.loc[frame["unavailable"], "start"]:
            axis.axvspan(point, point + settings.resolution_bp, color="0.85", zorder=0)
        axis.axvspan(region.region_start, region.region_end, color="#dbeafe", zorder=0)
        axis.plot(frame.start, frame.a, color="#2563eb", lw=0.8, label="A")
        axis.plot(frame.start, frame.b, color="#16a34a", lw=0.8, label="B")
        axis.plot(frame.start, frame.b - frame.a, color="#6b21a8", lw=0.8, label="B-A")
        for event in events.itertuples(index=False):
            if event.is_change_type == "broad_shift":
                axis.axvspan(event.event_start, event.event_end, color="#dc2626", alpha=0.20)
            elif pd.notna(event.summit):
                axis.axvline(event.summit, color="#dc2626", lw=1)
        axis.set_ylabel(f"{scale // 1000} kb")
        axis.grid(axis="y", alpha=0.2)
    is_axes[0].legend(loc="upper right", ncol=3, fontsize=8)
    is_axes[-1].set_xlabel(f"{region.chrom} coordinate (bp)")
    for axis in (heatmap_a, heatmap_b, *is_axes[:-1]):
        axis.tick_params(axis="x", labelbottom=False)
    return figure, (heatmap_a, heatmap_b), is_axes


def plot_region(
    region,
    events,
    track_a,
    track_b,
    bins_a,
    bins_b,
    coolers_a,
    coolers_b,
    sample_a,
    sample_b,
    settings,
    output,
):
    output = Path(output)
    figure = None
    temporary = None
    try:
        figure, _, _ = build_region_figure(
            region, events, track_a, track_b, bins_a, bins_b,
            coolers_a, coolers_b, sample_a, sample_b, settings,
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary = tempfile.mkstemp(
            dir=output.parent, prefix=f".{output.name}.", suffix=".png"
        )
        os.close(descriptor)
        figure.savefig(temporary, dpi=140)
        os.replace(temporary, output)
    except BaseException:
        if temporary is not None:
            Path(temporary).unlink(missing_ok=True)
        raise
    finally:
        if figure is not None:
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
    required = {
        "region_id", "chrom", "event_start", "event_end", "summit", "is_change_type",
        "confidence_tier", "rank_within_region", "max_amplitude_percentile",
        "max_transition_percentile", "max_pulse_percentile",
    }
    missing = required.difference(primary.columns)
    if missing:
        raise ValueError(f"primary loci missing columns: {', '.join(sorted(missing))}")
    if not primary["confidence_tier"].eq("primary").all():
        raise ValueError("--primary contains preliminary rows; provide primary-tier loci only")
    required_regions = {"region_id", "chrom", "region_start", "region_end"}
    if required_regions.difference(regions.columns):
        raise ValueError("candidate regions missing required coordinate columns")
    track_a, track_b = read_track(args.is_a, settings), read_track(args.is_b, settings)
    bins_a = read_balanced_bins(args.cool_a, settings.resolution_bp)
    bins_b = read_balanced_bins(args.cool_b, settings.resolution_bp)
    coolers_a = open_coolers_by_chrom(args.cool_a, "sample A")
    coolers_b = open_coolers_by_chrom(args.cool_b, "sample B")
    sample_a = Path(args.is_a).stem.removesuffix("_IS")
    sample_b = Path(args.is_b).stem.removesuffix("_IS")
    top = primary.assign(
        _support_band_count=primary["support_scale_bands"].fillna("").map(
            lambda bands: len({band.strip() for band in str(bands).split(",") if band.strip()})
        )
    ).sort_values(
        [
            "_support_band_count", "max_transition_percentile", "max_pulse_percentile",
            "max_amplitude_percentile", "region_id", "chrom", "event_start", "event_end",
            "is_change_type",
        ],
        ascending=[False, False, False, False, True, True, True, True, True],
        kind="stable",
        na_position="last",
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
        plot_region(
            region, events, track_a, track_b, bins_a, bins_b,
            coolers_a, coolers_b, sample_a, sample_b, settings, output,
        )
        manifest.append({"region_id": region_id, "chrom": region.chrom, "plot_path": str(output.resolve())})
    atomic_manifest(pd.DataFrame(manifest, columns=["region_id", "chrom", "plot_path"]), args.manifest)


if __name__ == "__main__":
    main()
