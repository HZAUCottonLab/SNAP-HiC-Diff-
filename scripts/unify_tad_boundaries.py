#!/usr/bin/env python3
"""Build one provenance-preserving raw/balanced boundary view per sample."""

import argparse
from pathlib import Path
import sys

import cooler
import numpy as np
import pandas as pd
import yaml

from cworld_is_config import parse_cworld_is_settings
from is_boundary_config import parse_is_boundary_settings
from is_boundary_errors import ISBoundaryError
from is_boundary_io import (
    BOUNDARY_VIEW_COLUMNS,
    TAD_UNEVALUABLE_COLUMNS,
    UNIFIED_BOUNDARY_COLUMNS,
    atomic_write_tsv,
    read_multires_boundaries,
    read_resolution_boundaries,
)
from is_boundary_unification import (
    JITTER_OFFSET_COLUMNS,
    build_jitter_offsets,
    build_sensitivity_views,
    unify_boundaries,
)
from tad_config import parse_tad_settings


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--raw-multires", required=True)
    parser.add_argument("--balanced-multires", required=True)
    parser.add_argument("--raw-resolution-boundaries", required=True, nargs=3)
    parser.add_argument("--balanced-resolution-boundaries", required=True, nargs=3)
    parser.add_argument("--mcool", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-unified", required=True)
    parser.add_argument("--output-offsets", required=True)
    parser.add_argument("--output-views", required=True)
    parser.add_argument("--output-tad-unevaluable", required=True)
    parser.add_argument("--log", required=True)
    return parser


def _settings(path):
    try:
        path = Path(path).resolve()
        config = yaml.safe_load(path.read_text())
        return config, parse_is_boundary_settings(
            config, parse_cworld_is_settings(config), parse_tad_settings(config)
        ), path.parent
    except ValueError as error:
        raise ISBoundaryError("IS_BOUNDARY_CONFIG_INVALID", reason=str(error)) from error


def _read_bad_regions(path):
    if not path:
        return []
    rows = []
    try:
        for line_number, line in enumerate(Path(path).read_text().splitlines(), 1):
            if not line or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.split()
            start, end = int(start), int(end)
            if start < 0 or end <= start:
                raise ValueError("coordinates")
            rows.append((chrom, start, end))
    except (OSError, ValueError) as error:
        raise ISBoundaryError(
            "IS_BOUNDARY_INPUT_SCHEMA_INVALID", path=path, reason=f"bad_region_{type(error).__name__}",
            line_number=locals().get("line_number"),
        ) from error
    return rows


def _mcool_metadata(path, resolution):
    try:
        clr = cooler.Cooler(f"{Path(path).resolve()}::/resolutions/{resolution}")
        bins = clr.bins()[:][["chrom", "start", "end", "weight"]]
    except Exception as error:
        raise ISBoundaryError(
            "IS_BOUNDARY_INPUT_SCHEMA_INVALID", path=path, reason=f"mcool_{type(error).__name__}"
        ) from error
    chrom_sizes = {str(chrom): int(length) for chrom, length in clr.chromsizes.items()}
    return list(chrom_sizes), chrom_sizes, bins


def _attach_eligibility(frame, chrom_sizes, bins, bad_regions, base_bin_bp):
    frame = frame.copy()
    if frame.empty:
        return frame.assign(
            touches_chromosome_end=pd.Series(dtype=bool),
            has_nonfinite_base_bins=pd.Series(dtype=bool),
            overlaps_bad_region=pd.Series(dtype=bool),
            terminal_partial_coordinate=pd.Series(dtype=bool),
        )
    columns = {chrom: group for chrom, group in bins.groupby("chrom", observed=True)}
    end_flags, missing_flags, bad_flags = [], [], []
    for row in frame.itertuples(index=False):
        if (
            row.chrom not in chrom_sizes
            or row.start < 0
            or row.end <= row.start
            or row.end > chrom_sizes[row.chrom]
            or row.start % base_bin_bp
            or (row.end % base_bin_bp and row.end != chrom_sizes[row.chrom])
        ):
            raise ISBoundaryError(
                "IS_BOUNDARY_INPUT_SCHEMA_INVALID", chrom=row.chrom,
                start=row.start, end=row.end, reason="mcool_coordinate",
            )
        chrom_bins = columns.get(row.chrom, pd.DataFrame())
        covered = chrom_bins[(chrom_bins["start"] < row.end) & (chrom_bins["end"] > row.start)]
        end_flags.append(row.start == 0 or row.end == chrom_sizes[row.chrom])
        missing_flags.append(covered.empty or not np.isfinite(covered["weight"]).all())
        bad_flags.append(any(
            chrom == row.chrom and max(start, row.start) < min(end, row.end)
            for chrom, start, end in bad_regions
        ))
    return frame.assign(
        touches_chromosome_end=end_flags,
        has_nonfinite_base_bins=missing_flags,
        overlaps_bad_region=bad_flags,
        terminal_partial_coordinate=[
            row.end == chrom_sizes[row.chrom] and row.end % base_bin_bp != 0
            for row in frame.itertuples(index=False)
        ],
    )


def _tad_unevaluable(sample, chrom_sizes, bins, bad_regions, base_bin_bp):
    """Publish intervals where absence of a boundary cannot be interpreted."""
    rows = []
    for chrom, chrom_size in chrom_sizes.items():
        first_end = min(base_bin_bp, chrom_size)
        last_start = max(0, chrom_size - base_bin_bp)
        rows.append((sample, chrom, 0, first_end, "terminal_insufficient_tad_coverage"))
        rows.append((sample, chrom, last_start, chrom_size, "terminal_insufficient_tad_coverage"))
    for row in bins.loc[~np.isfinite(bins["weight"])].itertuples(index=False):
        rows.append((sample, row.chrom, int(row.start), int(row.end), "nonfinite_base_bins"))
    for chrom, start, end in bad_regions:
        if chrom in chrom_sizes:
            rows.append((sample, chrom, max(0, start), min(chrom_sizes[chrom], end), "bad_region"))
    mask = pd.DataFrame(rows, columns=TAD_UNEVALUABLE_COLUMNS)
    return mask.loc[mask["end"] > mask["start"]].sort_values(
        ["chrom", "start", "end", "reason"], kind="stable"
    ).reset_index(drop=True)


def _sort_for_chromosomes(frame, chrom_order):
    if frame.empty:
        return frame
    ranks = {chrom: index for index, chrom in enumerate(chrom_order)}
    if set(frame["chrom"]).difference(ranks):
        raise ISBoundaryError("IS_BOUNDARY_INPUT_SCHEMA_INVALID", reason="unknown_chromosome")
    return frame.assign(_chrom_rank=frame["chrom"].map(ranks)).sort_values(
        ["_chrom_rank", "start", "end", "boundary_id"], kind="stable"
    ).drop(columns="_chrom_rank").reset_index(drop=True)


def _write_log(path, message):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("a") as handle:
        handle.write(f"{message}\n")


def run_from_args(args):
    config, settings, config_dir = _settings(args.config)
    chrom_order, chrom_sizes, bins = _mcool_metadata(args.mcool, settings.base_bin_bp)
    bad_region_bed = config["tad"].get("bad_region_bed")
    if bad_region_bed and not Path(bad_region_bed).is_absolute():
        bad_region_bed = config_dir / bad_region_bed
    bad_regions = _read_bad_regions(
        bad_region_bed if config["tad"].get("badregionfilter") else None
    )
    raw = read_multires_boundaries(args.raw_multires, "raw")
    balanced = read_multires_boundaries(args.balanced_multires, "balanced")
    resolutions = tuple(sorted(config["tad"]["resolutions"]))
    raw_resolution = read_resolution_boundaries(
        zip(resolutions, args.raw_resolution_boundaries), "raw"
    )
    balanced_resolution = read_resolution_boundaries(
        zip(resolutions, args.balanced_resolution_boundaries), "balanced"
    )
    raw = _sort_for_chromosomes(raw, chrom_order)
    balanced = _sort_for_chromosomes(balanced, chrom_order)
    per_resolution = _sort_for_chromosomes(
        pd.concat([raw_resolution, balanced_resolution], ignore_index=True), chrom_order
    )
    raw = _attach_eligibility(raw, chrom_sizes, bins, bad_regions, settings.base_bin_bp)
    balanced = _attach_eligibility(
        balanced, chrom_sizes, bins, bad_regions, settings.base_bin_bp
    )
    per_resolution = _attach_eligibility(
        per_resolution, chrom_sizes, bins, bad_regions, settings.base_bin_bp
    )
    raw = raw.assign(sample=args.sample)
    balanced = balanced.assign(sample=args.sample)
    unified = unify_boundaries(
        raw, balanced, settings.match_tie_tolerance_bp, chrom_order=chrom_order
    )
    offsets = build_jitter_offsets(
        raw, balanced, per_resolution, settings.base_bin_bp,
        settings.match_tie_tolerance_bp,
    )
    views = build_sensitivity_views(
        args.sample, raw, balanced, per_resolution, settings
    )
    tad_unevaluable = _tad_unevaluable(
        args.sample, chrom_sizes, bins, bad_regions, settings.base_bin_bp
    )
    if tuple(unified.columns) != UNIFIED_BOUNDARY_COLUMNS:
        raise ISBoundaryError("IS_BOUNDARY_INPUT_SCHEMA_INVALID", reason="unified_schema")
    if tuple(offsets.columns) != JITTER_OFFSET_COLUMNS:
        raise ISBoundaryError("IS_BOUNDARY_INPUT_SCHEMA_INVALID", reason="offset_schema")
    if tuple(views.columns) != BOUNDARY_VIEW_COLUMNS:
        raise ISBoundaryError("IS_BOUNDARY_INPUT_SCHEMA_INVALID", reason="view_schema")
    if tuple(tad_unevaluable.columns) != TAD_UNEVALUABLE_COLUMNS:
        raise ISBoundaryError("IS_BOUNDARY_INPUT_SCHEMA_INVALID", reason="tad_unevaluable_schema")
    atomic_write_tsv(unified, args.output_unified, UNIFIED_BOUNDARY_COLUMNS)
    atomic_write_tsv(offsets, args.output_offsets, JITTER_OFFSET_COLUMNS)
    atomic_write_tsv(views, args.output_views, BOUNDARY_VIEW_COLUMNS)
    atomic_write_tsv(tad_unevaluable, args.output_tad_unevaluable, TAD_UNEVALUABLE_COLUMNS)
    _write_log(
        args.log,
        f"sample={args.sample} unified={len(unified)} offsets={len(offsets)} views={len(views)}",
    )


def cli(argv=None):
    args = build_parser().parse_args(argv)
    try:
        run_from_args(args)
    except ISBoundaryError as error:
        print(error, file=sys.stderr)
        return 1
    except Exception as error:
        wrapped = ISBoundaryError(
            "IS_BOUNDARY_ANNOTATION_FAILED", sample=args.sample,
            stage="unification", error=type(error).__name__,
        )
        print(wrapped, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
