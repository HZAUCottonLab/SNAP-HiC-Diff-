#!/usr/bin/env python3
"""Annotate primary and preliminary IS-change loci with TAD-boundary events."""

import argparse
from pathlib import Path
import sys

import pandas as pd
import yaml

from call_is_change_loci import PRIMARY_COLUMNS, index_is_track, read_is_track
from cworld_is_config import parse_cworld_is_settings
from is_boundary_annotation import annotation_qc, annotate_loci, summarize_annotations
from is_boundary_config import parse_is_boundary_settings
from is_boundary_errors import ISBoundaryError
from is_boundary_io import (
    ANNOTATION_COLUMNS,
    BOUNDARY_EVENT_COLUMNS,
    UNIFIED_BOUNDARY_COLUMNS,
    STRENGTH_LONG_COLUMNS,
    STRENGTH_SUMMARY_COLUMNS,
    TAD_UNEVALUABLE_COLUMNS,
    atomic_write_tables,
    read_tad_unevaluable,
    read_unified_boundaries,
)
from is_boundary_point_state import (
    POINT_BOUNDARY_COLUMNS, POINT_STATE_COLUMNS, BoundaryIndex, attach_point_state,
)
from is_boundary_strength import compute_strength_long, summarize_strength
from is_change_config import parse_is_change_settings
from tad_config import parse_tad_settings


QC_COLUMNS = ("metric", "value", "detail")
SUMMARY_COLUMNS = (
    "confidence_tier", "position_state", "strength_state", "is_consistency_state",
    "evidence_tier", "n_loci", "denominator", "fraction",
)
MASK_COLUMNS = ("chrom", "start", "end", "reason")
SCALED_MASK_COLUMNS = ("chrom", "start", "end", "reason", "scale_bp")


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", required=True)
    parser.add_argument("--events", required=True)
    parser.add_argument("--primary", required=True)
    parser.add_argument("--preliminary", required=True)
    parser.add_argument("--masked", required=True)
    parser.add_argument("--is-a", required=True)
    parser.add_argument("--is-b", required=True)
    parser.add_argument("--tad-unevaluable-a")
    parser.add_argument("--tad-unevaluable-b")
    parser.add_argument("--boundaries-a", required=True)
    parser.add_argument("--boundaries-b", required=True)
    parser.add_argument(
        "--matching-qc",
        help="matching QC table; its jitter_threshold_bp overrides point.shift_tolerance_bp",
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--out-primary", required=True)
    parser.add_argument("--out-preliminary", required=True)
    parser.add_argument("--out-strength", required=True)
    parser.add_argument("--out-strength-summary", required=True)
    parser.add_argument("--out-qc", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-point-map", required=True)
    parser.add_argument("--log", required=True)
    return parser


def _settings(path):
    try:
        config = yaml.safe_load(Path(path).read_text())
        cworld = parse_cworld_is_settings(config)
        return (
            parse_is_boundary_settings(config, cworld, parse_tad_settings(config)),
            parse_is_change_settings(config, cworld),
        )
    except Exception as error:
        if isinstance(error, ISBoundaryError):
            raise
        raise ISBoundaryError("IS_BOUNDARY_CONFIG_INVALID", reason=type(error).__name__) from error


def _read_table(path, columns):
    frame = pd.read_csv(path, sep="\t")
    if tuple(frame.columns) != tuple(columns):
        raise ValueError(f"unexpected schema for {path}")
    return frame


def _read_loci(path):
    frame = pd.read_csv(path, sep="\t")
    missing = set(PRIMARY_COLUMNS).difference(frame.columns)
    if missing:
        raise ValueError(f"IS loci missing columns: {', '.join(sorted(missing))}")
    if frame.columns.duplicated().any():
        raise ValueError("IS loci contain duplicate columns")
    return frame


def _read_masked(path):
    """Read the IS-change mask table.

    Accepts the legacy four-column layout (no scale) and the current
    five-column layout that carries ``scale_bp``.  A legacy row means the
    position is masked at every scale, which preserves the old behaviour.
    """
    try:
        frame = pd.read_csv(path, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=SCALED_MASK_COLUMNS)
    if frame.empty:
        return pd.DataFrame(columns=SCALED_MASK_COLUMNS)
    header = tuple(frame.iloc[0].astype(str))
    if header in (MASK_COLUMNS, SCALED_MASK_COLUMNS):
        frame = frame.iloc[1:].reset_index(drop=True)
    if frame.shape[1] == len(MASK_COLUMNS):
        frame.columns = list(MASK_COLUMNS)
        frame["scale_bp"] = pd.NA
    elif frame.shape[1] == len(SCALED_MASK_COLUMNS):
        frame.columns = list(SCALED_MASK_COLUMNS)
    else:
        raise ValueError("masked IS table must have four or five columns")
    for column in ("start", "end"):
        values = pd.to_numeric(frame[column], errors="raise")
        if not values.eq(values.astype(int)).all():
            raise ValueError("masked IS coordinates must be integers")
        frame[column] = values.astype(int)
    scales = pd.to_numeric(frame["scale_bp"], errors="coerce")
    if scales.notna().any():
        if not scales.dropna().eq(scales.dropna().astype(int)).all():
            raise ValueError("masked IS scale values must be integers")
        frame["scale_bp"] = scales.astype("Int64")
    if (frame["start"] < 0).any() or (frame["end"] <= frame["start"]).any():
        raise ValueError("masked IS coordinates must be non-empty")
    return frame.loc[:, list(SCALED_MASK_COLUMNS)]


def _validate_events(events, comparison):
    if events["boundary_event_id"].duplicated().any():
        raise ValueError("boundary event IDs must be unique")
    comparisons = set(events["comparison"].dropna())
    if comparisons and comparisons != {comparison}:
        raise ValueError("boundary event comparison does not match --comparison")


def _canonicalize_current_locus_comparison(loci, events, comparison):
    """Map only this event map's legacy IS comparison name to its canonical key."""
    pairs = events[["sample_a", "sample_b"]].dropna().drop_duplicates()
    if len(pairs) == 1:
        sample_a, sample_b = pairs.iloc[0]
    elif pairs.empty and "_vs_" in comparison:
        sample_a, sample_b = comparison.split("_vs_", 1)
    else:
        return loci
    canonical = f"{sample_a}_vs_{sample_b}"
    if comparison != canonical:
        return loci
    legacy = f"{sample_a}_IS__vs__{sample_b}_IS"
    return loci.assign(comparison=loci["comparison"].where(loci["comparison"] != legacy, canonical))


def _pair_samples(events, comparison):
    """Sample names for the comparison, taken from the event map."""
    pairs = events[["sample_a", "sample_b"]].dropna().drop_duplicates()
    if len(pairs) == 1:
        return str(pairs.iloc[0]["sample_a"]), str(pairs.iloc[0]["sample_b"])
    if "_vs_" in comparison:
        sample_a, sample_b = comparison.split("_vs_", 1)
        return sample_a, sample_b
    raise ValueError("cannot resolve sample pair for the comparison")


def _jitter_threshold(path, fallback):
    """Empirical jitter threshold from the matching QC table.

    The matching stage already calibrates this value from the observed
    raw-vs-balanced offsets; reusing it keeps the shift rule tied to the data
    rather than to a second guess. A missing or unusable value falls back to
    the configured tolerance.
    """
    try:
        frame = pd.read_csv(path, sep="\t")
    except Exception:
        return fallback
    if not {"metric", "value"}.issubset(frame.columns):
        return fallback
    value = frame.loc[frame["metric"].eq("jitter_threshold_bp"), "value"]
    if value.empty:
        return fallback
    try:
        return int(float(value.iloc[0]))
    except (TypeError, ValueError):
        return fallback


def _write_log(path, message):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("a") as handle:
        handle.write(f"{message}\n")


def _raise(code, error):
    if isinstance(error, ISBoundaryError) and error.code == code:
        raise error
    raise ISBoundaryError(code, error=type(error).__name__) from error


def run_from_args(args):
    boundary_settings, is_settings = _settings(args.config)
    try:
        events = _read_table(args.events, BOUNDARY_EVENT_COLUMNS)
        sample_a, sample_b = _pair_samples(events, args.comparison)
        _validate_events(events, args.comparison)
        primary = _read_loci(args.primary)
        preliminary = _read_loci(args.preliminary)
        primary = _canonicalize_current_locus_comparison(primary, events, args.comparison)
        preliminary = _canonicalize_current_locus_comparison(preliminary, events, args.comparison)
        masked = _read_masked(args.masked)
        tad_unevaluable = pd.concat([
            read_tad_unevaluable(path) if path else pd.DataFrame(columns=TAD_UNEVALUABLE_COLUMNS)
            for path in (args.tad_unevaluable_a, args.tad_unevaluable_b)
        ], ignore_index=True)
    except Exception as error:
        _raise("IS_BOUNDARY_ANNOTATION_FAILED", error)
    try:
        track_a = read_is_track(args.is_a, is_settings)
        track_b = read_is_track(args.is_b, is_settings)
        strength_long = compute_strength_long(
            events, index_is_track(track_a), index_is_track(track_b), boundary_settings, is_settings
        )
        strength_summary = summarize_strength(
            strength_long, events, boundary_settings, is_settings
        )
        if tuple(strength_long.columns) != STRENGTH_LONG_COLUMNS:
            raise ValueError("invalid strength output schema")
        if tuple(strength_summary.columns) != STRENGTH_SUMMARY_COLUMNS:
            raise ValueError("invalid strength summary schema")
    except Exception as error:
        _raise("IS_BOUNDARY_STRENGTH_FAILED", error)
    try:
        annotated_primary = annotate_loci(
            primary, events, strength_summary, masked, boundary_settings, tad_unevaluable,
            is_change_settings=is_settings,
        )
        annotated_preliminary = annotate_loci(
            preliminary, events, strength_summary, masked, boundary_settings, tad_unevaluable,
            is_change_settings=is_settings,
        )
        boundaries_a = BoundaryIndex(read_unified_boundaries(args.boundaries_a))
        boundaries_b = BoundaryIndex(read_unified_boundaries(args.boundaries_b))
        tolerance = boundary_settings.point.shift_tolerance_bp
        if args.matching_qc:
            tolerance = _jitter_threshold(args.matching_qc, tolerance)
        annotated_primary, point_primary = attach_point_state(
            annotated_primary, boundaries_a, boundaries_b, sample_a, sample_b,
            args.comparison, boundary_settings.base_bin_bp, tolerance,
        )
        annotated_preliminary, point_preliminary = attach_point_state(
            annotated_preliminary, boundaries_a, boundaries_b, sample_a, sample_b,
            args.comparison, boundary_settings.base_bin_bp, tolerance,
        )
        point_map = pd.concat([point_primary, point_preliminary], ignore_index=True)
        summary = summarize_annotations(annotated_primary, annotated_preliminary)
        qc = annotation_qc(
            annotated_primary, annotated_preliminary, events, strength_long, strength_summary
        )
        for frame, source in (
            (annotated_primary, primary), (annotated_preliminary, preliminary)
        ):
            # ``locus_key`` is an internal join key added by the point layer; it
            # is not part of the published schema.
            published = [c for c in frame.columns if c != "locus_key"]
            expected = [*source.columns, *ANNOTATION_COLUMNS, *POINT_STATE_COLUMNS]
            if published != expected:
                raise ValueError(
                    "invalid annotation schema: "
                    f"missing={sorted(set(expected) - set(published))} "
                    f"extra={sorted(set(published) - set(expected))}"
                )
        if tuple(point_map.columns) != POINT_BOUNDARY_COLUMNS:
            raise ValueError("invalid point map schema")
        if tuple(summary.columns) != SUMMARY_COLUMNS or tuple(qc.columns) != QC_COLUMNS:
            raise ValueError("invalid annotation summary schema")
        atomic_write_tables([
            (annotated_primary, args.out_primary, annotated_primary.columns),
            (annotated_preliminary, args.out_preliminary, annotated_preliminary.columns),
            (strength_long, args.out_strength, STRENGTH_LONG_COLUMNS),
            (strength_summary, args.out_strength_summary, STRENGTH_SUMMARY_COLUMNS),
            (qc, args.out_qc, QC_COLUMNS),
            (summary, args.out_summary, SUMMARY_COLUMNS),
            (point_map, args.out_point_map, POINT_BOUNDARY_COLUMNS),
        ])
    except Exception as error:
        _raise("IS_BOUNDARY_ANNOTATION_FAILED", error)
    _write_log(
        args.log,
        f"comparison={args.comparison} primary={len(annotated_primary)} "
        f"preliminary={len(annotated_preliminary)} events={len(events)} "
        f"point_rows={len(point_map)}",
    )


def cli(argv=None):
    args = build_parser().parse_args(argv)
    try:
        run_from_args(args)
    except ISBoundaryError as error:
        print(error, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
