#!/usr/bin/env python3
"""Call ordered A-to-B TAD-boundary events for one comparison."""

import argparse
from pathlib import Path
import sys

import pandas as pd
import yaml

from cworld_is_config import parse_cworld_is_settings
from is_boundary_config import parse_is_boundary_settings
from is_boundary_errors import ISBoundaryError
from is_boundary_io import (
    BOUNDARY_EVENT_COLUMNS,
    BOUNDARY_STABILITY_COLUMNS,
    atomic_write_tsv,
    read_boundary_views,
    read_jitter_offsets,
    read_unified_boundaries,
)
from is_boundary_matching import call_boundary_events, compare_event_views, estimate_jitter
from tad_config import parse_tad_settings


QC_COLUMNS = ("metric", "value", "detail")


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", required=True)
    parser.add_argument("--sample-a", required=True)
    parser.add_argument("--sample-b", required=True)
    parser.add_argument("--boundaries-a", required=True)
    parser.add_argument("--boundaries-b", required=True)
    parser.add_argument("--offsets-a", required=True)
    parser.add_argument("--offsets-b", required=True)
    parser.add_argument("--views-a", required=True)
    parser.add_argument("--views-b", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-events", required=True)
    parser.add_argument("--output-stability", required=True)
    parser.add_argument("--output-qc", required=True)
    parser.add_argument("--log", required=True)
    return parser


def _settings(path):
    try:
        config = yaml.safe_load(Path(path).read_text())
        return parse_is_boundary_settings(
            config, parse_cworld_is_settings(config), parse_tad_settings(config)
        )
    except Exception as error:
        if isinstance(error, ISBoundaryError):
            raise
        raise ISBoundaryError("IS_BOUNDARY_CONFIG_INVALID", reason=type(error).__name__) from error


def _reference_mismatch(message):
    raise ISBoundaryError("IS_BOUNDARY_REFERENCE_MISMATCH", reason=message)


def _validate_reference(boundaries_a, boundaries_b, settings, sample_a, sample_b):
    if sample_a == sample_b:
        _reference_mismatch("sample_names_equal")
    if set(boundaries_a["chrom"]) != set(boundaries_b["chrom"]):
        _reference_mismatch("chromosome_sets_differ")
    coordinates = pd.concat([
        boundaries_a[["start", "end"]], boundaries_b[["start", "end"]],
    ], ignore_index=True)
    if not coordinates.empty and (coordinates.to_numpy(dtype=int) % settings.base_bin_bp).any():
        _reference_mismatch("coordinates_not_on_base_grid")


def _validate_sample(frame, sample, label):
    if not frame.empty and set(frame["sample"]) != {sample}:
        _reference_mismatch(f"{label}_sample_mismatch")


def _qc(events, stability, jitter, settings):
    rows = [
        {"metric": "jitter_threshold_bp", "value": jitter.threshold_bp, "detail": ""},
        {"metric": "jitter_source", "value": jitter.source, "detail": ""},
        {"metric": "n_offsets", "value": jitter.n_offsets, "detail": ""},
        {"metric": "raw_quantile_bp", "value": jitter.raw_quantile_bp, "detail": ""},
    ]
    for state, count in events["position_state"].value_counts(sort=False).sort_index().items():
        rows.append({"metric": f"position_state_{state}", "value": int(count), "detail": ""})
    for radius in settings.shift_sensitivity_bp:
        column = f"state_at_{radius}"
        rows.append({
            "metric": f"shift_radius_concordance_{radius}",
            "value": int(events[column].eq(events["position_state"]).sum()),
            "detail": f"of_{len(events)}",
        })
    for view in ("raw_only_input", "balanced_only_input"):
        subset = stability[stability["view"] == view]
        rows.append({
            "metric": f"{view}_concordance", "value": int(subset["concordant"].sum()),
            "detail": f"of_{len(subset)}",
        })
    leave_one_out = stability[stability["view"].str.startswith("leave_r", na=False)]
    rows.append({
        "metric": "leave_one_resolution_out_concordance",
        "value": int(leave_one_out["concordant"].sum()),
        "detail": f"of_{len(leave_one_out)}",
    })
    return pd.DataFrame(rows, columns=QC_COLUMNS)


def _write_log(path, message):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("a") as handle:
        handle.write(f"{message}\n")


def run_from_args(args):
    settings = _settings(args.config)
    boundaries_a = read_unified_boundaries(args.boundaries_a)
    boundaries_b = read_unified_boundaries(args.boundaries_b)
    offsets_a = read_jitter_offsets(args.offsets_a)
    offsets_b = read_jitter_offsets(args.offsets_b)
    views_a = read_boundary_views(args.views_a)
    views_b = read_boundary_views(args.views_b)
    _validate_reference(boundaries_a, boundaries_b, settings, args.sample_a, args.sample_b)
    _validate_sample(offsets_a, args.sample_a, "offsets_a")
    _validate_sample(offsets_b, args.sample_b, "offsets_b")
    _validate_sample(views_a, args.sample_a, "views_a")
    _validate_sample(views_b, args.sample_b, "views_b")
    try:
        jitter = estimate_jitter(offsets_a, offsets_b, settings)
    except Exception as error:
        raise ISBoundaryError(
            "IS_BOUNDARY_JITTER_CALIBRATION_FAILED", error=type(error).__name__
        ) from error
    try:
        events = call_boundary_events(
            boundaries_a, boundaries_b, jitter, settings, args.comparison,
            args.sample_a, args.sample_b,
        )
        stability = compare_event_views(
            events, views_a, views_b, jitter, settings, args.comparison,
            args.sample_a, args.sample_b,
        )
    except Exception as error:
        raise ISBoundaryError("IS_BOUNDARY_MATCH_FAILED", error=type(error).__name__) from error
    qc = _qc(events, stability, jitter, settings)
    atomic_write_tsv(events, args.output_events, BOUNDARY_EVENT_COLUMNS)
    atomic_write_tsv(stability, args.output_stability, BOUNDARY_STABILITY_COLUMNS)
    atomic_write_tsv(qc, args.output_qc, QC_COLUMNS)
    _write_log(
        args.log,
        f"comparison={args.comparison} events={len(events)} stability_rows={len(stability)} "
        f"jitter_bp={jitter.threshold_bp} jitter_source={jitter.source}",
    )


def cli(argv=None):
    args = build_parser().parse_args(argv)
    try:
        run_from_args(args)
    except ISBoundaryError as error:
        print(error, file=sys.stderr)
        return 1
    except Exception as error:
        print(
            ISBoundaryError("IS_BOUNDARY_MATCH_FAILED", error=type(error).__name__),
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
