#!/usr/bin/env python3
"""Run optional, constrained set-level IS/TAD-boundary permutations."""

import argparse
from pathlib import Path
import sys

import pandas as pd
import yaml

from call_is_change_loci import read_regions
from cworld_is_config import parse_cworld_is_settings
from is_boundary_config import parse_is_boundary_settings
from is_boundary_errors import ISBoundaryError
from is_boundary_io import BOUNDARY_EVENT_COLUMNS, atomic_write_tsv
from is_boundary_permutation import PERMUTATION_COLUMNS, run_permutations
from tad_config import parse_tad_settings


MASK_COLUMNS = ("chrom", "start", "end", "reason")


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loci", required=True)
    parser.add_argument("--candidate-regions", required=True)
    parser.add_argument("--events", required=True)
    parser.add_argument("--masked", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
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


def _read_masked(path):
    try:
        frame = pd.read_csv(path, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=MASK_COLUMNS)
    if frame.empty:
        return pd.DataFrame(columns=MASK_COLUMNS)
    if tuple(frame.iloc[0].astype(str)) == MASK_COLUMNS:
        frame = frame.iloc[1:].reset_index(drop=True)
    if frame.shape[1] != len(MASK_COLUMNS):
        raise ValueError("masked table must have four columns")
    frame.columns = MASK_COLUMNS
    frame[["start", "end"]] = frame[["start", "end"]].apply(pd.to_numeric, errors="raise")
    return frame


def _read_events(path):
    frame = pd.read_csv(path, sep="\t")
    if tuple(frame.columns) != BOUNDARY_EVENT_COLUMNS:
        raise ValueError("unexpected event schema")
    return frame


def _write_log(path, message):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("a") as handle:
        handle.write(f"{message}\n")


def run_from_args(args):
    settings = _settings(args.config)
    try:
        loci = pd.read_csv(args.loci, sep="\t")
        regions = read_regions(args.candidate_regions, settings.base_bin_bp)
        events = _read_events(args.events)
        masked = _read_masked(args.masked)
        result = run_permutations(loci, regions, events, masked, settings)
        atomic_write_tsv(result, args.output, PERMUTATION_COLUMNS)
    except Exception as error:
        if isinstance(error, ISBoundaryError):
            raise
        raise ISBoundaryError("IS_BOUNDARY_ANNOTATION_FAILED", stage="permutation", error=type(error).__name__) from error
    _write_log(
        args.log,
        f"comparison={result.iloc[0].comparison} tier={result.iloc[0].confidence_tier} "
        f"permutations={settings.permutation.n_permutations}",
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
