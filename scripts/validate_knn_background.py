"""Validate an external KNN empirical-p-value background table."""

import argparse
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = (
    "comparison", "chrom", "window_start", "window_end",
    "EC", "IB", "SP", "MSE", "SSIM", "Spearman", "tw_pvalue", "NCC",
)
NUMERIC_COLUMNS = REQUIRED_COLUMNS[2:]


def validate_background_table(path, resolution_bp, window_bins, step_bins, knn_k):
    """Return an audit for a valid scale-specific external background table.

    Raises:
        ValueError: If the background schema, coordinates, or finite-row count is
            incompatible with the requested KNN scale.
    """
    frame = pd.read_csv(path, sep="\t", usecols=lambda column: column in REQUIRED_COLUMNS)
    missing = sorted(set(REQUIRED_COLUMNS).difference(frame.columns))
    if missing:
        raise ValueError(f"background missing columns: {', '.join(missing)}")
    for column in NUMERIC_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    logical_keys = ["comparison", "chrom", "window_start", "window_end"]
    if frame.duplicated(logical_keys).any():
        raise ValueError("background has duplicate comparison/chrom/window coordinates")

    expected_span = resolution_bp * window_bins
    expected_step = resolution_bp * step_bins
    if (
        (frame["window_start"] < 0).any()
        or (frame["window_start"] % resolution_bp != 0).any()
    ):
        raise ValueError(f"background starts must align to {resolution_bp} bp")

    complete_count = 0
    partial_count = 0
    for keys, group in frame.groupby(["comparison", "chrom"], sort=False):
        group = group.sort_values(["window_start", "window_end"], kind="stable")
        differences = group["window_start"].diff().dropna()
        if (differences <= 0).any() or (differences % expected_step != 0).any():
            raise ValueError(f"background start grid mismatch for {keys}")
        spans = group["window_end"] - group["window_start"]
        complete = spans.eq(expected_span)
        partial = spans.lt(expected_span) & spans.ge(expected_span / 2)
        if (~complete & ~partial).any():
            raise ValueError(f"background window span mismatch for {keys}")
        if partial.any() and (partial.sum() != 1 or not partial.iloc[-1]):
            raise ValueError(f"partial background window must be terminal for {keys}")
        complete_count += int(complete.sum())
        partial_count += int(partial.sum())

    finite = np.isfinite(frame[list(NUMERIC_COLUMNS)].to_numpy(dtype=float)).all(axis=1)
    valid_rows = int(finite.sum())
    if valid_rows < knn_k:
        raise ValueError(f"background has {valid_rows} valid rows, fewer than knn_k={knn_k}")
    return {
        "path": str(Path(path).resolve()),
        "resolution_bp": resolution_bp,
        "window_bins": window_bins,
        "window_size_bp": expected_span,
        "step_bins": step_bins,
        "step_size_bp": expected_step,
        "total_rows": len(frame),
        "valid_rows": valid_rows,
        "complete_rows": complete_count,
        "terminal_partial_rows": partial_count,
        "comparisons": int(frame["comparison"].nunique()),
        "chromosomes": int(frame["chrom"].nunique()),
    }


def write_audit_atomically(audit, destination):
    """Write *audit* to *destination* without replacing an existing audit early."""
    destination = Path(destination)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            json.dump(audit, temporary, indent=2, sort_keys=True)
            temporary.write("\n")
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Validate a scale-specific external KNN background table."
    )
    parser.add_argument("--background", required=True)
    parser.add_argument("--resolution-bp", required=True, type=int)
    parser.add_argument("--window-bins", required=True, type=int)
    parser.add_argument("--step-bins", required=True, type=int)
    parser.add_argument("--knn-k", required=True, type=int)
    parser.add_argument("--audit", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        audit = validate_background_table(
            args.background,
            args.resolution_bp,
            args.window_bins,
            args.step_bins,
            args.knn_k,
        )
    except ValueError as error:
        raise SystemExit(str(error))
    write_audit_atomically(audit, args.audit)


if __name__ == "__main__":
    main()
