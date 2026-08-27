#!/usr/bin/env python3
"""Validate one scale-specific KNN background and write a JSON audit."""

import argparse
import csv
import json
from math import isfinite
import os
from pathlib import Path
import tempfile


IDENTITY_COLUMNS = ("comparison", "chrom", "window_start", "window_end")
NUMERIC_COLUMNS = ("EC", "IB", "SP", "MSE", "SSIM", "Spearman", "tw_pvalue", "NCC")
REQUIRED_COLUMNS = (*IDENTITY_COLUMNS, *NUMERIC_COLUMNS)


def positive_integer(value: str) -> int:
    try:
        integer = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a positive integer") from error
    if integer <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return integer


def validate_background(path: Path, scale, knn_k: int) -> dict:
    if type(knn_k) is not int or knn_k <= 0:
        raise ValueError("knn_k must be a positive integer")
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"background is missing or empty: {path}")

    resolution = int(scale["resolution_bp"])
    window_bins = int(scale["window_bins"])
    window_size = int(scale["window_size_bp"])
    step_bins = int(scale["step_bins"])
    step_size = int(scale["step_size_bp"])
    minimum_tail = ((window_bins // 2) - 1) * resolution + 1
    identities = set()
    row_count = 0
    valid_knn_rows = 0
    observed_lengths = set()

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(set(REQUIRED_COLUMNS).difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"background missing columns: {', '.join(missing)}")
        for line_number, record in enumerate(reader, start=2):
            row_count += 1
            try:
                start = int(record["window_start"])
                end = int(record["window_end"])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number} has non-integer coordinates"
                ) from error
            if start % step_size:
                raise ValueError(
                    f"{path}:{line_number} window_start {start} is not on the "
                    f"{step_size} bp step grid"
                )
            length = end - start
            observed_lengths.add(length)
            if not minimum_tail <= length <= window_size:
                raise ValueError(
                    f"{path}:{line_number} window length {length} is outside "
                    f"[{minimum_tail}, {window_size}]"
                )
            identity = (
                record["comparison"],
                record["chrom"],
                start,
                end,
            )
            if identity in identities:
                raise ValueError(
                    f"{path}:{line_number} duplicate logical background row: {identity}"
                )
            identities.add(identity)
            try:
                numeric = [float(record[column]) for column in NUMERIC_COLUMNS]
            except ValueError:
                numeric = []
            if numeric and all(isfinite(value) for value in numeric):
                valid_knn_rows += 1

    if row_count == 0:
        raise ValueError(f"background has no data rows: {path}")
    if valid_knn_rows < knn_k:
        raise ValueError(
            f"valid KNN rows {valid_knn_rows} is smaller than knn_k {knn_k}"
        )
    return {
        "background_file": str(path.resolve()),
        "resolution_bp": resolution,
        "window_bins": window_bins,
        "window_size_bp": window_size,
        "step_bins": step_bins,
        "step_size_bp": step_size,
        "row_count": row_count,
        "valid_knn_rows": valid_knn_rows,
        "observed_window_lengths_bp": sorted(observed_lengths),
    }


def atomic_json(data, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "w") as handle:
            json.dump(data, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, output)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--background", required=True, type=Path)
    parser.add_argument("--resolution-bp", required=True, type=positive_integer)
    parser.add_argument("--window-bins", required=True, type=positive_integer)
    parser.add_argument("--step-bins", required=True, type=positive_integer)
    parser.add_argument("--knn-k", required=True, type=int)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    scale = {
        "resolution_bp": args.resolution_bp,
        "window_bins": args.window_bins,
        "window_size_bp": args.resolution_bp * args.window_bins,
        "step_bins": args.step_bins,
        "step_size_bp": args.resolution_bp * args.step_bins,
    }
    audit = validate_background(args.background, scale, args.knn_k)
    atomic_json(audit, args.out)
    print(
        f"[OK] r{args.resolution_bp} w{args.window_bins}: "
        f"{audit['row_count']} rows"
    )


if __name__ == "__main__":
    main()
