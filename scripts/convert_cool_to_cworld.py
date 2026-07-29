#!/usr/bin/env python3
"""
Convert balanced cooler (.cool) to cworld-dekker matrix format.

Uses cooler built-in balance: clr.matrix(balance=True).fetch().
Returns balanced float64 dense matrix directly.

cworld matrix format:
    Header: [empty]\tHiC_1|aaa|chrom:start-end\t...
    Rows: HiC_i|aaa|chrom:start-end\tval\tval\t...

Usage:
    python scripts/convert_cool_to_cworld.py \
        --cool results/balanced/HC04_rep1/HC04_rep1__HC04_A01_iced.cool \
        --chrom HC04_A01 \
        --out /tmp/HC04_rep1__HC04_A01.matrix.gz \
        --log /tmp/convert.log
"""

import argparse
import gzip
import sys
from pathlib import Path

import cooler
import numpy as np


def log(msg, log_file):
    print(msg)
    sys.stdout.flush()
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def convert_cool_to_cworld(
    cool_path, chrom, output_path, log_file=None, expected_bin_size=None
):
    output_path = Path(output_path)
    log(f"[INFO] Opening cooler: {cool_path}", log_file)
    clr = cooler.Cooler(cool_path)

    resolution = int(clr.info["bin-size"])
    if expected_bin_size is not None and resolution != expected_bin_size:
        raise ValueError(
            f"Cooler {cool_path} has wrong bin size: expected {expected_bin_size} bp, "
            f"found {resolution} bp"
        )

    if chrom not in clr.chromnames:
        raise ValueError(f"Chromosome {chrom} not found in {cool_path}")

    bins = clr.bins()[:]
    chrom_bins = bins[bins["chrom"] == chrom].reset_index(drop=True)
    n_bins = len(chrom_bins)
    log(f"[INFO] Chrom {chrom}: {n_bins} bins @ {resolution}bp", log_file)

    log(f"[INFO] Fetching balanced matrix...", log_file)
    mat = clr.matrix(balance=True).fetch(chrom).astype(np.float64)
    log(f"[INFO] Matrix loaded: {n_bins}x{n_bins}, NaN cells: {int(np.isnan(mat).sum())}", log_file)

    mat = np.nan_to_num(mat, nan=0.0)

    bin_labels = [
        f"HiC_{i+1}|aaa|{chrom}:{chrom_bins.iloc[i]['start']}-{chrom_bins.iloc[i]['end']}"
        for i in range(n_bins)
    ]

    log(f"[INFO] Writing matrix: {output_path}", log_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, "wt") as f:
        f.write("\t" + "\t".join(bin_labels) + "\n")
        for i in range(n_bins):
            row = mat[i]
            parts = [bin_labels[i]] + [f"{v:.6g}" if v != 0.0 else "0" for v in row]
            f.write("\t".join(parts) + "\n")
            if (i + 1) % 2000 == 0:
                f.flush()
                log(f"[INFO]  row {i+1}/{n_bins}", log_file)
        f.flush()

    log(f"[INFO] Done: {output_path}", log_file)


def run():
    parser = argparse.ArgumentParser(
        description="Convert balanced cooler to cworld matrix."
    )
    parser.add_argument("--cool", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--log")
    parser.add_argument(
        "--expected-bin-size",
        type=int,
        help="Require this cooler bin size in bp before conversion",
    )
    args = parser.parse_args()

    if args.log:
        Path(args.log).parent.mkdir(parents=True, exist_ok=True)

    log(f"[convert_cool_to_cworld] cool={args.cool}", args.log)
    log(f"[convert_cool_to_cworld] chrom={args.chrom}", args.log)
    log(f"[convert_cool_to_cworld] out={args.out}", args.log)

    try:
        convert_cool_to_cworld(
            args.cool,
            args.chrom,
            args.out,
            log_file=args.log,
            expected_bin_size=args.expected_bin_size,
        )
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    run()
