#!/usr/bin/env python3
"""
Convert cooltools insulation TSV output to Meltron IS input format.

cooltools output columns:
    chrom, start, end,
    log2_insulation_score_100000, n_valid_pixels_100000, boundary_strength_100000, is_boundary_100000,
    log2_insulation_score_200000, ...

Meltron required format:
    chrom  viewpoint  IS_distance  IS
    Ghir_A01  4000    100000       0.123
    Ghir_A01  4000    200000       0.234

One TSV per sample (all chromosomes merged).

Usage:
    python scripts/convert_cooltools_IS_to_meltron_IS.py \
        --sample HC04_rep1 \
        --chroms "Ghir_A01 Ghir_A02 ..." \
        --outdir results \
        --is-distances 100000,200000,500000,1000000 \
        --log-file results/logs/compute_IS/HC04_rep1.log
"""

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd


def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def natural_sort_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]


def convert_chrom(insulation_tsv, distances):
    """
    Convert a single cooltools insulation TSV to Meltron long format.

    Parameters
    ----------
    insulation_tsv : str
        Path to cooltools output TSV (one chromosome).
    distances : list of int
        IS distance values.

    Returns
    -------
    pd.DataFrame
        Meltron format: chrom, viewpoint, IS_distance, IS
    """
    df = pd.read_csv(insulation_tsv, sep='\t')

    rows = []
    for dist in distances:
        # cooltools column name: log2_insulation_score_{DISTANCE}
        col = f"log2_insulation_score_{dist}"
        if col not in df.columns:
            log(f"[WARN] Column '{col}' not found in {insulation_tsv} — skipping distance {dist}", None)
            continue

        # Viewpoint = start of the bin (midpoint between start and end also acceptable,
        # but Meltron uses start/viewpoint; use midpoint for better precision)
        subset = df[["chrom", "start", "end", col]].copy()
        subset.rename(columns={col: "IS"}, inplace=True)
        subset["IS_distance"] = int(dist)
        subset["viewpoint"] = (subset["start"] + subset["end"]) // 2
        rows.append(subset[["chrom", "viewpoint", "IS_distance", "IS"]])

    return pd.concat(rows, ignore_index=True)


def run():
    parser = argparse.ArgumentParser(
        description="Convert cooltools insulation TSVs to Meltron IS format."
    )
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument(
        "--chroms", required=True,
        help="Space-separated chromosome list"
    )
    parser.add_argument("--outdir", required=True, help="Top-level results directory (OUTDIR)")
    parser.add_argument(
        "--is-distances", required=True,
        help="Comma-separated IS distances in bp"
    )
    parser.add_argument(
        "--log-file",
        help="Log file path"
    )
    args = parser.parse_args()

    distances = [int(d.strip()) for d in args.is_distances.split(",")]
    chrom_list = args.chroms.split()

    insulation_dir = os.path.join(args.outdir, "insulation", args.sample)
    is_out = os.path.join(args.outdir, "IS", args.sample, f"{args.sample}_IS.tsv")

    all_rows = []
    for chrom in chrom_list:
        base = f"{args.sample}__{chrom}"
        tsv = os.path.join(insulation_dir, f"{base}_insulation.tsv")

        if not os.path.exists(tsv):
            log(f"[WARN] Insulation TSV not found: {tsv} — skipping", args.log_file)
            continue

        log(f"[INFO] Converting {tsv}", args.log_file)
        chrom_df = convert_chrom(tsv, distances)
        all_rows.append(chrom_df)

    if not all_rows:
        log("[ERROR] No insulation TSVs found — cannot produce Meltron IS file", args.log_file)
        sys.exit(1)

    merged = pd.concat(all_rows, ignore_index=True)
    merged = merged.sort_values(
        by="chrom",
        key=lambda col: col.map(natural_sort_key)
    ).reset_index(drop=True)

    merged = merged[["chrom", "viewpoint", "IS_distance", "IS"]]

    os.makedirs(os.path.dirname(is_out), exist_ok=True)
    merged.to_csv(is_out, sep="\t", index=False)

    log(f"[INFO] Meltron IS file written: {is_out} ({len(merged)} rows)", args.log_file)


if __name__ == "__main__":
    run()
