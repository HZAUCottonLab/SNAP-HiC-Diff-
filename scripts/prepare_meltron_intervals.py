#!/usr/bin/env python3
"""
Prepare Meltron genomic intervals from KNN enriched windows.

Steps:
1. Read knn_enriched.tsv
2. Require the five canonical raw empirical p-value columns
3. Filter: windows with >= min_significant metrics raw P < raw_pvalue_threshold
4. Output strict 4-col Meltron interval format
5. Output separate interval_map.tsv for downstream join

Usage:
    python scripts/prepare_meltron_intervals.py \
        --enriched results/knn_empirical_pvalue/Mutant1_vs_Mutant2/Mutant1_vs_Mutant2_knn_enriched.tsv \
        --out results/meltron/Mutant1_vs_Mutant2/Mutant1_vs_Mutant2_intervals.tsv \
        --min-significant 1 \
        --raw-pvalue-threshold 0.05 \
        --log-file results/logs/meltron/Mutant1_vs_Mutant2.log
"""

import argparse
import os
import sys

import pandas as pd


CANONICAL_PVALUE_COLUMNS = (
    "MSE_empirical_pvalue",
    "d_SSIM_empirical_pvalue",
    "d_Spearman_empirical_pvalue",
    "tw_score_empirical_pvalue",
    "d_NCC_empirical_pvalue",
)

SCALE_COLUMNS = (
    "resolution_bp",
    "window_bins",
    "window_size_bp",
    "step_bins",
    "step_size_bp",
)

INTERVAL_MAP_COLUMNS = [
    *SCALE_COLUMNS,
    "chrom",
    "interval_start",
    "interval_end",
    "window_start",
    "window_end",
    "window_id",
    "ID",
    "n_significant_metrics",
]

def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")




def run():
    parser = argparse.ArgumentParser(
        description="Prepare Meltron genomic intervals from KNN enriched windows."
    )
    parser.add_argument(
        "--enriched", required=True,
        help="Path to knn_enriched.tsv"
    )
    parser.add_argument(
        "--out", required=True,
        help="Path to output Meltron intervals TSV (4 cols: chrom, start, end, ID)"
    )
    parser.add_argument(
        "--min-significant", type=int, default=1,
        help="Minimum number of raw empirical p-values below the threshold (default: 1)"
    )
    parser.add_argument(
        "--raw-pvalue-threshold", type=float, default=0.05,
        help="Raw empirical p-value threshold for significance (default: 0.05)"
    )
    parser.add_argument(
        "--log-file",
        help="Log file path"
    )
    args = parser.parse_args()

    log(f"[prepare_meltron_intervals] Input: {args.enriched}", args.log_file)

    if not os.path.exists(args.enriched):
        log(f"[ERROR] Input file not found: {args.enriched}", args.log_file)
        sys.exit(1)

    df = pd.read_csv(args.enriched, sep="\t")
    log(f"[INFO] knn_enriched rows: {len(df)}", args.log_file)

    canonical_coordinates = {"start", "end"}
    window_coordinates = {"window_start", "window_end"}
    has_canonical_coordinates = canonical_coordinates.issubset(df.columns)
    has_window_coordinates = window_coordinates.issubset(df.columns)
    if not (has_canonical_coordinates or has_window_coordinates):
        parser.error(
            "input coordinate schema requires either start/end or window_start/window_end"
        )
    if has_canonical_coordinates:
        if not has_window_coordinates:
            df["window_start"] = df["start"]
            df["window_end"] = df["end"]
    else:
        df["start"] = df["window_start"]
        df["end"] = df["window_end"]

    # Identify p-value columns (empirical p-values from KNN comparison)
    missing_scale_columns = [
        column for column in SCALE_COLUMNS if column not in df.columns
    ]
    if missing_scale_columns:
        parser.error(
            "input requires scale provenance columns: "
            + ", ".join(missing_scale_columns)
        )
    empirical_pvalue_columns = [
        column for column in df.columns if column.endswith("_empirical_pvalue")
    ]
    if set(empirical_pvalue_columns) != set(CANONICAL_PVALUE_COLUMNS):
        parser.error("input requires exactly the canonical empirical p-value columns")
    df[list(SCALE_COLUMNS)] = df[list(SCALE_COLUMNS)].apply(
        pd.to_numeric, errors="raise"
    )
    raw_pvalues = df[list(CANONICAL_PVALUE_COLUMNS)].apply(
        pd.to_numeric, errors="coerce"
    )
    df["n_significant_metrics"] = (
        raw_pvalues < args.raw_pvalue_threshold
    ).sum(axis=1)

    log(f"[INFO] Raw empirical p-value columns: {list(CANONICAL_PVALUE_COLUMNS)}", args.log_file)

    # Filter by min_significant
    n_before = len(df)
    df_filtered = df[df["n_significant_metrics"] >= args.min_significant].copy()
    n_after = len(df_filtered)
    log(f"[INFO] Filtered: {n_before} -> {n_after} (n_significant_metrics >= {args.min_significant})", args.log_file)

    if n_after == 0:
        log("[WARN] No windows pass filter — outputting empty intervals file", args.log_file)
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        interval_cols = ["chrom", "start", "end", "ID"]
        pd.DataFrame(columns=interval_cols).to_csv(args.out, sep="\t", index=False)
        map_out = args.out.replace("_intervals.tsv", "_interval_map.tsv")
        pd.DataFrame(columns=INTERVAL_MAP_COLUMNS).to_csv(
            map_out, sep="\t", index=False
        )
        sys.exit(0)

    # Assign ID = {chrom}_{window_id}
    df_filtered["ID"] = (
        "r" + df_filtered["resolution_bp"].astype(int).astype(str)
        + "_w" + df_filtered["window_bins"].astype(int).astype(str)
        + "_" + df_filtered["chrom"].astype(str)
        + "_" + df_filtered["window_id"].astype(str)
    )

    # Prepare interval_map: for joining meltron scores back to original windows
    interval_map = df_filtered[
        [
            *SCALE_COLUMNS,
            "chrom", "start", "end", "window_start", "window_end",
            "window_id", "ID", "n_significant_metrics",
        ]
    ].rename(
        columns={"start": "interval_start", "end": "interval_end"}
    )
    interval_map = interval_map.loc[:, INTERVAL_MAP_COLUMNS]

    # Prepare strict 4-col Meltron interval output
    intervals_out = df_filtered[["chrom", "start", "end", "ID"]].copy()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    intervals_out.to_csv(args.out, sep="\t", index=False)
    log(f"[INFO] Intervals written: {args.out} ({len(intervals_out)} rows)", args.log_file)

    # Write interval_map alongside
    map_out = args.out.replace("_intervals.tsv", "_interval_map.tsv")
    interval_map.to_csv(map_out, sep="\t", index=False)
    log(f"[INFO] Interval map written: {map_out}", args.log_file)


if __name__ == "__main__":
    run()
