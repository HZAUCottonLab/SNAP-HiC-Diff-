#!/usr/bin/env python3
"""
Prepare Meltron genomic intervals from KNN enriched windows.

Steps:
1. Read knn_enriched.tsv
2. Compute BH-FDR for 5 empirical p-value columns
3. Filter: windows with >= min_significant metrics FDR < fdr_threshold
4. Output strict 4-col Meltron interval format
5. Output separate interval_map.tsv for downstream join

Usage:
    python scripts/prepare_meltron_intervals.py \
        --enriched results/knn_empirical_pvalue/Mutant1_vs_Mutant2/Mutant1_vs_Mutant2_knn_enriched.tsv \
        --out results/meltron/Mutant1_vs_Mutant2/Mutant1_vs_Mutant2_intervals.tsv \
        --min-significant 1 \
        --fdr-threshold 0.05 \
        --log-file results/logs/meltron/Mutant1_vs_Mutant2.log
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def benjamini_hochberg(p_values):
    """
    Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    p_values : array-like
        Unadjusted p-values.

    Returns
    -------
    q_values : np.ndarray
        BH-adjusted q-values (FDR).
    """
    p = np.asarray(p_values)
    n = len(p)
    if n == 0:
        return np.array([])

    order = np.argsort(p, kind="quicksort")
    sorted_p = p[order]
    # BH formula: q_i = p_i * n / rank_i
    q = sorted_p * n / np.arange(1, n + 1)
    # Cumulative minimum from the largest q-value backward
    q_cummin = np.minimum.accumulate(q[::-1])[::-1]
    # Map back to original order
    q_values = np.empty(n)
    q_values[order] = q_cummin
    # Cap at 1
    q_values = np.minimum(q_values, 1.0)
    return q_values


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
        help="Minimum number of metrics with FDR < fdr-threshold (default: 1)"
    )
    parser.add_argument(
        "--fdr-threshold", type=float, default=0.05,
        help="FDR threshold for significance (default: 0.05)"
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

    # Identify p-value columns (empirical p-values from KNN comparison)
    pval_cols = [c for c in df.columns if "_empirical_pvalue" in c]
    if not pval_cols:
        log("[WARN] No p-value columns found — treating all windows as candidates", args.log_file)
        pval_cols = []
        df["n_significant"] = 0
    else:
        log(f"[INFO] P-value columns: {pval_cols}", args.log_file)

        # Compute BH-FDR for each p-value column
        for col in pval_cols:
            q_col = col.replace("_pvalue", "_FDR")
            if q_col not in df.columns:
                df[q_col] = benjamini_hochberg(df[col].values)

        fdr_cols = [c.replace("_pvalue", "_FDR") for c in pval_cols]
        # Count how many metrics are FDR < threshold
        fdr_mask = df[fdr_cols] < args.fdr_threshold
        df["n_significant"] = fdr_mask.sum(axis=1)

    # Filter by min_significant
    n_before = len(df)
    df_filtered = df[df["n_significant"] >= args.min_significant].copy()
    n_after = len(df_filtered)
    log(f"[INFO] Filtered: {n_before} -> {n_after} (n_significant >= {args.min_significant})", args.log_file)

    if n_after == 0:
        log("[WARN] No windows pass filter — outputting empty intervals file", args.log_file)
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        interval_cols = ["chrom", "start", "end", "ID"]
        pd.DataFrame(columns=interval_cols).to_csv(args.out, sep="\t", index=False)
        map_out = args.out.replace("_intervals.tsv", "_interval_map.tsv")
        pd.DataFrame(columns=["ID", "chrom", "start", "end", "window_id", "n_significant"]).to_csv(
            map_out, sep="\t", index=False
        )
        sys.exit(0)

    # Assign ID = {chrom}_{window_id}
    df_filtered["ID"] = df_filtered["chrom"] + "_" + df_filtered["window_id"].astype(str)

    # Prepare interval_map: for joining meltron scores back to original windows
    interval_map = df_filtered[["chrom", "start", "end", "window_start", "window_end",
                                 "window_id", "ID", "n_significant"]].copy()
    interval_map.columns = ["chrom", "interval_start", "interval_end",
                             "window_start", "window_end", "window_id", "ID", "n_significant"]

    # Prepare strict 4-col Meltron interval output
    intervals_out = df_filtered[["chrom", "window_start", "window_end", "ID"]].copy()
    intervals_out.columns = ["chrom", "start", "end", "ID"]

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    intervals_out.to_csv(args.out, sep="\t", index=False)
    log(f"[INFO] Intervals written: {args.out} ({len(intervals_out)} rows)", args.log_file)

    # Write interval_map alongside
    map_out = args.out.replace("_intervals.tsv", "_interval_map.tsv")
    interval_map.to_csv(map_out, sep="\t", index=False)
    log(f"[INFO] Interval map written: {map_out}", args.log_file)


if __name__ == "__main__":
    run()
