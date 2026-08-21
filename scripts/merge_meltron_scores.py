#!/usr/bin/env python3
"""
Merge Meltron decondensation and condensation results into one TSV.

Outputs:
    - meltron_merged.tsv: chrom, start, end, ID, meltron_melting_score,
      meltron_compaction_score, meltron_best_score, meltron_best_direction,
      meltron_is_significant, window_id, n_significant_metrics

Usage:
    python scripts/merge_meltron_scores.py \
        --comparison Mutant1_vs_Mutant2 \
        --decondensation results/meltron/Mutant1_vs_Mutant2/Mutant1_vs_Mutant2_decondensation.tsv \
        --condensation results/meltron/Mutant1_vs_Mutant2/Mutant1_vs_Mutant2_condensation.tsv \
        --interval-map results/meltron/Mutant1_vs_Mutant2/Mutant1_vs_Mutant2_interval_map.tsv \
        --out results/meltron/Mutant1_vs_Mutant2/Mutant1_vs_Mutant2_meltron_merged.tsv \
        --score-threshold 5.0 \
        --log-file results/logs/meltron/Mutant1_vs_Mutant2.log
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


SCALE_COLUMNS = [
    "resolution_bp",
    "window_bins",
    "window_size_bp",
    "step_bins",
    "step_size_bp",
]

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

MERGED_COLUMNS = [
    *INTERVAL_MAP_COLUMNS,
    "meltron_melting_score",
    "meltron_compaction_score",
    "meltron_best_score",
    "meltron_best_direction",
    "meltron_is_significant",
    "perc_na_or_missing_ref_decondensation",
    "perc_na_or_missing_inp_decondensation",
    "perc_na_or_missing_ref_condensation",
    "perc_na_or_missing_inp_condensation",
]

def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def normalize_directional_score(frame, source_column, normalized_column, direction, parser):
    if frame.empty:
        return frame
    if source_column not in frame.columns:
        parser.error(f"{direction} result requires {source_column}")
    normalized = frame.copy()
    normalized[normalized_column] = normalized[source_column]
    return normalized


def normalize_directional_na_columns(frame, direction):
    return frame.rename(
        columns={
            prefix: f"{prefix}_{direction}"
            for prefix in ("perc_na_or_missing_ref", "perc_na_or_missing_inp")
            if prefix in frame.columns
        }
    )


def read_interval_map(path, parser):
    if not path or not os.path.exists(path):
        return pd.DataFrame(columns=INTERVAL_MAP_COLUMNS)

    interval_map = pd.read_csv(path, sep="\t")
    metric_count_column = (
        "n_significant_metrics"
        if "n_significant_metrics" in interval_map.columns
        else "n_significant"
    )
    if metric_count_column not in interval_map.columns:
        parser.error("interval map requires n_significant or n_significant_metrics")
    interval_map = interval_map.rename(
        columns={metric_count_column: "n_significant_metrics"}
    ).copy()
    for column in INTERVAL_MAP_COLUMNS:
        if column not in interval_map.columns:
            interval_map[column] = np.nan
    return interval_map.loc[:, INTERVAL_MAP_COLUMNS]

def run():
    parser = argparse.ArgumentParser(
        description="Merge Meltron decondensation and condensation results."
    )
    parser.add_argument(
        "--comparison", required=True,
        help="Comparison name"
    )
    parser.add_argument(
        "--decondensation", required=True,
        help="Path to decondensation result TSV"
    )
    parser.add_argument(
        "--condensation", required=True,
        help="Path to condensation result TSV"
    )
    parser.add_argument(
        "--interval-map",
        help="Path to interval_map.tsv from prepare_meltron_intervals.py"
    )
    parser.add_argument(
        "--out", required=True,
        help="Path to merged output TSV"
    )
    parser.add_argument(
        "--score-threshold", type=float, default=5.0,
        help="Advisory threshold for meltron_is_significant flag (default: 5.0)"
    )
    parser.add_argument(
        "--log-file",
        help="Log file path"
    )
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # Read decondensation and condensation
    decond = pd.read_csv(args.decondensation, sep="\t") if os.path.exists(args.decondensation) else pd.DataFrame()
    condens = pd.read_csv(args.condensation, sep="\t") if os.path.exists(args.condensation) else pd.DataFrame()

    if decond.empty and condens.empty:
        log("[WARN] Both Meltron results are empty — writing empty output", args.log_file)
        pd.DataFrame(columns=MERGED_COLUMNS).to_csv(args.out, sep="\t", index=False)
        sys.exit(0)

    decond = normalize_directional_score(
        decond, "melting_score", "_meltron_decond_score", "decondensation", parser
    )
    condens = normalize_directional_score(
        condens, "compaction_score", "_meltron_condens_score", "condensation", parser
    )

    if decond.empty:
        log("[WARN] Decondensation result is empty — condensation only", args.log_file)
        merged = normalize_directional_na_columns(condens.copy(), "condensation")
        merged["meltron_melting_score"] = np.nan
        merged["meltron_compaction_score"] = merged["_meltron_condens_score"]
        merged["meltron_best_score"] = merged["meltron_compaction_score"]
        merged["meltron_best_direction"] = "condensation"
    elif condens.empty:
        log("[WARN] Condensation result is empty — decondensation only", args.log_file)
        merged = normalize_directional_na_columns(decond.copy(), "decondensation")
        merged["meltron_compaction_score"] = np.nan
        merged["meltron_melting_score"] = merged["_meltron_decond_score"]
        merged["meltron_best_score"] = merged["meltron_melting_score"]
        merged["meltron_best_direction"] = "decondensation"
    else:
        # Directional scores use unique internal columns before this outer merge.
        # The directional NA metrics still receive merge suffixes below.
        key_cols = ["chrom", "start", "end", "ID"]
        merged = decond.merge(condens, on=key_cols, how="outer", suffixes=("_decond", "_condens"))

        merged["meltron_melting_score"] = merged["_meltron_decond_score"]
        merged["meltron_compaction_score"] = merged["_meltron_condens_score"]

        melting_score = merged["meltron_melting_score"]
        compaction_score = merged["meltron_compaction_score"]
        merged["meltron_best_score"] = merged[["meltron_melting_score", "meltron_compaction_score"]].max(axis=1)
        has_scores = melting_score.notna() & compaction_score.notna()
        best_direction = pd.Series(pd.NA, index=merged.index, dtype="object")
        best_direction.loc[melting_score.notna() & compaction_score.isna()] = "decondensation"
        best_direction.loc[melting_score.isna() & compaction_score.notna()] = "condensation"
        best_direction.loc[has_scores & (melting_score >= compaction_score)] = "decondensation"
        best_direction.loc[has_scores & (melting_score < compaction_score)] = "condensation"
        merged["meltron_best_direction"] = best_direction

    # Advisory significance flag
    merged["meltron_is_significant"] = (merged["meltron_best_score"] >= args.score_threshold).astype(int)

    # The interval map is the authoritative adapter for source-window identity.
    provenance = read_interval_map(args.interval_map, parser)
    merged = merged.merge(
        provenance,
        on="ID",
        how="left",
        validate="one_to_one",
        suffixes=("_score", ""),
    )
    merged["chrom"] = merged["chrom"].combine_first(merged["chrom_score"])

    # Rename perc_na columns for clarity in output
    for suffix, direction in [("_decond", "decondensation"), ("_condens", "condensation")]:
        for prefix in ["perc_na_or_missing_ref", "perc_na_or_missing_inp"]:
            col = prefix + suffix
            if col in merged.columns:
                new_col = prefix + "_" + direction
                merged.rename(columns={col: new_col}, inplace=True)

    # Emit a stable schema even when one direction has no NA metrics.
    for column in MERGED_COLUMNS:
        if column not in merged.columns:
            merged[column] = np.nan
    merged = merged.loc[:, MERGED_COLUMNS]

    merged.to_csv(args.out, sep="\t", index=False)
    log(f"[INFO] Merged Meltron output: {args.out} ({len(merged)} rows)", args.log_file)


if __name__ == "__main__":
    run()
