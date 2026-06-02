#!/usr/bin/env python3
"""
Merge covariates, Gaussian metrics, and TW results into a single file.
Combines all window pair metrics for downstream analysis.

Usage:
    python scripts/merge_window_pair_metrics.py \
        --comparison H1398_rep1_vs_H1398_rep2 \
        --sample-a H1398_rep1 \
        --sample-b H1398_rep2 \
        --covariates results/window_compare/H1398_rep1_vs_H1398_rep2/H1398_rep1_vs_H1398_rep2_covariates.tsv \
        --gaussian results/window_compare/H1398_rep1_vs_H1398_rep2/H1398_rep1_vs_H1398_rep2_gaussian_metrics.tsv \
        --tw results/window_compare/H1398_rep1_vs_H1398_rep2/H1398_rep1_vs_H1398_rep2_tw.tsv \
        --out results/window_compare/H1398_rep1_vs_H1398_rep2/H1398_rep1_vs_H1398_rep2_window_pair_metrics.tsv \
        --log logs/window_compare/H1398_rep1_vs_H1398_rep2/merge_metrics.log
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd


def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def merge_window_pair_metrics(comparison, sample_a, sample_b, covariates_path, gaussian_path,
                             tw_path, output, log_file):
    """Merge all window pair metrics into a single DataFrame."""
    log(f"[INFO] Merging metrics for {comparison}", log_file)

    try:
        # Load all metrics
        cov_df = pd.read_csv(covariates_path, sep="\t")
        log(f"[INFO] Loaded covariates: {len(cov_df)} rows", log_file)

        gauss_df = pd.read_csv(gaussian_path, sep="\t")
        log(f"[INFO] Loaded Gaussian metrics: {len(gauss_df)} rows", log_file)

        tw_df = pd.read_csv(tw_path, sep="\t")
        log(f"[INFO] Loaded TW results: {len(tw_df)} rows", log_file)

        # Merge on common keys
        merge_keys = ["chrom", "window_id", "window_start", "window_end"]
        merged = cov_df.merge(gauss_df[merge_keys + ["MSE", "SSIM", "Spearman", "NCC"]],
                             on=merge_keys, how="left")
        merged = merged.merge(tw_df[merge_keys + ["tw_statistic", "tw_pvalue"]],
                            on=merge_keys, how="left")

        # Add comparison metadata
        merged["comparison"] = comparison
        merged["sample_a"] = sample_a
        merged["sample_b"] = sample_b

        # Sort by chromosome and position
        merged = merged.sort_values(["chrom", "window_start", "window_end"])

        Path(output).parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(output, sep="\t", index=False)

        log(f"[INFO] Saved merged metrics to {output}", log_file)
        log(f"[INFO] Total rows: {len(merged)}", log_file)

        # Summary statistics
        sig_count = (merged["tw_pvalue"] < 0.05).sum() if "tw_pvalue" in merged.columns else 0
        log(f"[INFO] Significant windows (p < 0.05): {sig_count}", log_file)

        return True

    except Exception as e:
        log(f"[ERROR] Failed to merge metrics: {e}", log_file)
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(description="Merge window pair metrics.")
    parser.add_argument("--comparison", required=True, help="Comparison name")
    parser.add_argument("--sample-a", required=True, help="Sample A name")
    parser.add_argument("--sample-b", required=True, help="Sample B name")
    parser.add_argument("--covariates", required=True, help="Covariates TSV file")
    parser.add_argument("--gaussian", required=True, help="Gaussian metrics TSV file")
    parser.add_argument("--tw", required=True, help="Tracy-Widom results TSV file")
    parser.add_argument("--out", required=True, help="Output merged TSV file")
    parser.add_argument("--log", required=True, help="Log file path")
    args = parser.parse_args()

    Path(args.log).parent.mkdir(parents=True, exist_ok=True)

    success = merge_window_pair_metrics(
        comparison=args.comparison,
        sample_a=args.sample_a,
        sample_b=args.sample_b,
        covariates_path=args.covariates,
        gaussian_path=args.gaussian,
        tw_path=args.tw,
        output=args.out,
        log_file=args.log
    )

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
