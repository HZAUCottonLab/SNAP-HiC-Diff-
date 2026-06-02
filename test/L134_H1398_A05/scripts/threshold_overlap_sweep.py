#!/usr/bin/env python3
"""
Threshold sweep — evaluate all 186 metric+threshold+mode combinations on A05.

Outputs:
  results/threshold_overlap_summary.tsv  — one row per threshold_id
  results/threshold_overlap_detail.tsv   — one row per (threshold_id, window)
  logs/threshold_overlap_test.log        — run log
"""

import argparse
import os
import sys
import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

warnings.filterwarnings("ignore")

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REPO_ROOT = os.path.dirname(PROJECT_ROOT)
TEST_ROOT = os.path.dirname(os.path.abspath(__file__))
TEST_ROOT = os.path.dirname(TEST_ROOT)


PVALUE_COLS = {
    "MSE":       "MSE_empirical_pvalue",
    "d_SSIM":    "d_SSIM_empirical_pvalue",
    "d_Spearman":"d_Spearman_empirical_pvalue",
    "tw_score":  "tw_score_empirical_pvalue",
    "d_NCC":     "d_NCC_empirical_pvalue",
}
FDR_COLS = {k: v.replace("pvalue", "FDR") for k, v in PVALUE_COLS.items()}

METRICS = list(PVALUE_COLS.keys())
THRESHOLDS = [0.01, 0.05]
FILTER_MODES = ["pvalue", "fdr_chrom", "fdr_global"]


def bh_fdr(pvalues):
    """Benjamini-Hochberg FDR correction."""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    idx = np.argsort(p, kind="quicksort")
    sorted_p = p[idx]
    adjusted = np.empty(n)
    for rank in range(n):
        adjusted[idx[rank]] = min(sorted_p[rank] * n / (rank + 1), 1.0)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    return np.clip(adjusted, 0.0, 1.0)


def setup_log(log_path):
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.FileHandler(log_path, mode="w"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger(__name__)


def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg FDR correction — returns adjusted array in original order."""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    idx = np.argsort(p, kind="quicksort")
    sorted_p = p[idx]
    adjusted = np.empty(n)
    for rank in range(n):
        bh_value = sorted_p[rank] * n / (rank + 1)
        adjusted[rank] = min(bh_value, 1.0)
    adjusted[::-1] = np.minimum.accumulate(adjusted[::-1])
    adjusted = np.clip(adjusted, 0.0, 1.0)
    reverse_idx = np.argsort(idx)
    return adjusted[reverse_idx]


def load_enriched(enriched_path, chrom):
    """Load and filter to target chromosome."""
    df = pd.read_csv(enriched_path, sep="\t")
    df = df[df["chrom"] == chrom].copy()
    df = df.sort_values(["chrom", "window_id"]).reset_index(drop=True)
    return df


def load_reference(ref_tsv_path):
    """Load reference windows as set of (chrom, window_start, window_end)."""
    ref = pd.read_csv(ref_tsv_path, sep="\t")
    ref_set = set(zip(ref["chrom"], ref["start"], ref["end"]))
    return ref_set, ref


def add_fdr_columns(df, global_enriched_path, logger=None):
    """Add two FDR columns: chrom-internal and global."""
    pval_matrix = df[list(PVALUE_COLS.values())].values.astype(float)

    # Chrom-internal FDR (A05 only)
    fdr_chrom = np.column_stack([
        bh_fdr(pval_matrix[:, i]) for i in range(pval_matrix.shape[1])
    ])
    for i, m in enumerate(PVALUE_COLS.keys()):
        df[FDR_COLS[m] + "_chrom"] = fdr_chrom[:, i]

    # Global FDR (all chromosomes)
    usecols = ["chrom", "window_id"] + list(PVALUE_COLS.values())
    global_df = pd.read_csv(global_enriched_path, sep="\t", usecols=usecols)
    if logger:
        logger.info("  Loaded global enriched (cols=%d, rows=%d)", len(usecols), len(global_df))

    global_pval_matrix = global_df[list(PVALUE_COLS.values())].values.astype(float)
    fdr_global = np.column_stack([
        bh_fdr(global_pval_matrix[:, i]) for i in range(global_pval_matrix.shape[1])
    ])

    global_fdr_df = global_df[["chrom", "window_id"]].copy()
    for i, m in enumerate(PVALUE_COLS.keys()):
        global_fdr_df[FDR_COLS[m] + "_global"] = fdr_global[:, i]

    df = df.merge(global_fdr_df, on=["chrom", "window_id"], how="left")
    return df


def load_grid_yaml(grid_yaml_path):
    with open(grid_yaml_path, "r") as f:
        cfg = yaml.safe_load(f)
    combos = cfg.get("combinations", [])
    if not combos:
        raise ValueError(f"No combinations found in {grid_yaml_path}")
    return cfg, combos


def validate_combo(combo):
    required = ["threshold_id", "metrics", "threshold", "filter_mode", "n_metrics_required"]
    for k in required:
        if k not in combo:
            raise ValueError(f"Missing '{k}' in combo: {combo}")
    for m in combo["metrics"]:
        if m not in METRICS:
            raise ValueError(f"Unknown metric '{m}' in combo {combo['threshold_id']}")
    if combo["filter_mode"] not in FILTER_MODES:
        raise ValueError(f"Unknown filter_mode '{combo['filter_mode']}'")


def get_filter_col(metric, mode):
    """Return the column name to use for filtering."""
    if mode == "pvalue":
        return PVALUE_COLS[metric]
    if mode == "fdr_chrom":
        return FDR_COLS[metric] + "_chrom"
    if mode == "fdr_global":
        return FDR_COLS[metric] + "_global"
    raise ValueError(f"Unknown filter_mode: {mode}")


def apply_threshold(df, metrics, threshold, filter_mode):
    """Return mask of windows passing threshold on ALL specified metrics."""
    cols = [get_filter_col(m, filter_mode) for m in metrics]
    mat = df[cols].values.astype(float)
    return np.all(mat < float(threshold), axis=1)


def compute_overlap_stats(pass_df, ref_set, ref_df, n_total_windows):
    """Compute per-threshold overlap statistics."""
    n_pass = len(pass_df)
    n_ref = len(ref_df)
    pass_windows = set(zip(pass_df["chrom"], pass_df["window_start"], pass_df["window_end"]))
    overlap_windows = pass_windows & ref_set
    n_overlap = len(overlap_windows)

    pass_overlap_ratio = n_overlap / n_pass if n_pass > 0 else 0.0
    reference_recall = n_overlap / n_ref if n_ref > 0 else 0.0

    precision = pass_overlap_ratio
    recall = reference_recall
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    expected = n_pass * n_ref / n_total_windows if n_total_windows > 0 else 0.0
    enrichment = n_overlap / expected if expected > 0 else float("inf")

    win_len = (pass_df["window_end"] - pass_df["window_start"]).astype(float)
    overlap_bp_list = win_len[pass_df.apply(
        lambda r: (r["chrom"], r["window_start"], r["window_end"]) in ref_set,
        axis=1
    )].values

    mean_overlap_bp = float(np.mean(overlap_bp_list)) if len(overlap_bp_list) else float("nan")
    median_overlap_bp = float(np.median(overlap_bp_list)) if len(overlap_bp_list) else float("nan")
    total_overlap_bp = float(np.sum(overlap_bp_list)) if len(overlap_bp_list) else 0.0

    return {
        "n_pass_windows": n_pass,
        "n_reference_windows": n_ref,
        "n_overlap_windows": n_overlap,
        "pass_overlap_ratio": round(pass_overlap_ratio, 6),
        "reference_recall": round(reference_recall, 6),
        "precision": round(precision, 6),
        "recall": round(recall, 6),
        "f1_score": round(f1, 6),
        "enrichment_score": round(enrichment, 4),
        "mean_overlap_bp": round(mean_overlap_bp, 1),
        "median_overlap_bp": round(median_overlap_bp, 1),
        "total_overlap_bp": round(total_overlap_bp, 1),
    }


def build_detail_rows(df, metrics, threshold, filter_mode, ref_set, threshold_id):
    """Build detail rows for each pass window."""
    mask = apply_threshold(df, metrics, threshold, filter_mode)
    pass_df = df[mask].copy()
    rows = []
    for _, row in pass_df.iterrows():
        is_overlap = (row["chrom"], row["window_start"], row["window_end"]) in ref_set
        row_data = {
            "threshold_id": threshold_id,
            "metrics": ",".join(metrics),
            "threshold": threshold,
            "filter_mode": filter_mode,
            "n_metrics_required": len(metrics),
            "chrom": row["chrom"],
            "window_id": row["window_id"],
            "window_start": row["window_start"],
            "window_end": row["window_end"],
            "overlaps_reference": 1 if is_overlap else 0,
        }
        for m in METRICS:
            if filter_mode == "pvalue":
                row_data[f"{m}_val"] = row.get(PVALUE_COLS[m], np.nan)
            elif filter_mode == "fdr_chrom":
                row_data[f"{m}_val"] = row.get(FDR_COLS[m] + "_chrom", np.nan)
            else:
                row_data[f"{m}_val"] = row.get(FDR_COLS[m] + "_global", np.nan)
        rows.append(row_data)
    return rows


def run_sweep(args):
    log = setup_log(args.log_file)
    log.info("=== Threshold overlap sweep ===")
    log.info("Project: %s | Chrom: %s", args.project, args.chrom)
    log.info("Enriched: %s", args.enriched)
    log.info("Reference: %s", args.reference)
    log.info("Output dir: %s", args.out_dir)

    # Load data
    log.info("Loading enriched windows (chrom=%s)...", args.chrom)
    df = load_enriched(args.enriched, args.chrom)
    log.info("  A05 windows: %d", len(df))

    log.info("Computing FDR columns (chrom + global)...")
    df = add_fdr_columns(df, args.global_enriched, logger=log)
    n_total_windows = len(df)

    log.info("Loading reference windows...")
    ref_set, ref_df = load_reference(args.reference)
    log.info("  Reference windows: %d", len(ref_set))

    # Load threshold grid
    cfg, combos = load_grid_yaml(args.grid_yaml)
    for c in combos:
        validate_combo(c)

    log.info("Total threshold combinations: %d", len(combos))

    summary_rows = []
    detail_rows = []

    for combo in combos:
        tid = int(combo["threshold_id"])
        metrics = list(combo["metrics"])
        thresh = float(combo["threshold"])
        mode = combo["filter_mode"]

        mask = apply_threshold(df, metrics, thresh, mode)
        pass_df = df[mask]
        stats = compute_overlap_stats(pass_df, ref_set, ref_df, n_total_windows=n_total_windows)
        detail_rows.extend(build_detail_rows(df, metrics, thresh, mode, ref_set, tid))

        summary_row = {
            "threshold_id": tid,
            "n_metrics_required": int(combo["n_metrics_required"]),
            "metrics": ",".join(metrics),
            "threshold": thresh,
            "filter_mode": mode,
            **stats,
        }
        summary_rows.append(summary_row)
        log.info(
            "[%03d] %s | thresh=%.2f | mode=%-11s | pass=%4d | overlap=%2d | recall=%.3f | f1=%.4f",
            tid, ",".join(metrics), thresh, mode,
            stats["n_pass_windows"], stats["n_overlap_windows"],
            stats["reference_recall"], stats["f1_score"],
        )

    # Write outputs
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_path = out_dir / "threshold_overlap_summary.tsv"
    detail_path = out_dir / "threshold_overlap_detail.tsv"

    summary_df = pd.DataFrame(summary_rows)
    # Reorder columns
    summary_cols = [
        "threshold_id", "n_metrics_required", "metrics", "threshold", "filter_mode",
        "n_pass_windows", "n_reference_windows", "n_overlap_windows",
        "pass_overlap_ratio", "reference_recall", "precision", "recall",
        "f1_score", "enrichment_score",
        "mean_overlap_bp", "median_overlap_bp", "total_overlap_bp",
    ]
    summary_df = summary_df[summary_cols]
    summary_df.to_csv(summary_path, sep="\t", index=False)
    log.info("Written: %s (%d rows)", summary_path, len(summary_df))

    detail_df = pd.DataFrame(detail_rows)
    detail_cols = [
        "threshold_id", "n_metrics_required", "metrics", "threshold", "filter_mode",
        "chrom", "window_id", "window_start", "window_end",
        "overlaps_reference",
    ] + [f"{m}_val" for m in METRICS]
    available_detail_cols = [c for c in detail_cols if c in detail_df.columns]
    detail_df = detail_df[available_detail_cols]
    detail_df.to_csv(detail_path, sep="\t", index=False)
    log.info("Written: %s (%d rows)", detail_path, len(detail_df))
    log.info("=== Done ===")


def main():
    parser = argparse.ArgumentParser(description="Threshold sweep — overlap evaluation")
    parser.add_argument("--project",    default="L134_H1398_A05")
    parser.add_argument("--chrom",       default="HC04_A05")
    parser.add_argument("--comparison",  default="H1398_vs_L134")
    parser.add_argument(
        "--grid-yaml",
        default="test/L134_H1398_A05/config/threshold_grid.yaml",
        help="Threshold grid YAML (relative to repo root or absolute)",
    )
    parser.add_argument(
        "--global-enriched",
        default="results/knn_empirical_pvalue/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv",
        help="Path to knn_enriched.tsv for global FDR (relative to repo root or absolute)",
    )
    parser.add_argument(
        "--enriched",
        default="results/knn_empirical_pvalue/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv",
        help="Path to knn_enriched.tsv (relative to repo root or absolute)",
    )
    parser.add_argument(
        "--reference",
        default="test/L134_H1398_A05/reference/reference_windows.tsv",
        help="Path to reference windows TSV",
    )
    parser.add_argument(
        "--out-dir",
        default="test/L134_H1398_A05/results",
    )
    parser.add_argument(
        "--log-file",
        default="test/L134_H1398_A05/logs/threshold_overlap_test.log",
    )
    args = parser.parse_args()

    # Resolve relative paths from repo root
    def resolve(p):
        if os.path.isabs(p):
            return p
        return os.path.join(REPO_ROOT, p)

    args.enriched = resolve(args.enriched)
    args.global_enriched = resolve(args.global_enriched)
    args.grid_yaml = resolve(args.grid_yaml)
    args.reference = resolve(args.reference)
    args.out_dir = resolve(args.out_dir)
    args.log_file = resolve(args.log_file)

    # Ensure log dir exists
    Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)

    run_sweep(args)


if __name__ == "__main__":
    main()
