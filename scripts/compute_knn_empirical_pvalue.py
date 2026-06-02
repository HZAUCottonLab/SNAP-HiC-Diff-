#!/usr/bin/env python3
"""
KNN-based empirical P-value computation for window pair metrics.

Usage:
    python scripts/compute_knn_empirical_pvalue.py \
        --comparison H1398_vs_L134 \
        --input results/window_compare/H1398_vs_L134/H1398_vs_L134_window_pair_metrics.tsv \
        --background /share/home/xhhuang24/project/PanTE/background_pool/results/window_compare/HiC_window_background.tsv \
        --out results/knn_empirical_pvalue/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv \
        --log logs/knn_empirical_pvalue/H1398_vs_L134/knn_empirical_pvalue.log
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from joblib import Parallel, delayed
from sklearn.neighbors import BallTree

warnings.filterwarnings("ignore", category=FutureWarning)


def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def derive_scores(df, derive_configs):
    """Derive score columns from raw metrics."""
    df = df.copy()
    for col, formula in derive_configs.items():
        if "max(tw_pvalue" in formula:
            df[col] = -np.log10(np.maximum(df["tw_pvalue"], 1e-30))
        elif "1 - " in formula:
            src = formula.split("1 - ")[1]
            df[col] = 1.0 - df[src]
    return df


def robust_scaler_fit(data):
    """Fit robust scaler: median + IQR."""
    med = np.nanmedian(data)
    q75, q25 = np.nanpercentile(data, [75, 25])
    iqr = q75 - q25
    return med, iqr


def robust_scaler_transform(data, med, iqr):
    """Transform using fitted scaler."""
    if iqr == 0:
        return np.zeros_like(data, dtype=float)
    return (data - med) / iqr


def _one_window(i, all_indices, all_dists, bg_scores, query_scores, score_cols):
    """Compute per-metric empirical p-values and bg medians for one window.
    
    Called in a thread-parallel context; all data arrays are numpy views on
    shared memory so no copy occurs under threading backend.
    """
    idx_i = all_indices[i]
    dists_i = all_dists[i]
    valid_k = len(idx_i)

    row = {
        "knn_valid_k": valid_k,
        "knn_dist_median": float(np.median(dists_i)) if valid_k > 0 else float("nan"),
        "knn_dist_max": float(dists_i[-1]) if valid_k > 0 else float("nan"),
    }
    for col in score_cols:
        bg_vals = bg_scores[col][idx_i]
        q_val = query_scores[col][i]
        count = int(np.sum(bg_vals >= q_val))
        pval = (count + 1) / (valid_k + 1)
        row[f"{col}_empirical_pvalue"] = pval
        row[f"bg_{col}_median"] = float(np.nanmedian(bg_vals))
    return row


def compute_knn_pvalues(query_df, bg_df, knn_k, n_jobs, score_cols, log_file):
    """Compute KNN-based empirical P-values, parallelised by window."""
    ec_med, ec_iqr = robust_scaler_fit(bg_df["EC"].values)
    ib_med, ib_iqr = robust_scaler_fit(bg_df["IB"].values)
    sp_med, sp_iqr = robust_scaler_fit(bg_df["SP"].values)

    log(f"[INFO] EC  median={ec_med:.6f}, IQR={ec_iqr:.6f}", log_file)
    log(f"[INFO] IB  median={ib_med:.6f}, IQR={ib_iqr:.6f}", log_file)
    log(f"[INFO] SP  median={sp_med:.6f}, IQR={sp_iqr:.6f}", log_file)

    bg_knn_features = np.column_stack([
        robust_scaler_transform(bg_df["EC"].values, ec_med, ec_iqr),
        robust_scaler_transform(bg_df["IB"].values, ib_med, ib_iqr),
        robust_scaler_transform(bg_df["SP"].values, sp_med, sp_iqr),
    ])
    tree = BallTree(bg_knn_features, metric="euclidean")
    log(f"[INFO] BallTree built: n={len(bg_knn_features)}, metric=euclidean, k={knn_k}", log_file)

    # Pre-compute scaled query features
    query_ec = robust_scaler_transform(query_df["EC"].values, ec_med, ec_iqr)
    query_ib = robust_scaler_transform(query_df["IB"].values, ib_med, ib_iqr)
    query_sp = robust_scaler_transform(query_df["SP"].values, sp_med, sp_iqr)
    query_features = np.column_stack([query_ec, query_ib, query_sp])

    # Batch KNN query (single call, not parallelised)
    n_windows = len(query_df)
    log(f"[INFO] Running {n_windows} KNN queries with k={knn_k}...", log_file)
    all_dists, all_indices = tree.query(query_features, k=knn_k)

    # Prepare shared numpy arrays for thread workers
    bg_scores = {col: bg_df[col].values for col in score_cols}
    query_scores = {col: query_df[col].values for col in score_cols}

    log(f"[INFO] Computing per-window p-values with {n_jobs} threads...", log_file)
    rows = Parallel(n_jobs=n_jobs, backend="threading")(
        delayed(_one_window)(
            i, all_indices, all_dists, bg_scores, query_scores, score_cols
        )
        for i in range(n_windows)
    )
    log(f"[INFO] Parallel p-value computation done.", log_file)

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Compute KNN-based empirical P-values for window pair metrics.")
    parser.add_argument("--comparison", required=True, help="Comparison name")
    parser.add_argument("--input", required=True, help="Input window_pair_metrics TSV file")
    parser.add_argument("--background", required=True, help="Background TSV file")
    parser.add_argument("--out", required=True, help="Output TSV file")
    parser.add_argument("--log", required=True, help="Log file path")
    args = parser.parse_args()

    Path(args.log).parent.mkdir(parents=True, exist_ok=True)

    with open("config/config.yaml") as f:
        knn_cfg = yaml.safe_load(f)["knn_empirical_pvalue"]

    knn_k = knn_cfg["knn_k"]
    score_cols = knn_cfg["score_columns"]
    n_jobs = knn_cfg.get("n_jobs", 80)

    log(f"[INFO] === KNN Empirical P-Value ===", args.log)
    log(f"[INFO] Comparison: {args.comparison}", args.log)
    log(f"[INFO] Input: {args.input}", args.log)
    log(f"[INFO] Background: {args.background}", args.log)
    log(f"[INFO] KNN k: {knn_k}", args.log)
    log(f"[INFO] Scaler: {knn_cfg['scaler']}", args.log)
    log(f"[INFO] Score columns: {score_cols}", args.log)
    log(f"[INFO] n_jobs: {n_jobs} (threading backend)", args.log)

    query_df = pd.read_csv(args.input, sep="\t")
    log(f"[INFO] Loaded input: {len(query_df)} rows", args.log)

    derive_configs = knn_cfg.get("derive_scores", {})
    query_df = derive_scores(query_df, derive_configs)
    log(f"[INFO] Derived scores: d_SSIM, d_Spearman, tw_score, d_NCC", args.log)

    bg_df = pd.read_csv(args.background, sep="\t")
    log(f"[INFO] Loaded background: {len(bg_df)} rows", args.log)
    bg_cols = ["EC", "IB", "SP", "MSE", "SSIM", "Spearman", "tw_pvalue", "NCC"]
    for col in bg_cols:
        if col not in bg_df.columns:
            raise ValueError(f"Background missing column: {col}")
    bg_df = bg_df[bg_cols].copy()
    for col, formula in derive_configs.items():
        if "max(tw_pvalue" in formula:
            bg_df[col] = -np.log10(np.maximum(bg_df["tw_pvalue"], 1e-30))
        elif "1 - " in formula:
            src = formula.split("1 - ")[1]
            bg_df[col] = 1.0 - bg_df[src]
    bg_df = bg_df.replace([np.inf, -np.inf], np.nan).dropna()
    log(f"[INFO] Background after cleaning: {len(bg_df)} rows", args.log)

    knn_results = compute_knn_pvalues(
        query_df, bg_df, knn_k, n_jobs=n_jobs, score_cols=score_cols, log_file=args.log
    )

    key_cols = ["chrom", "window_id", "window_start", "window_end", "comparison"]
    out_df = query_df[key_cols].copy()
    out_df = pd.concat([out_df.reset_index(drop=True), knn_results], axis=1)

    log(f"[INFO] Final output: {len(out_df)} rows, {len(out_df.columns)} columns", args.log)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    log(f"[INFO] Saved to {args.out}", args.log)


if __name__ == "__main__":
    main()
