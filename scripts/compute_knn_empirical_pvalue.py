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
import os
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from joblib import Parallel, delayed
from sklearn.neighbors import BallTree

warnings.filterwarnings("ignore", category=FutureWarning)


SCALE_COLUMNS = [
    "resolution_bp",
    "window_bins",
    "window_size_bp",
    "step_bins",
    "step_size_bp",
]
WINDOW_COLUMNS = ["chrom", "window_id", "window_start", "window_end", "comparison"]


def attach_scale_metadata(frame, resolution_bp, window_bins, step_bins):
    values = (resolution_bp, window_bins, step_bins)
    if any(type(value) is not int or value <= 0 for value in values):
        raise ValueError("resolution_bp, window_bins, and step_bins must be positive integers")
    result = frame.copy()
    result["resolution_bp"] = resolution_bp
    result["window_bins"] = window_bins
    result["window_size_bp"] = resolution_bp * window_bins
    result["step_bins"] = step_bins
    result["step_size_bp"] = resolution_bp * step_bins
    return result


def empty_knn_output(query_df, score_cols):
    result = query_df[[*SCALE_COLUMNS, *WINDOW_COLUMNS]].copy()
    for column in ("knn_valid_k", "knn_dist_median", "knn_dist_max"):
        result[column] = pd.Series(index=result.index, dtype=float)
    for score in score_cols:
        result[f"{score}_empirical_pvalue"] = pd.Series(
            index=result.index, dtype=float
        )
        result[f"bg_{score}_median"] = pd.Series(index=result.index, dtype=float)
    return result


def atomic_tsv(frame, output):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp"
    )
    os.close(descriptor)
    try:
        frame.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, output)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


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


def _knn_query_chunk(chunk_features, bg_features, knn_k):
    """Run BallTree.query on one chunk of query features.
    
    Each loky worker builds its own BallTree from bg_features, so we avoid
    passing the large tree object across process boundaries as an argument.
    Instead bg_features is captured via closure in the parent scope and each
    worker gets its own copy when the subprocess is forked/spawned.
    """
    tree = BallTree(bg_features, metric="euclidean")
    dists, indices = tree.query(chunk_features, k=knn_k)
    return dists, indices


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

    # Build background KNN features
    bg_knn_features = np.column_stack([
        robust_scaler_transform(bg_df["EC"].values, ec_med, ec_iqr),
        robust_scaler_transform(bg_df["IB"].values, ib_med, ib_iqr),
        robust_scaler_transform(bg_df["SP"].values, sp_med, sp_iqr),
    ])
    log(f"[INFO] BallTree feature matrix: n={len(bg_knn_features)}, k={knn_k}", log_file)

    # Pre-compute scaled query features
    query_ec = robust_scaler_transform(query_df["EC"].values, ec_med, ec_iqr)
    query_ib = robust_scaler_transform(query_df["IB"].values, ib_med, ib_iqr)
    query_sp = robust_scaler_transform(query_df["SP"].values, sp_med, sp_iqr)
    query_features = np.column_stack([query_ec, query_ib, query_sp])

    # Parallel KNN query: chunk windows and dispatch to loky workers
    n_windows = len(query_df)
    chunk_size = max(1, n_windows // n_jobs)
    n_chunks = max(1, (n_windows + chunk_size - 1) // chunk_size)

    log(
        f"[INFO] Parallel KNN query: {n_windows} windows in {n_chunks} chunks "
        f"(~{chunk_size}/chunk) with {n_jobs} loky workers...",
        log_file,
    )

    # Build chunks of query feature arrays (views are fine, loky copies on fork)
    chunk_results = Parallel(n_jobs=n_jobs, backend="loky", max_nbytes="1G")(
        delayed(_knn_query_chunk)(
            query_features[i : i + chunk_size],
            bg_knn_features,
            knn_k,
        )
        for i in range(0, n_windows, chunk_size)
    )

    # Concatenate per-chunk results
    all_dists = np.vstack([r[0] for r in chunk_results])
    all_indices = np.vstack([r[1] for r in chunk_results])
    log(f"[INFO] KNN query done. result shape: dists={all_dists.shape}, indices={all_indices.shape}", log_file)

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
    parser.add_argument("--config", required=True, help="Workflow YAML configuration")
    parser.add_argument("--resolution-bp", required=True, type=int)
    parser.add_argument("--window-bins", required=True, type=int)
    parser.add_argument("--step-bins", required=True, type=int)
    args = parser.parse_args()

    Path(args.log).parent.mkdir(parents=True, exist_ok=True)
    with open(args.config) as handle:
        knn_cfg = yaml.safe_load(handle)["knn_empirical_pvalue"]

    knn_k = knn_cfg["knn_k"]
    score_cols = knn_cfg["score_columns"]
    n_jobs = knn_cfg.get("n_jobs", 80)

    query_df = pd.read_csv(args.input, sep="\t")
    query_df = derive_scores(query_df, knn_cfg.get("derive_scores", {}))
    query_df = attach_scale_metadata(
        query_df, args.resolution_bp, args.window_bins, args.step_bins
    )
    if query_df.empty:
        atomic_tsv(empty_knn_output(query_df, score_cols), args.out)
        log("[INFO] Input contains zero windows; wrote a valid empty output", args.log)
        return

    bg_df = pd.read_csv(args.background, sep="\t")
    bg_cols = ["EC", "IB", "SP", "MSE", "SSIM", "Spearman", "tw_pvalue", "NCC"]
    missing = sorted(set(bg_cols).difference(bg_df.columns))
    if missing:
        raise ValueError(f"Background missing columns: {', '.join(missing)}")
    bg_df = derive_scores(
        bg_df[bg_cols].copy(), knn_cfg.get("derive_scores", {})
    )
    bg_df = bg_df.replace([np.inf, -np.inf], np.nan).dropna()
    if len(bg_df) < knn_k:
        raise ValueError(
            f"Background has {len(bg_df)} valid rows, smaller than knn_k {knn_k}"
        )

    knn_results = compute_knn_pvalues(
        query_df, bg_df, knn_k, n_jobs=n_jobs, score_cols=score_cols, log_file=args.log
    )
    out_df = pd.concat(
        [
            query_df[[*SCALE_COLUMNS, *WINDOW_COLUMNS]].reset_index(drop=True),
            knn_results.reset_index(drop=True),
        ],
        axis=1,
    )
    atomic_tsv(out_df, args.out)
    log(f"[INFO] Saved {len(out_df)} rows to {args.out}", args.log)


if __name__ == "__main__":
    main()
