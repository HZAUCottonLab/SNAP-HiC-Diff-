#!/usr/bin/env python3
"""
Compute Tracy-Widom test statistics for ICE window pairs.
Fetches ICE-balanced contact matrices from .cool files, computes
Diffmat = mat_A / mat_B, then applies Tracy-Widom test.

Usage:
    python scripts/compare_iced_windows_tw.py \
        --comparison H1398_rep1_vs_H1398_rep2 \
        --sample-a H1398_rep1 \
        --sample-b H1398_rep2 \
        --chroms HC04_A01 HC04_A02 \
        --covariates results/window_compare/H1398_rep1_vs_H1398_rep2/H1398_rep1_vs_H1398_rep2_covariates.tsv \
        --out results/window_compare/H1398_rep1_vs_H1398_rep2/H1398_rep1_vs_H1398_rep2_tw.tsv \
        --outdir results/window_compare/H1398_rep1_vs_H1398_rep2 \
        --resolution 4000 \
        --min-nbin 10 \
        --filter-threshold 0.5 \
        --n-jobs 80 \
        --log logs/window_compare/H1398_rep1_vs_H1398_rep2/compare_tw.log
"""

import argparse
import os
import sys
from pathlib import Path

import cooler
import numpy as np
import pandas as pd
from TracyWidom import TracyWidom
from joblib import Parallel, delayed

# Initialize Tracy-Widom distribution (beta=1 for real symmetric matrices)
TW_DIST = TracyWidom(beta=1)


def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def normdiff_by_mean_sd(mat):
    """
    Normalize off-diagonal entries by their row-wise mean and SD.
    For each off-diagonal element D[i,j] (i != j):
        D_norm[i,j] = (D[i,j] - mean_k) / std_k
    where k = |i - j| (distance from diagonal).
    """
    n = mat.shape[0]
    result = np.full_like(mat, np.nan)

    for k in range(1, n):
        off_diag = np.full(n - k, np.nan)
        valid = np.isfinite(mat[np.arange(k, n), np.arange(n - k)])
        off_diag[valid] = mat[np.arange(k, n), np.arange(n - k)][valid]

        valid_mask = np.isfinite(off_diag)
        if np.sum(valid_mask) < 2:
            continue

        vals = off_diag[valid_mask]
        mu = np.mean(vals)
        sigma = np.std(vals, ddof=1)

        if sigma > 0:
            normalized = (off_diag - mu) / sigma
            result[np.arange(k, n), np.arange(n - k)] = normalized
            result[np.arange(n - k), np.arange(k, n)] = normalized

    return result


def compute_tracy_widom_stat(mat):
    """
    Compute Tracy-Widom statistic for a normalized difference matrix.
    Steps:
        1. D = log(Diffmat)  -- off-diagonal already normalized by normdiff_by_mean_sd
        2. D_normalized = D / sqrt(n)
        3. TW_stat = n^(2/3) * (max_eigenvalue(D_normalized) - 2)
    """
    n = mat.shape[0]
    valid_mask = np.isfinite(mat)

    row_valid = np.sum(valid_mask, axis=1) >= n * 0.5
    col_valid = np.sum(valid_mask, axis=0) >= n * 0.5
    keep = row_valid & col_valid

    if np.sum(keep) < max(10, n * 0.3):
        return np.nan, 1.0, "INSUFFICIENT_VALID"

    sub_mat = mat[np.ix_(keep, keep)]
    sub_n = sub_mat.shape[0]

    # Replace any remaining NaN with 0 for eigenvalue computation
    nan_mask = np.isnan(sub_mat)
    if np.any(nan_mask):
        sub_mat = np.where(nan_mask, 0.0, sub_mat)

    try:
        eigenvalues = np.linalg.eigvalsh(sub_mat / np.sqrt(sub_n))
        lambda_max = np.max(np.abs(eigenvalues))
        tw_stat = (sub_n ** (2.0 / 3.0)) * (lambda_max - 2.0)
    except Exception:
        return np.nan, 1.0, "EIGENVALUE_ERROR"

    p_value = 1.0 - TW_DIST.cdf(tw_stat)

    return float(tw_stat), float(p_value), "OK"


def fetch_ice_matrix(sample, chrom, win_start, win_end, resolution, outdir):
    """
    Fetch ICE-balanced matrix from .cool file for a genomic window.
    Returns symmetric matrix with NaN for masked/gap bins.
    """
    cool_path = f"{outdir}/balanced/{sample}/{sample}__{chrom}_iced.cool"

    clr = cooler.Cooler(cool_path)
    mat = clr.matrix(balance=True, sparse=False).fetch(f"{chrom}:{win_start}-{win_end}")

    mat = np.asarray(mat, dtype=np.float64)
    mat = (mat + mat.T) / 2.0
    return mat


def process_single_window(row, chrom, sample_a, sample_b, resolution, outdir, min_nbin, na_fraction, log_file):
    """Process a single window: fetch ICE matrices, compute Diffmat, run TW."""
    window_id = row["window_id"]
    win_start = int(row["window_start"])
    win_end = int(row["window_end"])

    try:
        mat_a = fetch_ice_matrix(sample_a, chrom, win_start, win_end, resolution, outdir)
        mat_b = fetch_ice_matrix(sample_b, chrom, win_start, win_end, resolution, outdir)
    except Exception as e:
        log(f"[WARN] Failed to fetch matrices for window {window_id}: {e}", log_file)
        return {
            "chrom": chrom,
            "window_id": window_id,
            "window_start": win_start,
            "window_end": win_end,
            "tw_statistic": np.nan,
            "tw_pvalue": np.nan,
            "n_bins": 0,
            "status": "FETCH_ERROR",
        }

    if mat_a.shape != mat_b.shape:
        min_n = min(mat_a.shape[0], mat_b.shape[0])
        mat_a = mat_a[:min_n, :min_n]
        mat_b = mat_b[:min_n, :min_n]

    n = mat_a.shape[0]
    if n < min_nbin:
        return {
            "chrom": chrom,
            "window_id": window_id,
            "window_start": win_start,
            "window_end": win_end,
            "tw_statistic": np.nan,
            "tw_pvalue": np.nan,
            "n_bins": n,
            "status": "TOO_SMALL",
        }

    with np.errstate(divide="ignore", invalid="ignore"):
        diffmat = mat_a / mat_b

    nan_mask = np.isnan(diffmat) | np.isinf(diffmat)
    nan_row_frac = np.mean(nan_mask, axis=1)
    valid_rows = nan_row_frac <= na_fraction
    valid_cols = nan_row_frac <= na_fraction

    if np.sum(valid_rows) < min_nbin:
        return {
            "chrom": chrom,
            "window_id": window_id,
            "window_start": win_start,
            "window_end": win_end,
            "tw_statistic": np.nan,
            "tw_pvalue": np.nan,
            "n_bins": n,
            "status": "INSUFFICIENT_VALID",
        }

    sub_diff = diffmat[np.ix_(valid_rows, valid_cols)]
    sub_diff[~np.isfinite(sub_diff)] = 1.0

    normed = normdiff_by_mean_sd(sub_diff)

    # Fill diagonal with 0 (normalized difference on diagonal should be ~0)
    np.fill_diagonal(normed, 0.0)

    tw_stat, p_value, status = compute_tracy_widom_stat(normed)

    return {
        "chrom": chrom,
        "window_id": window_id,
        "window_start": win_start,
        "window_end": win_end,
        "tw_statistic": tw_stat,
        "tw_pvalue": p_value,
        "n_bins": n,
        "status": status,
    }


def compare_iced_windows_tw(
    comparison, sample_a, sample_b, chroms, covariates_path,
    output, outdir, resolution, min_nbin, na_fraction, n_jobs, log_file
):
    """Compute Tracy-Widom test for ICE window pairs."""
    log(f"[INFO] Computing Tracy-Widom for {comparison}", log_file)
    log(f"[INFO] Sample A: {sample_a}, Sample B: {sample_b}", log_file)
    log(f"[INFO] Resolution: {resolution} bp", log_file)
    log(f"[INFO] min_nbin: {min_nbin}, na_fraction: {na_fraction}", log_file)

    try:
        cov_df = pd.read_csv(covariates_path, sep="\t")
        cov_df = cov_df[cov_df["chrom"].isin(chroms)]
        log(f"[INFO] Loaded {len(cov_df)} windows from covariates", log_file)

        log(f"[INFO] Processing {len(cov_df)} windows with {n_jobs} jobs...", log_file)

        results = Parallel(n_jobs=n_jobs, backend="loky", verbose=10)(
            delayed(process_single_window)(
                row=row,
                chrom=row["chrom"],
                sample_a=sample_a,
                sample_b=sample_b,
                resolution=resolution,
                outdir=outdir,
                min_nbin=min_nbin,
                na_fraction=na_fraction,
                log_file=log_file,
            )
            for _, row in cov_df.iterrows()
        )

        result_df = pd.DataFrame(results)
        result_df = result_df.sort_values(["chrom", "window_id"]).reset_index(drop=True)

        Path(output).parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(output, sep="\t", index=False)

        log(f"[INFO] Saved TW results to {output}", log_file)
        log(f"[INFO] Total windows: {len(result_df)}", log_file)

        valid_results = result_df.dropna(subset=["tw_statistic"])
        log(f"[INFO] Valid TW results: {len(valid_results)}/{len(result_df)}", log_file)

        if len(valid_results) > 0:
            sig_count = (valid_results["tw_pvalue"] < 0.05).sum()
            log(f"[INFO] Significant (p < 0.05): {sig_count}", log_file)

        return True

    except Exception as e:
        log(f"[ERROR] Failed to compute TW statistics: {e}", log_file)
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Compute Tracy-Widom statistics for ICE window pairs."
    )
    parser.add_argument("--comparison", required=True, help="Comparison name")
    parser.add_argument("--sample-a", required=True, help="Sample A name")
    parser.add_argument("--sample-b", required=True, help="Sample B name")
    parser.add_argument("--chroms", nargs="+", required=True, help="Chromosome names")
    parser.add_argument("--covariates", required=True, help="Covariates TSV file")
    parser.add_argument("--out", required=True, help="Output TSV file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--resolution", type=int, required=True, help="Hi-C resolution (bp)")
    parser.add_argument("--min-nbin", type=int, default=10, help="Minimum valid bins for TW test")
    parser.add_argument("--filter-threshold", type=float, default=0.5, help="Max NaN fraction per row")
    parser.add_argument("--n-jobs", type=int, default=80, help="Number of parallel jobs")
    parser.add_argument("--log", required=True, help="Log file path")
    args = parser.parse_args()

    Path(args.log).parent.mkdir(parents=True, exist_ok=True)

    success = compare_iced_windows_tw(
        comparison=args.comparison,
        sample_a=args.sample_a,
        sample_b=args.sample_b,
        chroms=args.chroms,
        covariates_path=args.covariates,
        output=args.out,
        outdir=args.outdir,
        resolution=args.resolution,
        min_nbin=args.min_nbin,
        na_fraction=args.filter_threshold,
        n_jobs=args.n_jobs,
        log_file=args.log,
    )

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
