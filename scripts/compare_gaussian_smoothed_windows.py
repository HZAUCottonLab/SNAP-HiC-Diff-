#!/usr/bin/env python3
"""
Compute MSE, SSIM, Spearman correlation for Gaussian-smoothed window pairs.
Compares smoothed Hi-C matrices between two samples.

Usage:
    python scripts/compare_gaussian_smoothed_windows.py \
        --comparison H1398_rep1_vs_H1398_rep2 \
        --sample-a H1398_rep1 \
        --sample-b H1398_rep2 \
        --chroms HC04_A01 HC04_A02 \
        --covariates results/window_compare/H1398_rep1_vs_H1398_rep2/H1398_rep1_vs_H1398_rep2_covariates.tsv \
        --out results/window_compare/H1398_rep1_vs_H1398_rep2/H1398_rep1_vs_H1398_rep2_gaussian_metrics.tsv \
        --n-jobs 80 \
        --log logs/window_compare/H1398_rep1_vs_H1398_rep2/compare_gaussian.log
"""

import argparse
import h5py
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from skimage.metrics import structural_similarity as ssim
from joblib import Parallel, delayed


def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def load_window_matrix(h5_path: str, window_id: int) -> np.ndarray:
    """Load a window matrix from gaussian.h5 file."""
    with h5py.File(h5_path, 'r') as f:
        key = f"w{window_id:06d}"
        matrix = f[key]['matrix'][:]
    return matrix


def compute_window_metrics(matrix_a: np.ndarray, matrix_b: np.ndarray) -> dict:
    """Compute MSE, SSIM, Spearman, NCC for two matrices."""
    # Ensure matrices are the same shape
    if matrix_a.shape != matrix_b.shape:
        min_size = min(matrix_a.shape[0], matrix_b.shape[0])
        matrix_a = matrix_a[:min_size, :min_size]
        matrix_b = matrix_b[:min_size, :min_size]

    # Use upper triangle (Hi-C is symmetric)
    mask = np.triu(np.ones_like(matrix_a, dtype=bool), k=0)

    vals_a = matrix_a[mask].flatten()
    vals_b = matrix_b[mask].flatten()

    # Remove NaN/Inf
    valid_mask = np.isfinite(vals_a) & np.isfinite(vals_b)
    vals_a = vals_a[valid_mask]
    vals_b = vals_b[valid_mask]

    if len(vals_a) < 10:
        return {"MSE": np.nan, "SSIM": np.nan, "Spearman": np.nan, "NCC": np.nan}

    # MSE (Mean Squared Error)
    mse = np.mean((vals_a - vals_b) ** 2)

    # SSIM (Structural Similarity Index) - use 2D image method
    # Normalize matrices to [0, 1] for SSIM
    min_val = min(np.nanmin(vals_a), np.nanmin(vals_b))
    max_val = max(np.nanmax(vals_a), np.nanmax(vals_b))
    range_val = max(max_val - min_val, 1e-10)

    norm_a = (vals_a - min_val) / range_val
    norm_b = (vals_b - min_val) / range_val

    # Reshape to 2D for SSIM
    size = int(np.sqrt(len(vals_a)))
    if size * size == len(vals_a):
        img_a = norm_a.reshape(size, size)
        img_b = norm_b.reshape(size, size)
        ssim_val = ssim(img_a, img_b, data_range=1.0)
    else:
        # Not a perfect square, use simplified SSIM-like metric
        # Based on mean, variance, and covariance
        mean_a = np.mean(norm_a)
        mean_b = np.mean(norm_b)
        var_a = np.var(norm_a)
        var_b = np.var(norm_b)
        cov_ab = np.mean((norm_a - mean_a) * (norm_b - mean_b))
        c1 = 0.01 ** 2
        c2 = 0.03 ** 2
        luminance = (2 * mean_a * mean_b + c1) / (mean_a**2 + mean_b**2 + c1)
        contrast = (2 * np.sqrt(var_a) * np.sqrt(var_b) + c2) / (var_a + var_b + c2)
        structure = (cov_ab + c2/2) / (np.sqrt(var_a) * np.sqrt(var_b) + c2/2)
        ssim_val = luminance * contrast * structure

    # Spearman correlation
    spearman_r, _ = spearmanr(vals_a, vals_b)

    # NCC (Normalized Cross-Correlation)
    mean_a = np.mean(vals_a)
    mean_b = np.mean(vals_b)
    std_a = np.std(vals_a)
    std_b = np.std(vals_b)

    if std_a > 0 and std_b > 0:
        ncc = np.mean((vals_a - mean_a) * (vals_b - mean_b)) / (std_a * std_b)
    else:
        ncc = np.nan

    return {
        "MSE": float(mse),
        "SSIM": float(ssim_val),
        "Spearman": float(spearman_r) if not np.isnan(spearman_r) else np.nan,
        "NCC": float(ncc),
    }


def process_single_window(
    row: pd.Series,
    chroms: list,
    sample_a: str,
    sample_b: str,
    windows_dir: str,
    h5_a_files: list = None,
    h5_b_files: list = None,
) -> dict:
    """Process a single window and return metrics."""
    chrom = row["chrom"]
    window_id = row["window_id"]

    h5_a = f"{windows_dir}/{sample_a}__{chrom}_gaussian.h5"
    h5_b = f"{windows_dir}/{sample_b}__{chrom}_gaussian.h5"

    try:
        if h5_a_files and h5_b_files:
            # Use explicit file lists provided by Snakemake
            h5_a = None
            h5_b = None
            for f in h5_a_files:
                if f"/{sample_a}__{chrom}_gaussian.h5" in f:
                    h5_a = f
                    break
            for f in h5_b_files:
                if f"/{sample_b}__{chrom}_gaussian.h5" in f:
                    h5_b = f
                    break
            if not h5_a or not h5_b:
                raise FileNotFoundError(f"H5 files not found for {chrom}")
        else:
            h5_a = f"{windows_dir}/{sample_a}__{chrom}_gaussian.h5"
            h5_b = f"{windows_dir}/{sample_b}__{chrom}_gaussian.h5"

        matrix_a = load_window_matrix(h5_a, window_id)
        matrix_b = load_window_matrix(h5_b, window_id)
        metrics = compute_window_metrics(matrix_a, matrix_b)
    except Exception as e:
        metrics = {"MSE": np.nan, "SSIM": np.nan, "Spearman": np.nan, "NCC": np.nan}

    return {
        "chrom": chrom,
        "window_id": window_id,
        "window_start": row["window_start"],
        "window_end": row["window_end"],
        **metrics,
    }


def compare_gaussian_smoothed_windows(comparison, sample_a, sample_b, chroms, covariates_path,
                                     output, outdir, n_jobs, log_file,
                                     h5_a_files=None, h5_b_files=None):
    """Compute Gaussian-smoothed comparison metrics for window pairs."""
    log(f"[INFO] Comparing Gaussian-smoothed windows for {comparison}", log_file)
    log(f"[INFO] Sample A: {sample_a}, Sample B: {sample_b}", log_file)
    log(f"[INFO] Chroms: {chroms}", log_file)

    try:
        # Load covariates
        cov_df = pd.read_csv(covariates_path, sep="\t")
        log(f"[INFO] Loaded {len(cov_df)} windows from covariates", log_file)

        # Filter to relevant chroms
        cov_df = cov_df[cov_df["chrom"].isin(chroms)]
        log(f"[INFO] Filtered to {len(cov_df)} windows for specified chroms", log_file)

        # Process windows in parallel
        # outdir = results/window_compare/{comparison}/...; .h5 lives in results/windows
        windows_dir = f"{outdir}/../windows"

        log(f"[INFO] Processing windows with {n_jobs} jobs...", log_file)

        results = Parallel(n_jobs=n_jobs, verbose=10)(
            delayed(process_single_window)(
                row=row,
                chroms=chroms,
                sample_a=sample_a,
                sample_b=sample_b,
                windows_dir=windows_dir,
                h5_a_files=h5_a_files,
                h5_b_files=h5_b_files,
            )
            for _, row in cov_df.iterrows()
        )

        result_df = pd.DataFrame(results)

        # Sort by chrom and window_id
        result_df = result_df.sort_values(["chrom", "window_id"]).reset_index(drop=True)

        Path(output).parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(output, sep="\t", index=False)

        log(f"[INFO] Saved Gaussian metrics to {output}", log_file)
        log(f"[INFO] Total pairs: {len(result_df)}", log_file)

        # Summary statistics
        valid = result_df.dropna(subset=["SSIM", "Spearman"])
        if len(valid) > 0:
            log(f"[INFO] Summary (valid windows: {len(valid)}/{len(result_df)}):", log_file)
            log(f"  MSE:     mean={valid['MSE'].mean():.4f}, std={valid['MSE'].std():.4f}", log_file)
            log(f"  SSIM:    mean={valid['SSIM'].mean():.4f}, std={valid['SSIM'].std():.4f}", log_file)
            log(f"  Spearman: mean={valid['Spearman'].mean():.4f}, std={valid['Spearman'].std():.4f}", log_file)
            log(f"  NCC:     mean={valid['NCC'].mean():.4f}, std={valid['NCC'].std():.4f}", log_file)

        return True

    except Exception as e:
        log(f"[ERROR] Failed to compare windows: {e}", log_file)
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(description="Compare Gaussian-smoothed windows.")
    parser.add_argument("--comparison", required=True, help="Comparison name")
    parser.add_argument("--sample-a", required=True, help="Sample A name")
    parser.add_argument("--sample-b", required=True, help="Sample B name")
    parser.add_argument("--chroms", nargs="+", required=True, help="Chromosome names")
    parser.add_argument("--covariates", required=True, help="Covariates TSV file")
    parser.add_argument("--h5-a-files", nargs="+", default=None,
                        help="Explicit list of H5 files for sample A (from Snakemake input)")
    parser.add_argument("--h5-b-files", nargs="+", default=None,
                        help="Explicit list of H5 files for sample B (from Snakemake input)")
    parser.add_argument("--out", required=True, help="Output TSV file")
    parser.add_argument("--n-jobs", type=int, default=80, help="Number of parallel jobs")
    parser.add_argument("--log", required=True, help="Log file path")
    args = parser.parse_args()

    outdir = os.path.dirname(os.path.dirname(args.out))

    Path(args.log).parent.mkdir(parents=True, exist_ok=True)

    success = compare_gaussian_smoothed_windows(
        comparison=args.comparison,
        sample_a=args.sample_a,
        sample_b=args.sample_b,
        chroms=args.chroms,
        covariates_path=args.covariates,
        output=args.out,
        outdir=outdir,
        n_jobs=args.n_jobs,
        log_file=args.log,
        h5_a_files=args.h5_a_files,
        h5_b_files=args.h5_b_files
    )

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
