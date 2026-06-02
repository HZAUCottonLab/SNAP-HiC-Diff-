#!/usr/bin/env python3
"""
Sliding window processing for Hi-C matrices: gap fill, z-score, and Gaussian smoothing.
Reads from ln(O/E) matrix and pre-computed gap.tsv, produces per-window HDF5 files.

Usage:
    python sliding_window.py \
        --gap-tsv-files results/windows/{sample}/{sample}__{chrom}_gaps.tsv \
        --ln-oe-files results/oe/{sample}__{chrom}_oe.npz \
        --outdir results/windows/{sample} \
        --chrom {chrom} \
        --sample {sample} \
        --resolution 4000 \
        --window-size 300000 \
        --step-size 32000 \
        --sigma 1.0 \
        --log results/logs/windows/{sample}__{chrom}.log
"""

import argparse
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from joblib import Parallel, delayed
from scipy.ndimage import gaussian_filter


def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def global_robust_zscore(matrix):
    """Per-window global robust z-score: (x - median) / MAD over all finite values.
    Matches SNAP-HiC-Diff_v2 robust_zscore() strategy.
    NaN positions (d≤1 diagonal, gap hard-mask) stay NaN.
    """
    valid = np.isfinite(matrix)
    if not np.any(valid):
        return np.full_like(matrix, np.nan)

    center = float(np.nanmedian(matrix[valid]))
    spread = float(np.nanmedian(np.abs(matrix[valid] - center)))
    spread = max(spread, 1e-10)
    zscore_mat = np.where(valid, (matrix - center) / spread, np.nan)
    return zscore_mat


# ── Stage 2 helpers ────────────────────────────────────────────────────────────

def _fill_gaps_2d_nearest(sub_oe, gap_bin_mask):
    """
    Fill 1-2 consecutive gap bins with 2D nearest interpolation.
    Hard-mask 3+ consecutive gap bins (fill with NaN).
    Skip diagonal zone (|i-j| <= 1) from interpolation.

    Parameters
    ----------
    sub_oe       : np.ndarray[n, n]  Ln(O/E) matrix, may contain NaN and gap bins.
    gap_bin_mask : np.ndarray[n]     **1D** boolean mask of gap bins for the current window.

        IMPORTANT — dimensionality semantics:
        gap_bin_mask[k] = True means: row k AND column k of the Hi-C matrix
        are gap regions. In 2D matrix terms, position (i, j) is a gap iff
        gap_bin_mask[i] == True OR gap_bin_mask[j] == True.
        This is NOT directly indexable as a 2D array — use _build_2d_gap_masks()
        to convert to 2D masks.

    Returns
    -------
    stage_filled : np.ndarray[n, n]
    """
    n = sub_oe.shape[0]

    # Build 2D masks from consecutive run analysis on gap_bin_mask
    fill_cand_2d, hard_mask_2d = _build_2d_gap_masks(gap_bin_mask, n)

    # Diagonal skip zone
    idx = np.arange(n)
    diag_skip_2d = np.abs(idx[:, None] - idx[None, :]) <= 1
    fill_cand_2d[diag_skip_2d] = False   # gap 区域的对角线不参与插值

    # NaN mask (non-gap, non-diag-skip)
    gap_2d = np.broadcast_to(gap_bin_mask, (n, n)) | np.broadcast_to(gap_bin_mask, (n, n)).T
    nan_mask = np.isnan(sub_oe) & ~gap_2d & ~diag_skip_2d

    # Build interpolation input: preserve real values, zero only positions that
    # should NOT contribute to interpolation.
    # We only zero the diagonal zone IF there are gaps to fill; otherwise
    # preserve backward compatibility (original NaN stays NaN → zeroed in z-score).
    has_gaps = gap_bin_mask.any()
    fill_input = sub_oe.copy()
    fill_input[np.isnan(fill_input)] = 0.0
    if has_gaps:
        fill_input[diag_skip_2d] = 0.0
    fill_input[hard_mask_2d] = 0.0

    # Interpolate gap fill candidates and non-gap NaN positions
    target_mask = np.zeros((n, n), dtype=bool)
    target_mask[fill_cand_2d] = True
    target_mask[nan_mask] = True

    filled = _interpolate_2d_nearest(fill_input, target_mask, fallback=0.0)
    filled[hard_mask_2d] = np.nan
    filled[diag_skip_2d & ~hard_mask_2d] = np.nan  # diagonal-adjacent → NaN, not fill

    return filled


def _build_2d_gap_masks(gap_bin_mask, n):
    """
    Analyze consecutive gap runs on the **1D** gap_bin_mask, then construct 2D masks.

    IMPORTANT — dimensionality semantics:
    gap_bin_mask[k] = True  →  row k AND column k are gap regions.
    Position (i, j) belongs to run [start, end) iff i∈[start,end) OR j∈[start,end).

    Returns two n×n boolean masks:
        fill_cand_2d : positions in a 1-2 consecutive gap run
        hard_mask_2d : positions in a 3+ consecutive gap run
    """
    fill_cand_2d = np.zeros((n, n), dtype=bool)
    hard_mask_2d = np.zeros((n, n), dtype=bool)

    # Find consecutive runs in gap_bin_mask
    runs = []  # list of (start, end, length)
    j = 0
    while j < n:
        if gap_bin_mask[j]:
            start = j
            while j < n and gap_bin_mask[j]:
                j += 1
            runs.append((start, j, j - start))
        else:
            j += 1

    for start, end, run_len in runs:
        if run_len <= 2:
            target = fill_cand_2d
        else:
            target = hard_mask_2d

        # Row side: any column, rows in [start, end)
        target[start:end, :] = True
        # Column side: any row, columns in [start, end)
        target[:, start:end] = True

    return fill_cand_2d, hard_mask_2d


def _interpolate_2d_nearest(fill_input, target_mask, fallback=0.0):
    """
    Nearest-neighbor interpolation that finds the closest non-zero source for each target.

    For each target (yi, xi), searches all grid positions for the nearest one
    with a non-zero value (not in skip or hard_mask zones), breaking ties by
    preferring same-row > same-column > nearest grid distance.

    Parameters
    ----------
    fill_input   : np.ndarray[n, n]  Valid positions have real values; zeros are placeholder.
    target_mask  : bool[n, n]        Positions to interpolate.
    fallback    : float              Value when no source found.

    Returns
    -------
    result : np.ndarray[n, n]
    """
    n = fill_input.shape[0]
    yy, xx = np.where(target_mask)
    if len(yy) == 0:
        return fill_input.copy()

    # Build source mask: positions with valid (non-zero) values
    source_mask = (fill_input != 0.0)
    if not source_mask.any():
        result = fill_input.copy()
        result[target_mask] = fallback
        return result

    src_yy, src_xx = np.where(source_mask)

    result = fill_input.copy()
    for yi, xi in zip(yy, xx):
        best_val = fallback
        min_dist = np.inf

        for sy, sx in zip(src_yy, src_xx):
            # Same-row sources preferred (distance bias = 0.5)
            row_bias = 0.5 if sy == yi else 0.0
            # Same-column sources preferred (distance bias = 0.5)
            col_bias = 0.5 if sx == xi else 0.0
            dist = np.sqrt((yi - sy) ** 2 + (xi - sx) ** 2) - row_bias - col_bias

            if dist < min_dist:
                min_dist = dist
                best_val = fill_input[sy, sx]

        result[yi, xi] = best_val

    return result


def process_chromosome(gap_tsv_path, ln_oe_npz_path, outdir, chrom, sample,
                       resolution, window_size, step_size, sigma, log_file, n_jobs):
    """
    Process one chromosome: read pre-computed gap bins, build windows,
    apply gap fill → z-score → Gaussian smoothing. Writes 4 HDF5 files.
    """
    def _log(msg):
        print(msg)
        if log_file:
            with open(log_file, "a") as f:
                f.write(f"{msg}\n")

    _log(f"[{sample}/{chrom}] Starting sliding window processing")
    _log(f"  Gap TSV: {gap_tsv_path}")
    _log(f"  Ln(O/E): {ln_oe_npz_path}")
    _log(f"  Resolution: {resolution} bp, Window: {window_size} bp, Step: {step_size} bp")

    # ── Load gap bins ──────────────────────────────────────────────────────────
    gap_df = pd.read_csv(gap_tsv_path, sep="\t")
    gap_bins = set(gap_df["bin_start"].values) if len(gap_df) > 0 else set()
    _log(f"  Loaded {len(gap_bins)} gap bins from TSV")

    # ── Load ln(O/E) matrix ────────────────────────────────────────────────────
    sparse_oe = sp.load_npz(ln_oe_npz_path)
    ln_oe = sparse_oe.toarray().astype(np.float32)
    ln_oe = (ln_oe + ln_oe.T) / 2.0
    n_bins = ln_oe.shape[0]
    # Diagonal-adjacent: d ≤ 1 → NaN (no valid ln_oe signal, would explode z-score)
    diag_skip_1d = np.abs(np.arange(n_bins)[:, None] - np.arange(n_bins)[None, :]) <= 1
    ln_oe[diag_skip_1d] = np.nan
    _log(f"  Ln(O/E) shape: {ln_oe.shape}, nnz: {sparse_oe.nnz}")

    # ── Build sliding windows ──────────────────────────────────────────────────
    # Reconstruct bin boundaries from resolution
    bin_starts = np.arange(n_bins, dtype=np.int64) * resolution
    bin_ends = bin_starts + resolution

    chrom_length = int(bin_ends[-1])
    window_tasks = []
    window_id = 0
    bins_per_window = window_size // resolution

    for ws in range(0, chrom_length, step_size):
        we = min(ws + window_size, chrom_length)
        bin_mask = (bin_starts >= ws) & (bin_ends <= we)
        bin_idx = np.where(bin_mask)[0]
        if len(bin_idx) < bins_per_window // 2:
            continue
        window_tasks.append({
            "window_id": window_id,
            "genomic_start": ws,
            "genomic_end": we,
            "bin_start": int(bin_idx[0]),
            "bin_end": int(bin_idx[-1]),
        })
        window_id += 1

    _log(f"  Total windows: {len(window_tasks)}")

    if not window_tasks:
        _log(f"  No windows to process for {chrom}")
        return True

    # ── Per-window processing (parallel) ─────────────────────────────────────
    def process_window(task):
        w = task
        s, e = w["bin_start"], w["bin_end"]
        n = e - s
        if n <= 0:
            return None

        sub_oe = ln_oe[s:e, s:e].copy().astype(np.float32)

        # Gap bin mask
        gap_bin_mask = np.zeros(n, dtype=bool)
        for gbin_start in gap_bins:
            gbin = gbin_start // resolution
            local = gbin - s
            if 0 <= local < n:
                gap_bin_mask[local] = True

        # Stage 1: ln(O/E) — preserve NaN
        stage_ln_oe = sub_oe.copy()

        # Stage 2: gap fill — 2D nearest interpolation
        stage_filled = _fill_gaps_2d_nearest(sub_oe, gap_bin_mask)

        # Stage 3: z-score
        stage_zscore = global_robust_zscore(stage_filled)

        # Stage 4: Gaussian smooth
        # gaussian_filter does not support NaN, so:
        # 1. Save NaN positions (3+ gap regions)
        # 2. Fill NaN with 0 for smoothing
        # 3. Apply gaussian
        # 4. Restore NaN positions
        nan_mask = np.isnan(stage_zscore)
        stage_gaussian_temp = stage_zscore.copy()
        stage_gaussian_temp[nan_mask] = 0.0
        stage_gaussian = gaussian_filter(stage_gaussian_temp, sigma=sigma, mode='constant', cval=0.0)
        stage_gaussian[nan_mask] = np.nan

        return {
            "window_id": w["window_id"],
            "genomic_start": w["genomic_start"],
            "genomic_end": w["genomic_end"],
            "stage_ln_oe": stage_ln_oe,
            "stage_filled": stage_filled,
            "stage_zscore": stage_zscore,
            "stage_gaussian": stage_gaussian,
        }

    _log(f"  Processing {len(window_tasks)} windows with {n_jobs} workers...")

    results = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(process_window)(task) for task in window_tasks
    )
    results = [r for r in results if r is not None]
    _log(f"  Processed {len(results)} windows successfully")

    # ── Write HDF5 files ──────────────────────────────────────────────────────
    prefix = f"{outdir}/{sample}__{chrom}"
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)

    def write_h5(path, stage_name, results):
        with h5py.File(path, "w") as f:
            for r in results:
                wid = f"w{r['window_id']:06d}"
                f.create_dataset(f"{wid}/matrix", data=r[f"stage_{stage_name}"], compression="gzip")
                f.create_dataset(f"{wid}/window_start", data=r["genomic_start"])
                f.create_dataset(f"{wid}/window_end", data=r["genomic_end"])

    _log(f"  Writing HDF5 files...")
    write_h5(f"{prefix}_ln_oe.h5", "ln_oe", results)
    write_h5(f"{prefix}_gap_filled.h5", "filled", results)
    write_h5(f"{prefix}_zscore.h5", "zscore", results)
    write_h5(f"{prefix}_gaussian.h5", "gaussian", results)
    _log(f"  Wrote: ln_oe.h5, gap_filled.h5, zscore.h5, gaussian.h5")

    _log(f"[{sample}/{chrom}] Done")
    return True


def main(chroms, sample, gap_tsv_files, ln_oe_files, outdir,
         resolution, window_size, step_size, sigma, log_file, n_jobs):
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)

    results = []
    for i, chrom in enumerate(chroms):
        r = process_chromosome(
            gap_tsv_path=gap_tsv_files[i],
            ln_oe_npz_path=ln_oe_files[i],
            outdir=outdir,
            chrom=chrom,
            sample=sample,
            resolution=resolution,
            window_size=window_size,
            step_size=step_size,
            sigma=sigma,
            log_file=log_file,
            n_jobs=n_jobs,
        )
        results.append((chrom, r))

    failures = [chrom for chrom, ok in results if not ok]
    if failures:
        print(f"[ERROR] Failed chromosomes: {failures}")
        return False
    print(f"[INFO] Finished sample {sample}: all {len(chroms)} chromosomes done")
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sliding window: gap fill, zscore, gaussian."
    )
    parser.add_argument("--chroms", nargs="+", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--gap-tsv-files", nargs="+", required=True)
    parser.add_argument("--ln-oe-files", nargs="+", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--resolution", type=int, required=True)
    parser.add_argument("--window-size", type=int, required=True)
    parser.add_argument("--step-size", type=int, required=True)
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("--n-jobs", type=int, default=10)
    args = parser.parse_args()
    sys.exit(0 if main(**vars(args)) else 1)
