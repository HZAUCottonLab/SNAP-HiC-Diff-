#!/usr/bin/env python3
"""
Detect internal breakpoint in a single Hi-C window using 5x5 sliding window NaN density.

Algorithm:
  1. Slide a 5x5 sub-window along the diagonal (step=1 bin).
  2. At each position k (top-left at k,k), count NaN pixels in the 5x5 window.
  3. The 5x5 window with the HIGHEST NaN fraction → center bin K = k+2 is the breakpoint.
  4. With breakpoint K, partition the matrix:
     - upper_band: X[i, j] where i > K, j < K  (right of K in x, above K in y)
     - lower_band: X[i, j] where i > K, j > K  (right of K in x, below K in y)
  5. If upper_band NaN fraction > 50% → is_abnormal=True.

Usage:
    python scripts/detect_breakpoint_by_nan_density.py \
        --cool results/balanced/H1398_rep1__HC04_A01_iced.cool \
        --chrom HC04_A01 \
        --window-start 31200000 \
        --window-end 34400000 \
        --bin-size 4000 \
        --out results/test_breakpoint.tsv

    # With extra gap masking from gaps.tsv:
    python scripts/detect_breakpoint_by_nan_density.py \
        --cool results/balanced/H1398_rep1__HC04_A01_iced.cool \
        --gaps results/windows/H1398_rep1/H1398_rep1__HC04_A01_gaps.tsv \
        --chrom HC04_A01 \
        --window-start 31200000 \
        --window-end 34400000 \
        --bin-size 4000 \
        --out results/test_breakpoint.tsv
"""

import argparse
import os
import sys
from pathlib import Path

import cooler
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Core pure functions
# ---------------------------------------------------------------------------

def _scan_5x5_zero_density(X, gaps_mask=None):
    """
    Slide 5x5 window along the diagonal, find position with highest zero+NaN fraction.

    For a window of size n×n:
      - Valid k range: 0 to n-5 (top-left positions where 5x5 fits)
      - For each k: X[k:k+5, k:k+5] → count zero or NaN → zero_frac = count / valid_cells
      - best_k = argmax(zero_frac) → K = best_k + 2 (center of 5x5)

    Args:
        X: np.ndarray, n×n matrix
        gaps_mask: np.ndarray[bool], shape (n,n), True = gap row/col to exclude from count

    Returns:
        (best_k, zero_fracs) where:
          best_k: int, top-left position of best window
          zero_fracs: np.ndarray, shape (n-4,), zero+NaN fraction for each k
    """
    n = X.shape[0]
    if n < 5:
        raise ValueError(f"Matrix too small for 5x5 scan: n={n}")

    zero_fracs = np.zeros(n - 4)
    for k in range(n - 4):
        window = X[k:k + 5, k:k + 5]

        # Count both zeros and NaN as "missing/depleted"
        zero_mask = (window == 0) | np.isnan(window)

        # Exclude gap rows/cols if mask provided
        if gaps_mask is not None:
            gap_window = gaps_mask[k:k + 5, k:k + 5]
            # Count valid (non-gap) cells
            valid_cells = ~gap_window
            zero_and_valid = zero_mask & valid_cells
            valid_count = valid_cells.sum()
            if valid_count == 0:
                zero_fracs[k] = np.nan
            else:
                zero_fracs[k] = zero_and_valid.sum() / valid_count
        else:
            zero_fracs[k] = zero_mask.sum() / 25.0

    best_k = int(np.nanargmax(zero_fracs)) if not np.all(np.isnan(zero_fracs)) else -1
    return best_k, zero_fracs


def _compute_upper_zero_frac(X, K, gaps_mask=None):
    """
    Compute zero+NaN fraction in the upper_band region (i > K, j < K).

    For n×n matrix with breakpoint K:
      - Row indices: K+1 to n-1  (n - K - 1 rows)
      - Col indices: 0 to K-1    (K cols)
      - upper_band shape: (n-K-1, K)
      - upper_zero_frac = count(zero or NaN, excluding gaps) / valid_count

    Args:
        X: np.ndarray, n×n matrix
        K: int, breakpoint bin index
        gaps_mask: np.ndarray[bool], shape (n,n), True = gap row/col to exclude

    Returns:
        float, zero+NaN fraction in upper_band. np.nan if upper_band is empty.
    """
    n = X.shape[0]
    upper = X[K + 1:, :K]  # i > K, j < K
    if upper.size == 0:
        return np.nan

    # Count both zeros and NaN as "missing/depleted"
    zero_mask = (upper == 0) | np.isnan(upper)

    if gaps_mask is not None:
        gap_upper = gaps_mask[K + 1:, :K]
        valid_mask = ~gap_upper
        valid_count = valid_mask.sum()
        if valid_count == 0:
            return np.nan
        zero_frac = (zero_mask & valid_mask).sum() / valid_count
    else:
        zero_frac = zero_mask.sum() / upper.size

    return float(zero_frac)


def detect_breakpoint_by_nan_density(matrix, gaps_mask=None, bin_size=4000, window_start=0, abnormal_threshold=0.50):
    """
    Detect breakpoint using 5x5 sliding window zero+NaN density.

    Args:
        matrix: np.ndarray, square n×n matrix (already extracted from cool file)
        gaps_mask: np.ndarray[bool], shape (n,n), True=gap row/col to exclude from calculation
        bin_size: int, resolution in bp
        window_start: int, genomic start position in bp
        abnormal_threshold: float, upper band zero+NaN fraction threshold (> = abnormal)

    Returns dict:
        - best_K: int, bin index of breakpoint center (-1 if skip)
        - best_K_bp: int, genomic position of breakpoint (-1 if skip)
        - is_abnormal: bool
        - upper_zero_frac: float, zero+NaN fraction in upper_band (gap regions excluded)
        - status: str, "ok" or "skip"
    """
    M = np.asarray(matrix, dtype=float)

    n = M.shape[0]

    # Skip windows too small for 5x5
    if n < 5:
        return {
            "best_K": -1,
            "best_K_bp": -1,
            "is_abnormal": False,
            "upper_zero_frac": np.nan,
            "status": "skip",
        }

    # NO symmetrization - detect raw zeros/NaN (breakpoints create asymmetric depletion)
    # gaps_mask passed directly to exclude gap regions from counting
    X = M

    # 5x5 zero+NaN density scan (gap regions excluded from count)
    best_k, zero_fracs = _scan_5x5_zero_density(X, gaps_mask=gaps_mask)

    # All-NaN case: no valid 5x5 windows
    if best_k < 0:
        return {
            "best_K": -1,
            "best_K_bp": -1,
            "is_abnormal": False,
            "upper_zero_frac": np.nan,
            "status": "skip",
        }

    K = best_k + 2  # center of 5x5

    # Upper band zero+NaN fraction (gap regions excluded)
    upper_zero_frac = _compute_upper_zero_frac(X, K, gaps_mask=gaps_mask)

    # Abnormal detection: upper band > threshold zero+NaN
    is_abnormal = bool(upper_zero_frac > abnormal_threshold)

    # Genomic position of breakpoint
    best_K_bp = int(window_start + K * bin_size)

    return {
        "best_K": K,
        "best_K_bp": best_K_bp,
        "is_abnormal": is_abnormal,
        "upper_zero_frac": upper_zero_frac,
        "status": "ok",
    }


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def build_gaps_mask(gaps_path, window_start, window_end, n, bin_size=4000):
    """
    Build n×n boolean mask from gaps.tsv.

    gaps.tsv columns: bin_start, bin_end (genomic coordinates in bp)
    The window spans [window_start, window_end) in genomic coords.
    Map each global gap bin to window-local index: local_idx = (gap_bin - window_start) // bin_size

    For each gap bin in the window range:
        mask[local_idx, :] = True; mask[:, local_idx] = True

    Returns:
        np.ndarray[bool], shape (n, n), True = gap/NaN. All False if gaps.tsv absent.
    """
    if not os.path.exists(gaps_path):
        return np.zeros((n, n), dtype=bool)

    df = pd.read_csv(gaps_path, sep="\t", usecols=["bin_start", "bin_end"])

    mask = np.zeros((n, n), dtype=bool)

    for _, row in df.iterrows():
        gap_start = int(row["bin_start"])
        gap_end = int(row["bin_end"])

        # Intersect gap range with window range
        intersect_start = max(gap_start, window_start)
        intersect_end = min(gap_end, window_end)
        if intersect_start >= intersect_end:
            continue

        # Convert intersecting genomic bins to window-local indices
        first_local = (intersect_start - window_start) // bin_size
        last_local = (intersect_end - 1 - window_start) // bin_size  # inclusive
        first_local = max(0, first_local)
        last_local = min(n - 1, last_local)

        for local_idx in range(first_local, last_local + 1):
            mask[local_idx, :] = True
            mask[:, local_idx] = True

    return mask


def fetch_window_matrix(cool_path, chrom, window_start, window_end):
    """
    Fetch window sub-matrix from a .cool file using sparse region fetching.

    Uses cooler's region fetch to get only the window region, avoiding
    loading the entire chromosome into memory.

    Returns:
        (matrix, bin_size) where matrix is np.ndarray (n×n)
    """
    clr = cooler.Cooler(cool_path)

    # Infer bin size from first bin
    bins_df = clr.bins().fetch(chrom)
    bin_size = int(bins_df["end"].iloc[0] - bins_df["start"].iloc[0])

    # Use cooler region fetch - returns sparse matrix for the window region
    region = f"{chrom}:{window_start}-{window_end}"
    mat_sparse = clr.matrix(balance=True, sparse=True).fetch(region)

    # Convert to dense
    mat = mat_sparse.toarray().astype(np.float64)

    return mat, bin_size


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Detect breakpoint by 5x5 NaN density in Hi-C window."
    )
    parser.add_argument("--cool", required=True, help="Path to iced.cool file.")
    parser.add_argument("--chrom", required=True, help="Chromosome name (e.g. HC04_A01).")
    parser.add_argument("--window-start", required=True, type=int,
                        help="Genomic start position (bp).")
    parser.add_argument("--window-end", required=True, type=int,
                        help="Genomic end position (bp).")
    parser.add_argument("--bin-size", type=int, default=4000,
                        help="Resolution in bp (default: 4000).")
    parser.add_argument("--gaps", default=None,
                        help="Path to gaps.tsv for extra masking.")
    parser.add_argument("--out", required=True, help="Output TSV path.")
    return parser.parse_args()


def main():
    args = parse_args()

    # Fetch window matrix from cool file
    mat, inferred_bin_size = fetch_window_matrix(
        args.cool, args.chrom, args.window_start, args.window_end
    )
    bin_size = args.bin_size if args.bin_size else inferred_bin_size
    n = mat.shape[0]

    # Compute window bin range
    window_bin_start = args.window_start // bin_size
    window_bin_end = args.window_end // bin_size

    # Build gaps mask if provided
    gaps_mask = None
    if args.gaps:
        gaps_mask = build_gaps_mask(args.gaps, args.window_start, args.window_end, n, bin_size)

    # Run detection
    result = detect_breakpoint_by_nan_density(
        matrix=mat,
        gaps_mask=gaps_mask,
        bin_size=bin_size,
        window_start=args.window_start,
    )

    # Write TSV output
    fields = ["window_id", "best_K", "best_K_bp", "is_abnormal", "upper_zero_frac", "status"]
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # window_id derived from output filename if possible
    window_id = out_path.stem.replace("_nan_density", "").replace("_breakpoint", "")

    with open(out_path, "w") as fout:
        fout.write("\t".join(fields) + "\n")
        vals = []
        for f in fields:
            if f == "window_id":
                vals.append(window_id)
            else:
                v = result.get(f, "")
                if isinstance(v, float) and np.isnan(v):
                    vals.append("nan")
                else:
                    vals.append(str(v))
        fout.write("\t".join(vals) + "\n")

    sys.stderr.write(
        "[DONE] best_K={} best_K_bp={} is_abnormal={} "
        "upper_zero_frac={:.4f} status={}\n".format(
            result["best_K"],
            result["best_K_bp"],
            result["is_abnormal"],
            result["upper_zero_frac"],
            result["status"],
        )
    )


if __name__ == "__main__":
    main()
