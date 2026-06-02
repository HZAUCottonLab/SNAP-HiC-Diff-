#!/usr/bin/env python3
"""
Compute O/E matrix from a single balanced per-chromosome .cool.
Accepts parameters via argparse to avoid pickle round-trip issues with
Snakemake's dynamically created NamedList subclasses.

Usage:
    python compute_oe_matrix.py \
        --input-balanced /path/L128_rep1__HC04_A01_iced.cool \
        --input-expected /path/L128_rep1__HC04_A01_expected.npz \
        --output /path/L128_rep1__HC04_A01_oe.npz
"""

import argparse
import os
import sys
import numpy as np
import scipy.sparse as sp
from pathlib import Path
import cooler


def parse_chrom_from_balanced_path(input_path):
    basename = os.path.basename(input_path)
    parts = basename.replace("_iced.cool", "").split("__")
    if len(parts) == 2:
        return parts[0], parts[1]
    return None, None


def load_expected(expected_path):
    data = np.load(expected_path)
    return data["expected"], data.get("counts", None)


def compute_oe_sparse(mat_bal, expected_vec):
    n = mat_bal.shape[0]
    max_dist = min(len(expected_vec), n)
    rows, cols, data = [], [], []
    mat_coo = mat_bal.tocoo()
    for i, j, val in zip(mat_coo.row, mat_coo.col, mat_coo.data):
        if i < j:
            dist = j - i
            if dist < max_dist and expected_vec[dist] > 0:
                oe = val / expected_vec[dist]
                if oe > 0:
                    ln_oe = np.log(oe)
                    rows.extend([i, j])
                    cols.extend([j, i])
                    data.extend([ln_oe, ln_oe])
    if rows:
        return sp.csr_matrix((data, (rows, cols)), shape=(n, n), dtype=np.float32)
    return sp.csr_matrix((n, n), dtype=np.float32)


def main(input_balanced, input_expected, output):
    sample, chrom = parse_chrom_from_balanced_path(input_balanced)
    if not chrom:
        print(f"ERROR: Could not parse chromosome from {input_balanced}")
        return False

    print(f"Computing O/E matrix")
    print(f"  Sample: {sample}, Chromosome: {chrom}")
    print(f"  Input balanced: {input_balanced}")
    print(f"  Input expected: {input_expected}")
    print(f"  Output: {output}")

    Path(output).parent.mkdir(parents=True, exist_ok=True)

    expected_vec, _ = load_expected(input_expected)
    print(f"  Loaded expected values: {len(expected_vec)} distance bins")

    clr = cooler.Cooler(input_balanced)
    chroms_in_file = list(clr.chromnames)
    chrom = chroms_in_file[0] if chroms_in_file else chrom
    mat_bal = clr.matrix(balance=True, sparse=True).fetch(chrom).tocsr()
    print(f"  Matrix shape: {mat_bal.shape}, nnz: {mat_bal.nnz}")

    oe_sparse = compute_oe_sparse(mat_bal, expected_vec)
    print(f"  O/E matrix nnz: {oe_sparse.nnz}")

    total_upper = mat_bal.shape[0] * (mat_bal.shape[0] - 1) // 2
    nan_frac = (total_upper - oe_sparse.nnz) / max(total_upper, 1)
    print(f"  NaN fraction: {nan_frac:.3f}")

    sp.save_npz(output, oe_sparse)
    print(f"Saved to {output}")
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute O/E matrix from balanced .cool.")
    parser.add_argument("--input-balanced", required=True, help="Path to balanced .cool file")
    parser.add_argument("--input-expected", required=True, help="Path to expected .npz file")
    parser.add_argument("--output", required=True, help="Path to output O/E .npz file")
    args = parser.parse_args()
    sys.exit(0 if main(**vars(args)) else 1)
