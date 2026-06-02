#!/usr/bin/env python3
"""
Compute EC, IB, SP, and weighted features for all window pairs of a comparison.
Covariates are computed per-window pair (sample_a vs sample_b).

Usage:
    python scripts/compute_window_pair_covariates.py \
        --comparison H1398_rep1_vs_H1398_rep2 \
        --sample-a H1398_rep1 \
        --sample-b H1398_rep2 \
        --chroms HC04_A01 HC04_A02 \
        --out results/window_compare/H1398_rep1_vs_H1398_rep2/H1398_rep1_vs_H1398_rep2_covariates.tsv \
        --cpm-column raw_sum_cpm \
        --nz-column nz_frac \
        --pseudocount 0.01 \
        --weight-ec 0.33 \
        --weight-ib 0.33 \
        --weight-sp 0.33 \
        --n-jobs 80 \
        --log logs/window_compare/H1398_rep1_vs_H1398_rep2/compute_covariates.log
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed


def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def compute_covariate(window_row, sample_a_qc, sample_b_qc, cpm_col, nz_col, pseudocount):
    """Compute covariates for a single window pair."""
    sample = window_row["sample"]
    chrom = window_row["chrom"]
    window_id = window_row["window_id"]

    # Get QC metrics for both samples
    qc_a = sample_a_qc[
        (sample_a_qc["sample"] == sample.replace("_vs_", "_rep1_vs_").split("_rep1")[0] + "_rep1") &
        (sample_a_qc["chrom"] == chrom) &
        (sample_a_qc["window_id"] == window_id)
    ]
    qc_b = sample_b_qc[
        (sample_b_qc["sample"] == sample.replace("_vs_", "_rep2_vs_").split("_rep2")[0] + "_rep2") &
        (sample_b_qc["chrom"] == chrom) &
        (sample_b_qc["window_id"] == window_id)
    ]

    if qc_a.empty or qc_b.empty:
        return None

    cpm_a = max(float(qc_a[cpm_col].iloc[0]), 1e-10)
    cpm_b = max(float(qc_b[cpm_col].iloc[0]), 1e-10)
    nz_a = float(qc_a[nz_col].iloc[0])
    nz_b = float(qc_b[nz_col].iloc[0])

    # EC: Effective Coverage = harmonic mean of CPMs (matches SNAP-HiC-Diff_v2)
    ec = 2.0 * cpm_a * cpm_b / (cpm_a + cpm_b)

    # IB: Interaction Balance = |log2(CPM ratio)| (matches SNAP-HiC-Diff_v2)
    ib = abs(np.log2(cpm_a / cpm_b))

    # SP: SParity = mean of nz_fracs (matches SNAP-HiC-Diff_v2)
    sp = (nz_a + nz_b) / 2.0

    return {
        "chrom": chrom,
        "window_id": window_id,
        "window_start": window_row["window_start"],
        "window_end": window_row["window_end"],
        "EC": ec,
        "IB": ib,
        "SP": sp,
        "cpm_a": cpm_a,
        "cpm_b": cpm_b,
        "nz_a": nz_a,
        "nz_b": nz_b,
    }


def compute_window_pair_covariates(comparison, sample_a, sample_b, chroms, output,
                                   cpm_column, nz_column, pseudocount,
                                   weight_ec, weight_ib, weight_sp, outdir, n_jobs, log_file,
                                   qc_a_files=None, qc_b_files=None):
    """Compute covariates for all window pairs."""
    log(f"[INFO] Computing covariates for {comparison}", log_file)
    log(f"[INFO] Sample A: {sample_a}, Sample B: {sample_b}", log_file)
    log(f"[INFO] Chroms: {len(chroms)}", log_file)

    try:
        sample_a_qc = pd.DataFrame()
        sample_b_qc = pd.DataFrame()

        # Use explicitly provided file lists if available, otherwise construct from outdir
        if qc_a_files:
            for f in qc_a_files:
                if os.path.exists(f):
                    sample_a_qc = pd.concat([sample_a_qc, pd.read_csv(f, sep="\t")])
        else:
            for chrom in chroms:
                qc_a_chrom = f"{outdir}/windows/{sample_a}__{chrom}_window_qc.tsv"
                if os.path.exists(qc_a_chrom):
                    sample_a_qc = pd.concat([sample_a_qc, pd.read_csv(qc_a_chrom, sep="\t")])

        if qc_b_files:
            for f in qc_b_files:
                if os.path.exists(f):
                    sample_b_qc = pd.concat([sample_b_qc, pd.read_csv(f, sep="\t")])
        else:
            for chrom in chroms:
                qc_b_chrom = f"{outdir}/windows/{sample_b}__{chrom}_window_qc.tsv"
                if os.path.exists(qc_b_chrom):
                    sample_b_qc = pd.concat([sample_b_qc, pd.read_csv(qc_b_chrom, sep="\t")])

        if sample_a_qc.empty or sample_b_qc.empty:
            log(f"[WARNING] No QC data found, creating empty output", log_file)
            Path(output).parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(columns=["chrom", "window_id", "window_start", "window_end",
                                "EC", "IB", "SP"]).to_csv(output, sep="\t", index=False)
            return True

        log(f"[INFO] Sample A QC: {len(sample_a_qc)} windows", log_file)
        log(f"[INFO] Sample B QC: {len(sample_b_qc)} windows", log_file)

        # Filter to PASS windows in both samples before merge
        sample_a_qc = sample_a_qc[sample_a_qc["status"] == "PASS"]
        sample_b_qc = sample_b_qc[sample_b_qc["status"] == "PASS"]

        log(f"[INFO] Sample A PASS windows: {len(sample_a_qc)}", log_file)
        log(f"[INFO] Sample B PASS windows: {len(sample_b_qc)}", log_file)

        # Merge QC data on window_id
        merged = sample_a_qc.merge(
            sample_b_qc,
            on=["chrom", "window_id", "window_start", "window_end"],
            suffixes=("_a", "_b")
        )

        log(f"[INFO] Merged PASS windows: {len(merged)}", log_file)

        # Compute covariates (vectorized)
        cov_df = pd.DataFrame({
            "chrom": merged["chrom"],
            "window_id": merged["window_id"],
            "window_start": merged["window_start"],
            "window_end": merged["window_end"],
        })

        cpm_a = np.maximum(merged[f"{cpm_column}_a"].values, 1e-10)
        cpm_b = np.maximum(merged[f"{cpm_column}_b"].values, 1e-10)
        nz_a = merged[f"{nz_column}_a"].values
        nz_b = merged[f"{nz_column}_b"].values

        # EC: Effective Coverage = harmonic mean of CPMs (matches SNAP-HiC-Diff_v2)
        cov_df["EC"] = 2.0 * cpm_a * cpm_b / (cpm_a + cpm_b)

        # IB: Interaction Balance = |log2(CPM ratio)| (matches SNAP-HiC-Diff_v2)
        cov_df["IB"] = np.abs(np.log2(cpm_a / cpm_b))

        # SP: SParity = mean of nz_fracs (matches SNAP-HiC-Diff_v2)
        cov_df["SP"] = (nz_a + nz_b) / 2.0

        # Add raw values
        cov_df[f"{cpm_column}_a"] = cpm_a
        cov_df[f"{cpm_column}_b"] = cpm_b
        cov_df[f"{nz_column}_a"] = nz_a
        cov_df[f"{nz_column}_b"] = nz_b

        result_df = cov_df
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(output, sep="\t", index=False)

        log(f"[INFO] Saved covariates to {output}", log_file)
        log(f"[INFO] Total pairs: {len(result_df)}", log_file)
        return True

    except Exception as e:
        log(f"[ERROR] Failed to compute covariates: {e}", log_file)
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(description="Compute window pair covariates.")
    parser.add_argument("--comparison", required=True, help="Comparison name")
    parser.add_argument("--sample-a", required=True, help="Sample A name")
    parser.add_argument("--sample-b", required=True, help="Sample B name")
    parser.add_argument("--chroms", nargs="+", required=True, help="Chromosome names")
    parser.add_argument("--qc-a-files", nargs="+", default=None,
                        help="Explicit list of QC files for sample A (from Snakemake input)")
    parser.add_argument("--qc-b-files", nargs="+", default=None,
                        help="Explicit list of QC files for sample B (from Snakemake input")
    parser.add_argument("--out", required=True, help="Output TSV file")
    parser.add_argument("--cpm-column", default="raw_sum_cpm", help="CPM column name")
    parser.add_argument("--nz-column", default="nz_frac", help="NZ fraction column name")
    parser.add_argument("--pseudocount", type=float, default=0.01, help="Pseudocount for log transforms")
    parser.add_argument("--weight-ec", type=float, default=0.33, help="Weight for EC")
    parser.add_argument("--weight-ib", type=float, default=0.33, help="Weight for IB")
    parser.add_argument("--weight-sp", type=float, default=0.33, help="Weight for SP")
    parser.add_argument("--n-jobs", type=int, default=80, help="Number of parallel jobs")
    parser.add_argument("--log", required=True, help="Log file path")
    args = parser.parse_args()

    # Get outdir from the output path
    outdir = os.path.dirname(os.path.dirname(os.path.dirname(args.out)))

    Path(args.log).parent.mkdir(parents=True, exist_ok=True)

    success = compute_window_pair_covariates(
        comparison=args.comparison,
        sample_a=args.sample_a,
        sample_b=args.sample_b,
        chroms=args.chroms,
        output=args.out,
        cpm_column=args.cpm_column,
        nz_column=args.nz_column,
        pseudocount=args.pseudocount,
        weight_ec=args.weight_ec,
        weight_ib=args.weight_ib,
        weight_sp=args.weight_sp,
        outdir=outdir,
        n_jobs=args.n_jobs,
        log_file=args.log,
        qc_a_files=args.qc_a_files,
        qc_b_files=args.qc_b_files
    )

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
