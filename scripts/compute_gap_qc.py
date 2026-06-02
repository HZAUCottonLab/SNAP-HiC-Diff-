#!/usr/bin/env python3
"""
Compute gap detection and per-window QC for a single (sample, chrom) pair.
- Gap detection: balanced coverage (ICE) gradient analysis
- Per-window QC: raw_sum_cpm, nz_frac, bad_bin_frac, internal breakpoint scores
- Output columns raw_sum_cpm and nz_frac use ICE matrix (upper-triangle mean / finite-nonzero fraction)
- Dynamic P20 threshold applied per chromosome (raw_sum_cpm + nz_frac)
- Internal breakpoint: near-diagonal break + cross-quadrant depletion on ln(O/E) matrix
- Outputs: gap.tsv (per-bin) and window_qc.tsv (per-window with 8 breakpoint columns)

Usage:
    python compute_gap_qc.py \
        --chroms HC04_A01 \
        --sample H1398_rep1 \
        --balanced-cool-files results/balanced/H1398_rep1__HC04_A01_iced.cool \
        --ln-oe-files results/oe/H1398_rep1__HC04_A01_oe.npz \
        --outdir results/windows/H1398_rep1 \
        --resolution 4000 \
        --window-size 300000 \
        --step-size 32000 \
        --gradient-threshold 10.0 \
        --min-low-ratio 0.1 \
        --dynamic-raw-sum-cpm-pct 20 \
        --no-dynamic \
        --max-bad-bin-frac 0.5 \
        --short-bp 50000 \
        --flank-bp 80000 \
        --broad-bp 150000 \
        --min-margin-bp 30000 \
        --diag-drop-threshold 1.386 \
        --cross-median-threshold -1.050 \
        --total-score-threshold 1.733 \
        --log results/logs/gap_qc/H1398_rep1__HC04_A01.log \
        --n-jobs 10
"""

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cooler
import numpy as np
import pandas as pd
import scipy.sparse as sp
from joblib import Parallel, delayed

from detect_breakpoint_by_nan_density import detect_breakpoint_by_nan_density, build_gaps_mask


def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def find_breakpoint_by_gradient(sorted_signal, min_gradient_ratio):
    """Find breakpoint in sorted signal by gradient discontinuity."""
    if len(sorted_signal) < 2:
        return 0, 0.0, False
    ratios = []
    for i in range(len(sorted_signal) - 1):
        v_high = sorted_signal[i]
        v_low = sorted_signal[i + 1]
        if v_low <= 0:
            r = float("inf") if v_high > 0 else 1.0
        else:
            r = v_high / v_low
        ratios.append(r)
    ratios = np.array(ratios)
    max_ratio_idx = int(np.argmax(ratios))
    max_ratio_val = float(ratios[max_ratio_idx])
    breakpoint_idx = max_ratio_idx + 1
    is_valid = max_ratio_val >= min_gradient_ratio
    return breakpoint_idx, max_ratio_val, is_valid


def detect_gaps_in_window(bin_starts, bin_ends, cov_bal, win_start, win_end,
                           gradient_threshold, min_low_ratio):
    """Detect gap bins in a window using balanced coverage (ICE) gradient analysis."""
    mask = (bin_starts >= win_start) & (bin_ends <= win_end)
    idx = np.where(mask)[0]
    if len(idx) < 10:
        return np.array([], dtype=int)

    starts_win = bin_starts[idx]
    cov_win = cov_bal[idx]

    signal = np.nan_to_num(cov_win, nan=0.0)
    if len(signal) < 10 or np.all(signal == signal[0]) or np.all(signal == 0):
        return np.array([], dtype=int)

    sorted_indices = np.argsort(signal)[::-1]
    sorted_signal = signal[sorted_indices]

    bkp_idx, grad_val, is_grad_valid = find_breakpoint_by_gradient(
        sorted_signal, gradient_threshold
    )
    if not is_grad_valid:
        return np.array([], dtype=int)

    low_bin_ratio = (len(sorted_signal) - bkp_idx) / len(sorted_signal)
    if low_bin_ratio < min_low_ratio:
        return np.array([], dtype=int)

    low_idx = sorted_indices[bkp_idx:]
    return np.sort(idx[low_idx])


def process_chromosome(balanced_cool_path, outdir, chrom, sample,
                       resolution, window_size, step_size,
                       gradient_threshold, min_low_ratio,
                       dynamic_raw_sum_cpm_pct,
                       max_bad_bin_frac,
                       nan_abnormal_threshold,
                       log_file, n_jobs):
    """
    Process one chromosome: gap detection + per-window QC.
    Writes gap.tsv and window_qc.tsv.
    """
    log(f"[{sample}/{chrom}] Starting gap QC", log_file)
    log(f"  Balanced cool: {balanced_cool_path}", log_file)
    log(f"  Resolution: {resolution}, Window: {window_size}, Step: {step_size}", log_file)

    # ── Load balanced cool ────────────────────────────────────────────────────
    clr = cooler.Cooler(balanced_cool_path)
    mat_bal = clr.matrix(balance=True, sparse=True).fetch(chrom)
    ice_dense = mat_bal.toarray().astype(np.float64)
    cov_bal = np.nansum(ice_dense, axis=1).astype(np.float64)

    bins_df = clr.bins().fetch(chrom)
    bin_starts = bins_df["start"].to_numpy(dtype=np.int64)
    bin_ends = bins_df["end"].to_numpy(dtype=np.int64)
    n_bins = len(bin_starts)
    log(f"  Bins: {n_bins}, Balanced coverage range: [{cov_bal.min():.4f}, {cov_bal.max():.4f}]", log_file)

    # ── Build sliding windows ─────────────────────────────────────────────────
    chrom_length = int(bin_ends[-1])
    window_tasks = []
    window_id = 0
    for ws in range(0, chrom_length, step_size):
        we = min(ws + window_size, chrom_length)
        bin_mask = (bin_starts >= ws) & (bin_ends <= we)
        bin_idx = np.where(bin_mask)[0]
        if len(bin_idx) < (window_size // resolution) // 2:
            continue
        window_tasks.append({
            "window_id": window_id,
            "genomic_start": ws,
            "genomic_end": we,
            "bin_start": int(bin_idx[0]),
            "bin_end": int(bin_idx[-1]),
            "bin_count": len(bin_idx),
        })
        window_id += 1

    log(f"  Total windows: {len(window_tasks)}", log_file)

    if not window_tasks:
        log(f"  No windows to process for {chrom}", log_file)
        return True

    # ── NaN-density breakpoint params ─────────────────────────────────────────
    nan_density_params = {
        "resolution": resolution,
        "bin_starts": bin_starts,
        "bin_ends": bin_ends,
        "abnormal_threshold": nan_abnormal_threshold,
    }

    # ── Per-window processing (parallel) ──────────────────────────────────────
    def process_window(task, nan_density_params=nan_density_params):
        w = task
        s, e = w["bin_start"], w["bin_end"]
        n = e - s
        if n <= 0:
            return None

        # Extract submatrix (ICE — used for QC metrics)
        sub_ice = ice_dense[s:e, s:e].copy()

        # ── Gap detection ──────────────────────────────────────────────────────
        gap_bin_global = detect_gaps_in_window(
            bin_starts, bin_ends, cov_bal,
            w["genomic_start"], w["genomic_end"],
            gradient_threshold, min_low_ratio,
        )
        n_gaps = len(gap_bin_global)

        # ── Build gaps mask for NaN-density breakpoint detection ────────────────
        # Convert global gap bin indices to window-local indices
        gaps_mask = np.zeros((n, n), dtype=bool)
        for gbin in gap_bin_global:
            local_idx = gbin - s
            if 0 <= local_idx < n:
                gaps_mask[local_idx, :] = True
                gaps_mask[:, local_idx] = True

        # ── QC metrics (ICE upper triangle) ─────────────────────────────────────
        upper_mask = np.triu(np.ones((n, n), dtype=bool), k=1)
        upper_vals_ice = sub_ice[upper_mask]

        # raw_sum_cpm: ICE upper-triangle mean (always positive)
        raw_sum_cpm = float(np.nanmean(upper_vals_ice))

        # nz_frac: fraction of pixels that are both finite and non-zero
        total_pixels = upper_vals_ice.size
        nz_count = int(np.sum(np.isfinite(upper_vals_ice) & (upper_vals_ice != 0)))
        nz_frac = nz_count / max(total_pixels, 1)

        bad_bin_frac = n_gaps / max(n, 1)

        # ── NaN-density breakpoint detection (ICE matrix with gap mask) ──────────
        nan_bp_result = detect_breakpoint_by_nan_density(
            matrix=sub_ice,
            gaps_mask=gaps_mask,
            bin_size=nan_density_params["resolution"],
            window_start=w["genomic_start"],
            abnormal_threshold=nan_density_params["abnormal_threshold"],
        )

        # Status: FAIL logic
        status = "PASS"
        fail_reason = ""

        if bad_bin_frac > max_bad_bin_frac:
            status = "FAIL_BAD_BIN"
            fail_reason = "bad_bin_frac={:.3f}>{:.3f}".format(bad_bin_frac, max_bad_bin_frac)

        if nan_bp_result["is_abnormal"]:
            status = "FAIL_NAN_ABNORMAL"
            fail_reason = "nan_upper_zero_frac={:.3f}>{:.3f}".format(
                nan_bp_result["upper_zero_frac"], nan_density_params["abnormal_threshold"]
            )

        return {
            "sample": sample,
            "chrom": chrom,
            "window_id": w["window_id"],
            "window_start": w["genomic_start"],
            "window_end": w["genomic_end"],
            "bin_start": w["bin_start"],
            "bin_end": w["bin_end"],
            "n_gaps": n_gaps,
            "gap_bin_global": gap_bin_global,
            "raw_sum_cpm": raw_sum_cpm,
            "nz_frac": nz_frac,
            "bad_bin_frac": bad_bin_frac,
            # NaN-density breakpoint detection (ICE matrix with gap mask)
            "nan_is_abnormal": bool(nan_bp_result["is_abnormal"]),
            "nan_best_K": int(nan_bp_result["best_K"]),
            "nan_upper_zero_frac": float(nan_bp_result["upper_zero_frac"]) if not np.isnan(nan_bp_result["upper_zero_frac"]) else np.nan,
            "nan_status": nan_bp_result["status"],
            # QC thresholds
            "status": status,
            "fail_reason": fail_reason,
            "raw_sum_cpm_p20_threshold": None,
            "nz_frac_p20_threshold": None,
        }

    log(f"  Processing {len(window_tasks)} windows with {n_jobs} workers...", log_file)

    results = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(process_window)(task) for task in window_tasks
    )
    results = [r for r in results if r is not None]
    log(f"  Processed {len(results)} windows successfully", log_file)

    # ── Dynamic P20 threshold per chromosome (raw_sum_cpm + nz_frac) ───────────
    if dynamic_raw_sum_cpm_pct > 0:
        all_raw_sum = np.array([r["raw_sum_cpm"] for r in results], dtype=np.float64)
        all_nz_frac = np.array([r["nz_frac"] for r in results], dtype=np.float64)
        p20_raw_sum = float(np.percentile(all_raw_sum, dynamic_raw_sum_cpm_pct))
        p20_nz_frac = float(np.percentile(all_nz_frac, dynamic_raw_sum_cpm_pct))
        log(f"  Dynamic P20 threshold -- raw_sum_cpm: {p20_raw_sum:.6f}, nz_frac: {p20_nz_frac:.6f} (pct={dynamic_raw_sum_cpm_pct})", log_file)
        for r in results:
            r["raw_sum_cpm_p20_threshold"] = p20_raw_sum
            r["nz_frac_p20_threshold"] = p20_nz_frac
            if r["raw_sum_cpm"] < p20_raw_sum:
                r["status"] = "FAIL_RAW_SUM"
                r["fail_reason"] = "raw_sum_cpm={:.3f}<{:.3f}".format(r["raw_sum_cpm"], p20_raw_sum)
            if r["nz_frac"] < p20_nz_frac:
                r["status"] = "FAIL_NZ_FRAC"
                r["fail_reason"] = "nz_frac={:.3f}<{:.3f}".format(r["nz_frac"], p20_nz_frac)
    else:
        for r in results:
            r["raw_sum_cpm_p20_threshold"] = None
            r["nz_frac_p20_threshold"] = None

    # ── Write outputs ────────────────────────────────────────────────────────
    Path(outdir).mkdir(parents=True, exist_ok=True)
    prefix = f"{outdir}/{sample}__{chrom}"
    out_gaps = f"{prefix}_gaps.tsv"
    out_qc = f"{prefix}_window_qc.tsv"

    # gap.tsv
    all_gap_rows = []
    for r in results:
        for gbin in r["gap_bin_global"]:
            all_gap_rows.append({
                "sample": sample,
                "chrom": chrom,
                "bin_start": int(bin_starts[gbin]),
                "bin_end": int(bin_ends[gbin]),
                "window_id": r["window_id"],
                "window_start": r["window_start"],
                "window_end": r["window_end"],
            })

    if all_gap_rows:
        gap_df = pd.DataFrame(all_gap_rows).drop_duplicates(subset=["chrom", "bin_start", "bin_end", "window_id"])
        gap_df = gap_df.sort_values(["window_id", "bin_start"])
        gap_df.to_csv(out_gaps, index=False, sep="\t")
        log(f"  Wrote {out_gaps} ({len(gap_df)} gap bins from {len(results)} windows)", log_file)
    else:
        pd.DataFrame(columns=["sample", "chrom", "bin_start", "bin_end",
                              "window_id", "window_start", "window_end"]) \
            .to_csv(out_gaps, index=False, sep="\t")
        log(f"  Wrote {out_gaps} (no gaps detected)", log_file)

    # window_qc.tsv
    qc_df = pd.DataFrame([{k: v for k, v in r.items() if k != "gap_bin_global"} for r in results])
    pass_count = (qc_df["status"] == "PASS").sum()
    log(f"  Window QC: PASS={pass_count}/{len(qc_df)}, FAIL={len(qc_df)-pass_count}/{len(qc_df)}", log_file)

    qc_df.to_csv(out_qc, index=False, sep="\t")
    log(f"  Wrote {out_qc}", log_file)

    log(f"[{sample}/{chrom}] Done", log_file)
    return True


def main(chroms, sample, balanced_cool_files, outdir,
         resolution, window_size, step_size,
         gradient_threshold, min_low_ratio,
         dynamic_raw_sum_cpm_pct,
         max_bad_bin_frac,
         nan_abnormal_threshold,
         log_file, n_jobs):
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)

    for i, chrom in enumerate(chroms):
        ok = process_chromosome(
            balanced_cool_path=balanced_cool_files[i],
            outdir=outdir,
            chrom=chrom,
            sample=sample,
            resolution=resolution,
            window_size=window_size,
            step_size=step_size,
            gradient_threshold=gradient_threshold,
            min_low_ratio=min_low_ratio,
            dynamic_raw_sum_cpm_pct=dynamic_raw_sum_cpm_pct,
            max_bad_bin_frac=max_bad_bin_frac,
            nan_abnormal_threshold=nan_abnormal_threshold,
            log_file=log_file,
            n_jobs=n_jobs,
        )
        if not ok:
            print(f"[ERROR] Failed chromosome: {chrom}")
            return False

    print(f"[INFO] Finished sample {sample}: all {len(chroms)} chromosomes done")
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute gap detection and per-window QC."
    )
    parser.add_argument("--chroms", nargs="+", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--balanced-cool-files", nargs="+", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--resolution", type=int, required=True)
    parser.add_argument("--window-size", type=int, required=True)
    parser.add_argument("--step-size", type=int, required=True)
    parser.add_argument("--gradient-threshold", type=float, default=10.0)
    parser.add_argument("--min-low-ratio", type=float, default=0.04)
    parser.add_argument("--dynamic-raw-sum-cpm-pct", type=int, default=20,
                        help="Percentile for dynamic raw_sum_cpm + nz_frac threshold. 0 disables dynamic mode.")
    parser.add_argument("--no-dynamic", action="store_true",
                        help="Disable dynamic threshold, revert to simple PASS for raw_sum_cpm and nz_frac.")
    parser.add_argument("--max-bad-bin-frac", type=float, default=0.5)
    parser.add_argument("--nan-abnormal-threshold", type=float, default=0.50,
                        help="Upper band zero+NaN fraction threshold for is_abnormal (default: 0.50).")
    parser.add_argument("--log-file", required=True)
    parser.add_argument("--n-jobs", type=int, default=10)
    args = parser.parse_args()

    # Resolve dynamic pct: --no-dynamic overrides the flag value
    dynamic_pct = 0 if args.no_dynamic else args.dynamic_raw_sum_cpm_pct
    sys.exit(0 if main(
        chroms=args.chroms,
        sample=args.sample,
        balanced_cool_files=args.balanced_cool_files,
        outdir=args.outdir,
        resolution=args.resolution,
        window_size=args.window_size,
        step_size=args.step_size,
        gradient_threshold=args.gradient_threshold,
        min_low_ratio=args.min_low_ratio,
        dynamic_raw_sum_cpm_pct=dynamic_pct,
        max_bad_bin_frac=args.max_bad_bin_frac,
        nan_abnormal_threshold=args.nan_abnormal_threshold,
        log_file=args.log_file,
        n_jobs=args.n_jobs,
    ) else 1)
