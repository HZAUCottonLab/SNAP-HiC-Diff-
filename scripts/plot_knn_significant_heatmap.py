#!/usr/bin/env python3
"""
2x3 heatmap plotting for significant KNN windows (raw p-value based).

Scan mode writes three output files to --out-dir:
  1. significant_windows.tsv  - manifest (chrom/window_id/start/end)
  2. pval_table.tsv           - all windows, per-metric raw p-values + pass/fail
  3. pval_summary.tsv         - per-metric summary statistics

Usage:
    # Scan mode: discover significant windows, write manifest + pval summary/table
    python scripts/plot_knn_significant_heatmap.py \
        --scan-only \
        --enriched results/knn_empirical_pvalue/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv \
        --out-dir results/logs/knn_significant_heatmap/H1398_vs_L134/

    # Chrom mode: plot all significant windows on one chromosome (2x3: iced / ln(O/E) / gaussian)
    python scripts/plot_knn_significant_heatmap.py \
        --chrom-mode \
        --comparison H1398_vs_L134 \
        --chrom HC04_A01 \
        --enriched results/knn_empirical_pvalue/H1398_vs_L134/H1398_vs_L134_knn_enriched.tsv \
        --out-dir results/knn_significant_heatmap/H1398_vs_L134/
"""

import argparse
import os
import sys
import warnings
from pathlib import Path

import cooler
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter

os.environ["NUMEXPR_MAX_THREADS"] = "64"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
warnings.filterwarnings("ignore")

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

METRICS = [
    "MSE",
    "d_SSIM",
    "d_Spearman",
    "tw_score",
    "d_NCC",
]

PVALUE_COLS = [f"{m}_empirical_pvalue" for m in METRICS]


def filter_significant_windows(df, pvalue_threshold=0.05, pass_mode="all", metrics=None):
    """Filter significant windows using raw empirical p-values.

    pass_mode:
      - all:  all selected metrics pass pvalue_threshold
      - any:  any selected metric passes pvalue_threshold
      - subset: treat 'metrics' as subset; rule is 'all' over that subset
    """
    if metrics is None:
        metrics = list(METRICS)

    cols = [f"{m}_empirical_pvalue" for m in metrics]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"Missing p-value columns in enriched TSV: {missing}")

    pval_vals = df[cols].values.astype(float)
    pass_mat = pval_vals < pvalue_threshold

    mode = (pass_mode or "all").lower()
    if mode == "all" or mode == "subset":
        keep = np.all(pass_mat, axis=1)
    elif mode == "any":
        keep = np.any(pass_mat, axis=1)
    else:
        raise ValueError(f"Unknown pass_mode: {mode}")

    return df.loc[keep].copy()


def build_pval_df(df, pvalue_threshold=0.05, metrics=None):
    """Build per-window p-value table with pass/fail columns."""
    if metrics is None:
        metrics = list(METRICS)

    result = df[["chrom", "window_id", "window_start", "window_end", "comparison"]].copy()
    for m in metrics:
        col = f"{m}_empirical_pvalue"
        if col not in df.columns:
            raise KeyError(f"Missing p-value column: {col}")
        result[col] = df[col].values.astype(float)
        result[f"{m}_pass"] = (result[col] < pvalue_threshold).astype(int)
    return result


def build_summary_table(df, pvalue_threshold=0.05, metrics=None):
    """Build per-metric summary statistics from raw p-value columns."""
    if metrics is None:
        metrics = list(METRICS)

    rows = []
    for m in metrics:
        pcol = f"{m}_empirical_pvalue"
        if pcol not in df.columns:
            raise KeyError(f"Missing p-value column: {pcol}")

        raw_pvals = df[pcol].values.astype(float)
        pass_mask = raw_pvals < pvalue_threshold
        n_pass = int(np.sum(pass_mask))
        n_total = len(raw_pvals)
        pct_pass = n_pass / n_total * 100

        pass_pvals = raw_pvals[pass_mask]
        fail_pvals = raw_pvals[~pass_mask]

        rows.append({
            "metric": m,
            "n_total": n_total,
            "n_pass": n_pass,
            "pct_pass": round(pct_pass, 2),
            "pass_median_p": round(float(np.median(pass_pvals)), 6) if len(pass_pvals) > 0 else float("nan"),
            "fail_median_p": round(float(np.median(fail_pvals)), 6) if len(fail_pvals) > 0 else float("nan"),
            "pass_min_p": round(float(np.min(pass_pvals)), 8) if len(pass_pvals) > 0 else float("nan"),
            "pass_max_p": round(float(np.max(pass_pvals)), 6) if len(pass_pvals) > 0 else float("nan"),
        })

    summary_df = pd.DataFrame(rows)

    # Joint pass counts
    cols = [f"{m}_empirical_pvalue" for m in metrics]
    pval_matrix = df[cols].values.astype(float)
    pass_counts = np.sum(pval_matrix < pvalue_threshold, axis=1)
    summary_df["n_pass_0_metrics"] = int(np.sum(pass_counts == 0))
    summary_df["n_pass_1_metric"] = int(np.sum(pass_counts == 1))
    summary_df["n_pass_2_metrics"] = int(np.sum(pass_counts == 2))
    summary_df["n_pass_3_metrics"] = int(np.sum(pass_counts == 3))
    summary_df["n_pass_4_metrics"] = int(np.sum(pass_counts == 4))
    summary_df["n_pass_5_metrics"] = int(np.sum(pass_counts == 5))
    summary_df["n_pass_any"] = int(np.sum(pass_counts >= 1))
    summary_df["n_pass_all5"] = int(np.sum(pass_counts == len(metrics)))
    summary_df["pvalue_threshold"] = pvalue_threshold

    return summary_df


def get_pair_samples(comparison):
    """Split 'A_vs_B' -> (control='A', treatment='B')."""
    parts = comparison.split("_vs_")
    if len(parts) != 2:
        raise ValueError(f"Comparison '{comparison}' must be in format 'A_vs_B'")
    return parts[0], parts[1]


def load_iced_matrix(sample, chrom, window_start, window_end):
    """Load ICE O/E matrix from per-chrom .cool file."""
    cool_path = os.path.join(
        PROJECT_ROOT, "results", "balanced", sample,
        f"{sample}__{chrom}_iced.cool"
    )
    if not os.path.exists(cool_path):
        raise FileNotFoundError(f"Cool file not found: {cool_path}")
    clr = cooler.Cooler(cool_path)
    chrom_names = list(clr.chromnames)
    actual_chrom = chrom if chrom in chrom_names else chrom_names[0]
    mat = (
        clr.matrix(balance=True, sparse=True)
        .fetch(f"{actual_chrom}:{window_start}-{window_end}")
        .toarray()
    )
    return (mat + mat.T) / 2.0


def load_gaussian_matrix(sample, chrom, window_id):
    """Load Gaussian-smoothed matrix from per-chrom HDF5 file."""
    h5_path = os.path.join(
        PROJECT_ROOT, "results", "windows", sample,
        f"{sample}__{chrom}_gaussian.h5"
    )
    if not os.path.exists(h5_path):
        raise FileNotFoundError(f"H5 file not found: {h5_path}")
    key = f"w{window_id:06d}"
    with h5py.File(h5_path, "r") as f:
        if key not in f:
            raise KeyError(f"Key '{key}' not found in {h5_path}")
        mat = f[key]["matrix"][:]
    return (mat + mat.T) / 2.0


def _quantile_vmin_vmax(mats, q_low=0.01, q_high=0.99):
    """Compute robust vmin/vmax from one or more matrices using quantiles.

    Ignores NaN/Inf. Returns (vmin, vmax). Falls back to min/max if
    quantile computation is not possible.
    """
    vals = []
    for m in mats:
        a = np.asarray(m, dtype=float).ravel()
        a = a[np.isfinite(a)]
        if a.size:
            vals.append(a)
    if not vals:
        return (0.0, 1.0)
    allv = np.concatenate(vals)
    if allv.size == 0:
        return (0.0, 1.0)
    try:
        vmin = float(np.quantile(allv, q_low))
        vmax = float(np.quantile(allv, q_high))
    except Exception:
        vmin = float(np.nanmin(allv))
        vmax = float(np.nanmax(allv))

    if not np.isfinite(vmin) or not np.isfinite(vmax):
        vmin = float(np.nanmin(allv))
        vmax = float(np.nanmax(allv))

    if vmin == vmax:
        eps = 1e-9 if vmin == 0 else abs(vmin) * 1e-6
        vmin -= eps
        vmax += eps

    if vmin > vmax:
        vmin, vmax = vmax, vmin

    return vmin, vmax


# ------------------------------------------------------------------------------
def run_scan(enriched_path, out_dir, pvalue_threshold=0.05, pass_mode="all", metrics=None):
    """
    Scan knn_enriched.tsv, filter significant windows using raw empirical p-values.
    Writes to --out-dir:
      1. significant_windows.tsv - manifest (significant windows only)
      2. pval_table.tsv          - all windows, per-metric p-values + pass/fail
      3. pval_summary.tsv        - per-metric summary statistics
    """
    df = pd.read_csv(enriched_path, sep="\t")
    n_total = len(df)

    sig_df = filter_significant_windows(
        df,
        pvalue_threshold=pvalue_threshold,
        pass_mode=pass_mode,
        metrics=metrics,
    )
    n_sig = len(sig_df)

    comparison = os.path.basename(os.path.dirname(enriched_path)).replace(
        "_knn_empirical_pvalue.tsv", ""
    )

    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # 1. manifest TSV (significant windows only)
    manifest_df = sig_df[["chrom", "window_id", "window_start", "window_end"]].copy()
    manifest_df = manifest_df.sort_values(["chrom", "window_id"]).reset_index(drop=True)
    manifest_df.to_csv(
        os.path.join(out_dir, "significant_windows.tsv"), sep="\t", index=False
    )

    # 2. p-value full table (all windows)
    pval_df = build_pval_df(df, pvalue_threshold=pvalue_threshold, metrics=metrics)
    pval_df.sort_values(["chrom", "window_id"]).reset_index(drop=True).to_csv(
        os.path.join(out_dir, "pval_table.tsv"), sep="\t", index=False
    )

    # 3. Summary table (per-metric stats)
    summary_df = build_summary_table(df, pvalue_threshold=pvalue_threshold, metrics=metrics)
    summary_df.to_csv(
        os.path.join(out_dir, "pval_summary.tsv"), sep="\t", index=False
    )

    chroms = sorted(sig_df["chrom"].unique())

    print(
        f"[scan] {comparison}: {n_total} total, {n_sig} significant "
        f"(pass_mode={pass_mode}, threshold={pvalue_threshold}), {len(chroms)} chromosomes -> {out_dir}",
        flush=True,
    )


def mask_iced_diagonal(mat):
    """Set diagonal + adjacent bins to NaN."""
    m = mat.copy().astype(float)
    n = m.shape[0]
    for i in range(n):
        for j in range(max(0, i - 1), min(n, i + 2)):
            m[i, j] = np.nan
    return m


def load_ln_oe_matrix(sample, chrom, window_id):
    """Load ln(O/E) matrix from per-window HDF5."""
    h5_path = os.path.join(
        PROJECT_ROOT, "results", "windows", sample,
        f"{sample}__{chrom}_ln_oe.h5"
    )
    if not os.path.exists(h5_path):
        raise FileNotFoundError(f"ln_oe h5 not found: {h5_path}")
    key = f"w{window_id:06d}"
    with h5py.File(h5_path, "r") as f:
        if key not in f:
            raise KeyError(f"Key '{key}' not found in {h5_path}")
        mat = f[key]["matrix"][:]
    return (mat + mat.T) / 2.0


def plot_2x3_heatmap(
    m_ctrl_iced, m_treat_iced,
    m_ctrl_lnoe, m_treat_lnoe,
    m_ctrl_gauss, m_treat_gauss,
    window_start, window_end, chrom,
    control_name, treatment_name,
    out_png, out_pdf,
):
    """
    Draw a 2x3 heatmap.
    Row 0: control; Row 1: treatment.
    Col 0: iced (masked diagonal+adjacent); Col 1: ln(O/E); Col 2: gaussian.
    Each column has its own colorbar.
    """
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    vmin_iced, vmax_iced = _quantile_vmin_vmax([m_ctrl_iced, m_treat_iced], 0.01, 0.99)
    norm_iced = mcolors.Normalize(vmin=vmin_iced, vmax=vmax_iced)
    axes[0, 0].imshow(
        m_ctrl_iced, norm=norm_iced, cmap="YlOrRd",
        extent=[window_start, window_end, window_end, window_start],
        origin="upper", aspect="auto",
    )
    axes[0, 0].set_title(f"{control_name}\niced")
    axes[1, 0].imshow(
        m_treat_iced, norm=norm_iced, cmap="YlOrRd",
        extent=[window_start, window_end, window_end, window_start],
        origin="upper", aspect="auto",
    )
    axes[1, 0].set_title(f"{treatment_name}\niced")

    vmin_ln, vmax_ln = _quantile_vmin_vmax([m_ctrl_lnoe, m_treat_lnoe], 0.01, 0.99)
    norm_ln = mcolors.Normalize(vmin=vmin_ln, vmax=vmax_ln)
    axes[0, 1].imshow(
        m_ctrl_lnoe, norm=norm_ln, cmap="RdBu_r",
        extent=[window_start, window_end, window_end, window_start],
        origin="upper", aspect="auto",
    )
    axes[0, 1].set_title(f"{control_name}\nln(O/E)")
    axes[1, 1].imshow(
        m_treat_lnoe, norm=norm_ln, cmap="RdBu_r",
        extent=[window_start, window_end, window_end, window_start],
        origin="upper", aspect="auto",
    )
    axes[1, 1].set_title(f"{treatment_name}\nln(O/E)")

    vmin_gauss = min(np.nanmin(m_ctrl_gauss), np.nanmin(m_treat_gauss))
    vmax_gauss = max(np.nanmax(m_ctrl_gauss), np.nanmax(m_treat_gauss))
    norm_gauss = mcolors.Normalize(vmin=vmin_gauss, vmax=vmax_gauss)
    axes[0, 2].imshow(
        m_ctrl_gauss, norm=norm_gauss, cmap="RdBu_r",
        extent=[window_start, window_end, window_end, window_start],
        origin="upper", aspect="auto",
    )
    axes[0, 2].set_title(f"{control_name}\ngaussian")
    axes[1, 2].imshow(
        m_treat_gauss, norm=norm_gauss, cmap="RdBu_r",
        extent=[window_start, window_end, window_end, window_start],
        origin="upper", aspect="auto",
    )
    axes[1, 2].set_title(f"{treatment_name}\ngaussian")

    bp_fmt = FuncFormatter(lambda x, _: f"{int(x/1000):,}kb")
    for ax in axes.flat:
        ax.set_xlabel(f"{chrom} position (bp)")
        ax.set_ylabel(f"{chrom} position (bp)")
        ax.locator_params(axis="both", nbins=5)
        ax.xaxis.set_major_formatter(bp_fmt)
        ax.yaxis.set_major_formatter(bp_fmt)

    for col in range(3):
        im = axes[0, col].images[0]
        label = ["O/E (iced, masked)", "ln(O/E)", "z-score (gaussian)"][col]
        fig.colorbar(im, ax=axes[:, col], location="right",
                     fraction=0.02, pad=0.04, label=label)

    fig.suptitle(f"Chromosome {chrom}  window {window_start}-{window_end}", fontsize=14)
    fig.tight_layout()

    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def run_chrom(
    comparison,
    chrom,
    enriched_path,
    out_dir,
    manifest_path=None,
    pvalue_threshold=0.05,
    pass_mode="all",
    metrics=None,
):
    """Plot all significant windows on one chromosome."""
    control, treatment = get_pair_samples(comparison)

    if manifest_path is not None:
        chrom_df = pd.read_csv(manifest_path, sep="\t")
        chrom_df = chrom_df[chrom_df["chrom"] == chrom].sort_values("window_id")
    else:
        df = pd.read_csv(enriched_path, sep="\t")
        sig_df = filter_significant_windows(
            df,
            pvalue_threshold=pvalue_threshold,
            pass_mode=pass_mode,
            metrics=metrics,
        )
        chrom_df = sig_df[sig_df["chrom"] == chrom].sort_values("window_id")

    n_win = len(chrom_df)
    if n_win == 0:
        print(
            f"[chrom] {comparison} | {chrom}: no significant windows, skipping",
            flush=True,
        )
        return

    print(
        f"[chrom] {comparison} | {chrom}: {n_win} windows -> {out_dir}",
        flush=True,
    )

    os.makedirs(out_dir, exist_ok=True)

    for _, row in chrom_df.iterrows():
        wid = int(row["window_id"])
        ws = int(row["window_start"])
        we = int(row["window_end"])

        png_path = os.path.join(out_dir, f"{chrom}_{wid}.png")
        pdf_path = os.path.join(out_dir, f"{chrom}_{wid}.pdf")

        m_ci_raw = load_iced_matrix(control, chrom, ws, we)
        m_ti_raw = load_iced_matrix(treatment, chrom, ws, we)
        m_ci = mask_iced_diagonal(m_ci_raw)
        m_ti = mask_iced_diagonal(m_ti_raw)
        m_cl = load_ln_oe_matrix(control, chrom, wid)
        m_tl = load_ln_oe_matrix(treatment, chrom, wid)
        m_cg = load_gaussian_matrix(control, chrom, wid)
        m_tg = load_gaussian_matrix(treatment, chrom, wid)

        plot_2x3_heatmap(
            m_ci, m_ti, m_cl, m_tl, m_cg, m_tg,
            ws, we, chrom, control, treatment,
            png_path, pdf_path,
        )

    print(
        f"[chrom] {comparison} | {chrom}: done ({n_win} heatmaps)",
        flush=True,
    )


# ------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="2x3 heatmap for significant KNN windows (raw p-value based)."
    )
    parser.add_argument(
        "--scan-only",
        action="store_true",
        help="Scan mode: write manifest + p-value summary/table",
    )
    parser.add_argument(
        "--chrom-mode",
        action="store_true",
        help="Chrom mode: plot all significant windows on one chromosome",
    )
    parser.add_argument("--comparison", help="Comparison name (e.g. H1398_vs_L134)")
    parser.add_argument("--chrom", help="Chromosome identifier (e.g. HC04_A01)")
    parser.add_argument(
        "--manifest",
        help="Optional manifest TSV from scan mode (significant_windows.tsv)",
    )
    parser.add_argument(
        "--enriched",
        required=True,
        help="Path to knn_enriched.tsv (relative to project root or absolute)",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--pvalue-threshold",
        type=float,
        default=0.05,
        help="Raw empirical p-value threshold for calling significant windows",
    )
    parser.add_argument(
        "--pass-mode",
        choices=["all", "any", "subset"],
        default="all",
        help="How to combine per-metric pass/fail (default all)",
    )
    parser.add_argument(
        "--metrics",
        default=",".join(METRICS),
        help="Comma-separated metrics list (default all 5)",
    )
    args = parser.parse_args()

    metrics = [m.strip() for m in (args.metrics or "").split(",") if m.strip()]
    if len(metrics) == 0:
        metrics = list(METRICS)

    # Resolve enriched path relative to project root
    if os.path.isabs(args.enriched):
        enriched_path = args.enriched
    else:
        enriched_path = os.path.join(PROJECT_ROOT, args.enriched)

    manifest_path = None
    if args.manifest:
        if os.path.isabs(args.manifest):
            manifest_path = args.manifest
        else:
            manifest_path = os.path.join(PROJECT_ROOT, args.manifest)

    if args.scan_only:
        run_scan(
            enriched_path,
            args.out_dir,
            pvalue_threshold=args.pvalue_threshold,
            pass_mode=args.pass_mode,
            metrics=metrics,
        )
    elif args.chrom_mode:
        if not args.comparison or not args.chrom:
            parser.error("--comparison and --chrom required in chrom-mode")
        run_chrom(
            args.comparison,
            args.chrom,
            enriched_path,
            args.out_dir,
            manifest_path=manifest_path,
            pvalue_threshold=args.pvalue_threshold,
            pass_mode=args.pass_mode,
            metrics=metrics,
        )
    else:
        parser.error("Specify --scan-only or --chrom-mode")


if __name__ == "__main__":
    main()
