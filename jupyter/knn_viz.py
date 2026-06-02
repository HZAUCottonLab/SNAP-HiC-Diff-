#!/usr/bin/env python3
"""KNN empirical P-value KDE visualization — minimal, memory-safe."""

import os, sys, glob, argparse, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
warnings.filterwarnings("ignore")

# ── Detect project root from file location ────────────────────────────────────
_JUPYTER_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(_JUPYTER_DIR)   # parent of jupyter/


class KNNVizPlotter:
    """KDE-density plots for KNN empirical P-values."""

    SCORE_COLS = [
        "MSE_empirical_pvalue",
        "d_SSIM_empirical_pvalue",
        "d_Spearman_empirical_pvalue",
        "tw_score_empirical_pvalue",
        "d_NCC_empirical_pvalue",
    ]
    SCORE_LABELS = {
        "MSE_empirical_pvalue": "MSE",
        "d_SSIM_empirical_pvalue": "d_SSIM",
        "d_Spearman_empirical_pvalue": "d_Spearman",
        "tw_score_empirical_pvalue": "tw_score",
        "d_NCC_empirical_pvalue": "d_NCC",
    }
    BG_COLS = {
        "MSE_empirical_pvalue": "bg_MSE_median",
        "d_SSIM_empirical_pvalue": "bg_d_SSIM_median",
        "d_Spearman_empirical_pvalue": "bg_d_Spearman_median",
        "tw_score_empirical_pvalue": "bg_tw_score_median",
        "d_NCC_empirical_pvalue": "bg_d_NCC_median",
    }

    def __init__(self, df, comparison, outdir):
        self.df = df
        self.comparison = comparison
        self.outdir = os.path.join(outdir, comparison)
        self.saved_paths = []
        os.makedirs(self.outdir, exist_ok=True)
        plt.rcParams.update({
            "font.size": 12,
            "figure.dpi": 150,
            "axes.labelsize": 11,
            "axes.titlesize": 13,
        })

    def _save(self, name, fig):
        png = os.path.join(self.outdir, f"{name}.png")
        svg = os.path.join(self.outdir, f"{name}.svg")
        fig.savefig(png, dpi=150, bbox_inches="tight")
        fig.savefig(svg, bbox_inches="tight")
        self.saved_paths.append((name, png, svg))
        plt.close(fig)
        plt.close("all")

    def _log10(self, x, floor=1e-30):
        return -np.log10(np.maximum(np.asarray(x, dtype=float), floor))

    # ── 1. KDE per score ──────────────────────────────────────────────────────
    def kde_per_score(self):
        fig, axes = plt.subplots(1, 5, figsize=(20, 4))
        fig.suptitle(f"KDE: P-value Density — {self.comparison}", fontsize=14)
        for i, (ax, col) in enumerate(zip(axes, self.SCORE_COLS)):
            data = self.df[col].dropna().values
            sns.kdeplot(data, ax=ax, fill=True, alpha=0.4, color=f"C{i}", linewidth=2)
            ax.set_xlabel(self.SCORE_LABELS[col])
            ax.set_ylabel("Density")
            ax.set_xscale("log")
            ax.axvline(np.median(data), color="crimson", lw=1.5, ls="--",
                       label=f"med={np.median(data):.2g}")
            ax.legend(fontsize=8)
        plt.tight_layout()
        self._save("01_kde_per_score", fig)

    # ── 2. Chrom KDE heatmap ────────────────────────────────────────────────
    def chrom_kde_heatmap(self, score_col=None):
        if score_col is None:
            score_col = self.SCORE_COLS[0]
        label = self.SCORE_LABELS[score_col]
        chroms = sorted(self.df["chrom"].unique())
        n_bins = 30
        matrix = []
        for chrom in chroms:
            vals = self.df.loc[self.df["chrom"] == chrom, score_col].dropna().values
            if len(vals) < 5:
                matrix.append(np.zeros(n_bins))
            else:
                hist, _ = np.histogram(self._log10(vals), bins=n_bins, density=True)
                matrix.append(hist)
        fig, ax = plt.subplots(figsize=(10, max(5, len(chroms) * 0.3)))
        mat = np.array(matrix)
        im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", interpolation="nearest")
        ax.set_yticks(range(len(chroms)))
        ax.set_yticklabels(chroms, fontsize=8)
        ax.set_xlabel("-log10(p) bin")
        ax.set_ylabel("Chromosome")
        ax.set_title(f"Chrom KDE Heatmap — {label} — {self.comparison}")
        plt.colorbar(im, ax=ax, label="density")
        plt.tight_layout()
        self._save(f"02_chrom_kde_heatmap_{label}", fig)

    # ── 3. Score-pair KDE hexbin ──────────────────────────────────────────
    def hexbin_pairs(self):
        fig, axes = plt.subplots(5, 5, figsize=(16, 15))
        fig.suptitle(f"Pairwise Hexbin Density — {self.comparison}", fontsize=14)
        n_sub = min(5000, len(self.df))
        df_sub = self.df.sample(n=n_sub, random_state=42)
        for i, ci in enumerate(self.SCORE_COLS):
            for j, cj in enumerate(self.SCORE_COLS):
                ax = axes[i, j]
                xi = self._log10(df_sub[ci].values)
                yi = self._log10(df_sub[cj].values)
                m = np.isfinite(xi) & np.isfinite(yi)
                if i == j:
                    sns.kdeplot(df_sub[ci].dropna().values, ax=ax,
                                fill=True, alpha=0.4, color=f"C{i}")
                else:
                    ax.hexbin(xi[m], yi[m], gridsize=30, cmap="YlOrRd", mincnt=1, bins="log")
                if i == 0:
                    ax.set_title(self.SCORE_LABELS[cj], fontsize=9)
                if j == 0:
                    ax.set_ylabel(self.SCORE_LABELS[ci], fontsize=9)
                ax.tick_params(labelsize=7)
        plt.tight_layout()
        self._save("03_hexbin_pairs", fig)

    # ── 4. knn_dist vs score: hexbin ────────────────────────────────────
    def knn_vs_score(self):
        n_sub = min(10000, len(self.df))
        df_s = self.df.sample(n=n_sub, random_state=42)
        fig, axes = plt.subplots(1, 5, figsize=(20, 4))
        fig.suptitle(f"KNN Distance vs P-value — {self.comparison}", fontsize=14)
        dist = df_s["knn_dist_median"].values
        for i, (ax, col) in enumerate(zip(axes, self.SCORE_COLS)):
            pv = self._log10(df_s[col].values)
            m = np.isfinite(dist) & np.isfinite(pv)
            ax.hexbin(dist[m], pv[m], gridsize=40, cmap="YlOrRd", mincnt=1, bins="log")
            ax.set_xlabel("knn_dist_median")
            ax.set_ylabel(r"$-\log_{10}(p)$")
            ax.set_title(self.SCORE_LABELS[col], fontsize=10)
        plt.tight_layout()
        self._save("04_knn_vs_score", fig)

    # ── 5. PCA: scatter + density contour ─────────────────────────────────
    def pca_contour(self):
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        data = self.df[self.SCORE_COLS].dropna()
        if len(data) > 15000:
            data = data.sample(15000, random_state=42)
        X = StandardScaler().fit_transform(data)
        Xp = PCA(n_components=2).fit_transform(X)
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(Xp[:, 0], Xp[:, 1], s=1.5, alpha=0.3, c="steelblue")
        try:
            sns.kdeplot(x=Xp[:, 0], y=Xp[:, 1], ax=ax, levels=5, color="crimson", linewidths=1.5)
        except Exception:
            pass
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_title(f"PCA — {self.comparison}")
        plt.tight_layout()
        self._save("05_pca_contour", fig)

    # ── 6. Correlation heatmap ────────────────────────────────────────────
    def corr_heatmap(self):
        corr = self.df[self.SCORE_COLS].corr()
        labels = [self.SCORE_LABELS[c] for c in self.SCORE_COLS]
        fig, ax = plt.subplots(figsize=(7, 6))
        sns.heatmap(corr, annot=True, fmt=".2f", cmap="RdBu_r",
                    vmin=-1, vmax=1, xticklabels=labels, yticklabels=labels,
                    ax=ax, linewidths=0.5, annot_kws={"size": 11})
        ax.set_title(f"P-value Correlation — {self.comparison}")
        plt.tight_layout()
        self._save("06_corr_heatmap", fig)

    # ── 7. BG deviation KDE ────────────────────────────────────────────────
    def bg_dev_kde(self):
        fig, axes = plt.subplots(1, 5, figsize=(20, 4))
        fig.suptitle(f"Background Deviation KDE — {self.comparison}", fontsize=14)
        for i, (ax, col) in enumerate(zip(axes, self.SCORE_COLS)):
            bg = self.BG_COLS[col]
            dev = ((self.df[col] - self.df[bg]) /
                   np.maximum(self.df[bg], 1e-10)).dropna().values
            dev = dev[np.isfinite(dev)]
            sns.kdeplot(dev, ax=ax, fill=True, alpha=0.4, color=f"C{i}", linewidth=2)
            ax.axvline(0, color="black", lw=1.5, ls="--")
            ax.set_xlabel(f"({self.SCORE_LABELS[col]} - bg) / bg")
            ax.set_ylabel("Density")
            ax.set_title(f"{self.SCORE_LABELS[col]} dev", fontsize=10)
        plt.tight_layout()
        self._save("07_bg_dev_kde", fig)

    # ── 8. Summary stats table ────────────────────────────────────────────
    def summary_table(self):
        n_win = len(self.df)
        n_chr = self.df["chrom"].nunique()
        lines = [
            f"Comparison: {self.comparison}",
            f"Windows: {n_win:,}  |  Chroms: {n_chr}",
            "",
            "Metric     %p<0.05  %p<0.01  %p<0.001  median_p",
            ":---       ------    ------    -------   --------",
        ]
        for col in self.SCORE_COLS:
            pv = self.df[col]
            lbl = self.SCORE_LABELS[col]
            lines.append(f"{lbl:<10} "
                         f"{(pv<0.05).mean()*100:>6.1f}%  "
                         f"{(pv<0.01).mean()*100:>6.1f}%  "
                         f"{(pv<0.001).mean()*100:>7.1f}%  "
                         f"{pv.median():>8.3g}")
        lines += [
            "",
            f"knn_valid_k: [{self.df['knn_valid_k'].min():,} – {self.df['knn_valid_k'].max():,}]",
            f"knn_dist_median: [{self.df['knn_dist_median'].min():.4f} – {self.df['knn_dist_median'].max():.4f}]",
        ]
        fig, ax = plt.subplots(figsize=(12, 5))
        ax.axis("off")
        ax.text(0.05, 0.95, "\n".join(lines), transform=ax.transAxes,
                fontsize=12, fontfamily="monospace", verticalalignment="top")
        ax.set_title(f"Summary — {self.comparison}", fontsize=14)
        plt.tight_layout()
        self._save("08_summary_table", fig)


def auto_discover(data_dir=None):
    if data_dir is None:
        data_dir = os.path.join(PROJECT_ROOT, "results/window_compare")
    if not os.path.isabs(data_dir):
        data_dir = os.path.join(PROJECT_ROOT, data_dir)
    # Files are named: {comparison}_knn_enriched.tsv
    pat = os.path.join(data_dir, "*", "*_knn_enriched.tsv")
    dirs = sorted(set(os.path.basename(os.path.dirname(p)) for p in glob.glob(pat)))
    return dirs


def load_data(comparison, data_dir=None):
    if data_dir is None:
        data_dir = os.path.join(PROJECT_ROOT, "results/window_compare")
    if not os.path.isabs(data_dir):
        data_dir = os.path.join(PROJECT_ROOT, data_dir)
    path = os.path.join(data_dir, comparison, f"{comparison}_knn_enriched.tsv")
    df = pd.read_csv(path, sep="\t")
    for col in KNNVizPlotter.SCORE_COLS + ["knn_dist_median"]:
        if col in df.columns:
            df = df[np.isfinite(df[col])]
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--comparison", default=None)
    ap.add_argument("-d", "--data-dir",
                    default=os.path.join(PROJECT_ROOT, "results/window_compare"))
    ap.add_argument("-o", "--outdir",
                    default=os.path.join(PROJECT_ROOT, "results/figures"))
    args = ap.parse_args()

    comps = [args.comparison] if args.comparison else auto_discover(args.data_dir)
    print(f"Comparisons: {comps}")

    for comp in comps:
        print(f"\n-- {comp}")
        df = load_data(comp, args.data_dir)
        print(f"  {len(df):,} rows × {df['chrom'].nunique()} chroms")
        p = KNNVizPlotter(df, comp, args.outdir)

        p.kde_per_score()
        print("  ✓ kde_per_score")

        for col in KNNVizPlotter.SCORE_COLS:
            p.chrom_kde_heatmap(col)
        print("  ✓ chrom_kde_heatmap (5)")

        p.hexbin_pairs()
        print("  ✓ hexbin_pairs")

        p.knn_vs_score()
        print("  ✓ knn_vs_score")

        p.pca_contour()
        print("  ✓ pca_contour")

        p.corr_heatmap()
        print("  ✓ corr_heatmap")

        p.bg_dev_kde()
        print("  ✓ bg_dev_kde")

        p.summary_table()
        print("  ✓ summary_table")

        print(f"  → {len(p.saved_paths)} figures saved")

    print("\nDone.")


if __name__ == "__main__":
    main()
