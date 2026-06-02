# KNN Significant Window 2x3 Heatmap Module
# ==============================================================================
# Architecture: 3-tier
#   checkpoint: scans knn_enriched.tsv (raw p-value based), writes:
#     - significant_windows.tsv   (manifest of all-pass windows)
#     - pval_table.tsv            (all windows, per-metric raw p-values + pass/fail)
#     - pval_summary.tsv          (per-metric summary statistics)
#   rule plot_knn_chrom: per-chrom done marker + 2x3 heatmap PDF/PNG
#
# Data sources:
#   iced .cool:   results/balanced/{sample}/{sample}__{chrom}_iced.cool   (per-chrom)
#   ln_oe h5:     results/windows/{sample}/{sample}__{chrom}_ln_oe.h5     (per-window)
#   gauss h5:     results/windows/{sample}/{sample}__{chrom}_gaussian.h5  (per-window)
#   knn enriched: results/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv
# ==============================================================================

checkpoint plot_knn_significant_heatmap:
    input:
        script = SCRIPT_DIR + "/scripts/plot_knn_significant_heatmap.py",
        enriched = OUTDIR + "/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
    output:
        manifest = OUTDIR + "/logs/knn_significant_heatmap/{comparison}/significant_windows.tsv",
        pval_table = OUTDIR + "/logs/knn_significant_heatmap/{comparison}/pval_table.tsv",
        pval_summary = OUTDIR + "/logs/knn_significant_heatmap/{comparison}/pval_summary.tsv",
    log:
        OUTDIR + "/logs/knn_significant_heatmap/{comparison}/scan.log",
    benchmark:
        OUTDIR + "/benchmark/knn_significant_heatmap/{comparison}/scan.tsv",
    params:
        pvalue_threshold = config["knn_significant_heatmap"]["pvalue_threshold"],
        pass_mode = config["knn_significant_heatmap"]["pass_mode"],
        metrics = ",".join(config["knn_significant_heatmap"]["metrics"]),
    conda:
        "../env/SNAP-HiC-Diff.yaml",
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --scan-only \
            --enriched "{input.enriched}" \
            --out-dir "$(dirname "{output.manifest}")" \
            --pvalue-threshold "{params.pvalue_threshold}" \
            --pass-mode "{params.pass_mode}" \
            --metrics "{params.metrics}" \
            > "{log}" 2>&1
        """


# Per-chrom done marker: one Snakemake job per chromosome (always scheduled);
# script will no-op + write done if no significant windows for that chrom.
rule plot_knn_chrom:
    input:
        manifest = rules.plot_knn_significant_heatmap.output.manifest,
        script = SCRIPT_DIR + "/scripts/plot_knn_significant_heatmap.py",
        enriched = OUTDIR + "/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
    output:
        done = OUTDIR + "/logs/knn_significant_heatmap/{comparison}/{chrom}.done",
    log:
        OUTDIR + "/logs/knn_significant_heatmap/{comparison}/{chrom}.log",
    benchmark:
        OUTDIR + "/benchmark/knn_significant_heatmap/{comparison}/{chrom}.tsv",
    params:
        pvalue_threshold = config["knn_significant_heatmap"]["pvalue_threshold"],
        pass_mode = config["knn_significant_heatmap"]["pass_mode"],
        metrics = ",".join(config["knn_significant_heatmap"]["metrics"]),
        out_dir = OUTDIR + "/knn_significant_heatmap/{comparison}",
    threads:
        config["parallel"]["threads_per_job"],
    conda:
        "../env/SNAP-HiC-Diff.yaml",
    shell:
        r"""
        set -euo pipefail
        COMPARISON="{wildcards.comparison}"
        CHROM="{wildcards.chrom}"
        python {input.script} \
            --chrom-mode \
            --comparison "$COMPARISON" \
            --chrom "$CHROM" \
            --enriched "{input.enriched}" \
            --manifest "{input.manifest}" \
            --out-dir "{params.out_dir}" \
            --pvalue-threshold "{params.pvalue_threshold}" \
            --pass-mode "{params.pass_mode}" \
            --metrics "{params.metrics}" \
            > "{log}" 2>&1
        touch "{output.done}"
        """
