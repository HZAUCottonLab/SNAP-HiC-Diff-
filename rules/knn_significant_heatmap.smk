# Significant KNN heatmaps stay within the same scale branch.

checkpoint plot_knn_significant_heatmap:
    input:
        script=SCRIPT_DIR + "/scripts/plot_knn_significant_heatmap.py",
        enriched=SCALE_ROOT + "/{scale}/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
    output:
        manifest=SCALE_ROOT + "/{scale}/logs/knn_significant_heatmap/{comparison}/significant_windows.tsv",
        pval_table=SCALE_ROOT + "/{scale}/logs/knn_significant_heatmap/{comparison}/pval_table.tsv",
        pval_summary=SCALE_ROOT + "/{scale}/logs/knn_significant_heatmap/{comparison}/pval_summary.tsv",
    log:
        SCALE_ROOT + "/{scale}/logs/knn_significant_heatmap/{comparison}/scan.log",
    benchmark:
        SCALE_ROOT + "/{scale}/benchmark/knn_significant_heatmap/{comparison}/scan.tsv",
    params:
        pvalue_threshold=config["knn_significant_heatmap"]["pvalue_threshold"],
        pass_mode=config["knn_significant_heatmap"]["pass_mode"],
        metrics=",".join(config["knn_significant_heatmap"]["metrics"]),
        common_outdir=OUTDIR,
        scale_outdir=get_scale_outdir,
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.manifest}")" "$(dirname "{log}")"
        python {input.script} \
            --scan-only \
            --enriched "{input.enriched}" \
            --out-dir "$(dirname "{output.manifest}")" \
            --common-outdir "{params.common_outdir}" \
            --scale-outdir "{params.scale_outdir}" \
            --pvalue-threshold "{params.pvalue_threshold}" \
            --pass-mode "{params.pass_mode}" \
            --metrics "{params.metrics}" \
            > "{log}" 2>&1
        """


rule plot_knn_chrom:
    """Produce a marker for every chromosome, including no-op chromosomes."""
    input:
        manifest=rules.plot_knn_significant_heatmap.output.manifest,
        script=SCRIPT_DIR + "/scripts/plot_knn_significant_heatmap.py",
        enriched=SCALE_ROOT + "/{scale}/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
    output:
        done=SCALE_ROOT + "/{scale}/logs/knn_significant_heatmap/{comparison}/{chrom}.done",
    log:
        SCALE_ROOT + "/{scale}/logs/knn_significant_heatmap/{comparison}/{chrom}.log",
    benchmark:
        SCALE_ROOT + "/{scale}/benchmark/knn_significant_heatmap/{comparison}/{chrom}.tsv",
    params:
        pvalue_threshold=config["knn_significant_heatmap"]["pvalue_threshold"],
        pass_mode=config["knn_significant_heatmap"]["pass_mode"],
        metrics=",".join(config["knn_significant_heatmap"]["metrics"]),
        out_dir=SCALE_ROOT + "/{scale}/knn_significant_heatmap/{comparison}",
        common_outdir=OUTDIR,
        scale_outdir=get_scale_outdir,
    threads:
        config["parallel"]["threads_per_job"]
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.out_dir}" "$(dirname "{output.done}")" "$(dirname "{log}")"
        python {input.script} \
            --chrom-mode \
            --comparison "{wildcards.comparison}" \
            --chrom "{wildcards.chrom}" \
            --enriched "{input.enriched}" \
            --manifest "{input.manifest}" \
            --out-dir "{params.out_dir}" \
            --common-outdir "{params.common_outdir}" \
            --scale-outdir "{params.scale_outdir}" \
            --pvalue-threshold "{params.pvalue_threshold}" \
            --pass-mode "{params.pass_mode}" \
            --metrics "{params.metrics}" \
            > "{log}" 2>&1
        touch "{output.done}"
        """
