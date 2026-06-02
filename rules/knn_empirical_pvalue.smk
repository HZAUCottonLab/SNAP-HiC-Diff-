# KNN Empirical P-value Module
# Computes KNN-based empirical P-values for window pair metrics.

# ==============================================================================
# Knn_empirical_pvalue
# ==============================================================================

rule knn_empirical_pvalue:
    """Compute KNN-based empirical P-values for window pair metrics."""
    input:
        script = SCRIPT_DIR + "/scripts/compute_knn_empirical_pvalue.py",
        merged = OUTDIR + "/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv",
    output:
        enriched = OUTDIR + "/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
    log:
        OUTDIR + "/logs/knn_empirical_pvalue/{comparison}/knn_empirical_pvalue.log"
    benchmark:
        OUTDIR + "/benchmark/knn_empirical_pvalue/{comparison}/knn_empirical_pvalue.tsv"
    threads: config["parallel"]["threads_per_job"]
    params:
        background_file = config["knn_empirical_pvalue"]["background_file"],
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        COMPARISON="{wildcards.comparison}"

        mkdir -p "$(dirname "{output.enriched}")" "$(dirname "{log}")"

        python {input.script} \
            --comparison "$COMPARISON" \
            --input "{input.merged}" \
            --background "{params.background_file}" \
            --out "{output.enriched}" \
            --log "{log}" \
            > "{log}" 2>&1
        """
