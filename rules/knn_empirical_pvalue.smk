# KNN Empirical P-value Module
# Computes KNN-based empirical P-values for window pair metrics.

# ==============================================================================
# Knn_empirical_pvalue
# ==============================================================================

rule validate_knn_background:
    """Audit the configured KNN background for one active window scale."""
    input:
        script = SCRIPT_DIR + "/scripts/validate_knn_background.py",
        background = lambda w: str(get_background_preset(w.window_bins).file),
    output:
        audit = SCALE_ROOT + "/logs/knn_background/background_validation.json",
    log:
        SCALE_ROOT + "/logs/knn_background/background_validation.log",
    params:
        resolution_bp = RESOLUTION,
        window_bins = lambda w: int(w.window_bins),
        step_bins = lambda w: get_window_scale(w.window_bins).step_bins,
        knn_k = config["knn_empirical_pvalue"]["knn_k"],
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.audit}")" "$(dirname "{log}")"

        python {input.script} \
            --background "{input.background}" \
            --resolution-bp {params.resolution_bp} \
            --window-bins {params.window_bins} \
            --step-bins {params.step_bins} \
            --knn-k {params.knn_k} \
            --audit "{output.audit}" \
            > "{log}" 2>&1
        """

rule knn_empirical_pvalue:
    """Compute KNN-based empirical P-values for window pair metrics."""
    input:
        script = SCRIPT_DIR + "/scripts/compute_knn_empirical_pvalue.py",
        merged = SCALE_ROOT + "/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv",
        background = lambda w: str(get_background_preset(w.window_bins).file),
        audit = rules.validate_knn_background.output.audit,
    output:
        enriched = SCALE_ROOT + "/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
    log:
        SCALE_ROOT + "/logs/knn_empirical_pvalue/{comparison}/knn_empirical_pvalue.log"
    benchmark:
        SCALE_ROOT + "/benchmark/knn_empirical_pvalue/{comparison}/knn_empirical_pvalue.tsv"
    threads: config["parallel"]["threads_per_job"]
    params:
        config = ACTIVE_CONFIGFILE,
        resolution_bp = RESOLUTION,
        window_bins = lambda w: int(w.window_bins),
        window_size_bp = lambda w: get_window_scale(w.window_bins).window_size_bp,
        step_bins = lambda w: get_window_scale(w.window_bins).step_bins,
        step_size_bp = lambda w: get_window_scale(w.window_bins).step_size_bp,
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
            --background "{input.background}" \
            --config "{params.config}" \
            --resolution-bp {params.resolution_bp} \
            --window-bins {params.window_bins} \
            --window-size-bp {params.window_size_bp} \
            --step-bins {params.step_bins} \
            --step-size-bp {params.step_size_bp} \
            --out "{output.enriched}" \
            --log "{log}" \
            > "{log}" 2>&1
        """
