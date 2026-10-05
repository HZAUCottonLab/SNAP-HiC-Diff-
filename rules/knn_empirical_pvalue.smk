# KNN empirical p-values use the exact validated background for each scale.

rule validate_knn_background:
    input:
        script=SCRIPT_DIR + "/scripts/validate_knn_background.py",
        background=lambda w: get_scale_value(w, "background_file"),
    output:
        audit=SCALE_ROOT + "/{scale}/audits/knn_background.json",
    log:
        SCALE_ROOT + "/{scale}/logs/knn_empirical_pvalue/background_validation.log",
    params:
        resolution=lambda w: get_scale_value(w, "resolution_bp"),
        window_bins=lambda w: get_scale_value(w, "window_bins"),
        step_bins=lambda w: get_scale_value(w, "step_bins"),
        knn_k=config["knn_empirical_pvalue"]["knn_k"],
    threads:
        4
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.audit}")" "$(dirname "{log}")"
        python "{input.script}" --background "{input.background}" \
          --resolution-bp {params.resolution} \
          --window-bins {params.window_bins} \
          --step-bins {params.step_bins} \
          --knn-k {params.knn_k} --out "{output.audit}" > "{log}" 2>&1
        """


rule knn_empirical_pvalue:
    """Compute empirical p-values for one comparison and scale."""
    input:
        script=SCRIPT_DIR + "/scripts/compute_knn_empirical_pvalue.py",
        audit=rules.validate_knn_background.output.audit,
        background=lambda w: get_scale_value(w, "background_file"),
        merged=SCALE_ROOT + "/{scale}/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv",
    output:
        enriched=SCALE_ROOT + "/{scale}/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
    log:
        SCALE_ROOT + "/{scale}/logs/knn_empirical_pvalue/{comparison}/knn_empirical_pvalue.log",
    benchmark:
        SCALE_ROOT + "/{scale}/benchmark/knn_empirical_pvalue/{comparison}/knn_empirical_pvalue.tsv",
    threads:
        config["parallel"]["threads_per_job"]
    params:
        resolution=lambda w: get_scale_value(w, "resolution_bp"),
        window_bins=lambda w: get_scale_value(w, "window_bins"),
        step_bins=lambda w: get_scale_value(w, "step_bins"),
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.enriched}")" "$(dirname "{log}")"
        python {input.script} \
            --comparison "{wildcards.comparison}" \
            --input "{input.merged}" \
            --background "{input.background}" \
            --config config/config.yaml \
            --resolution-bp {params.resolution} \
            --window-bins {params.window_bins} \
            --step-bins {params.step_bins} \
            --out "{output.enriched}" \
            --log "{log}" \
            > "{log}" 2>&1
        """
