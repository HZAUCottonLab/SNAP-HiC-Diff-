# Meltron validation runs independently for every window scale.

rule meltron_per_comparison:
    input:
        is_a=lambda w: get_comparison_is_file(w.comparison, 0),
        is_b=lambda w: get_comparison_is_file(w.comparison, 1),
        knn_enriched=SCALE_ROOT + "/{scale}/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
        script_prepare=SCRIPT_DIR + "/scripts/prepare_meltron_intervals.py",
        script_run=SCRIPT_DIR + "/scripts/run_meltron.py",
        script_merge=SCRIPT_DIR + "/scripts/merge_meltron_scores.py",
    output:
        manifest=SCALE_ROOT + "/{scale}/logs/meltron/{comparison}.done",
        intervals=SCALE_ROOT + "/{scale}/meltron/{comparison}/{comparison}_intervals.tsv",
        interval_map=SCALE_ROOT + "/{scale}/meltron/{comparison}/{comparison}_interval_map.tsv",
        decond=SCALE_ROOT + "/{scale}/meltron/{comparison}/{comparison}_decondensation.tsv",
        condens=SCALE_ROOT + "/{scale}/meltron/{comparison}/{comparison}_condensation.tsv",
        merged=SCALE_ROOT + "/{scale}/meltron/{comparison}/{comparison}_meltron_merged.tsv",
    log:
        SCALE_ROOT + "/{scale}/logs/meltron/{comparison}.log",
    benchmark:
        SCALE_ROOT + "/{scale}/benchmark/meltron/{comparison}.tsv",
    threads:
        config["meltron"]["cores"]
    params:
        scale_outdir=get_scale_outdir,
        cutoff_na=config["meltron"]["cutoff_na"],
        score_threshold=config["meltron"]["score_threshold"],
        min_significant=config["meltron"]["min_significant"],
        meltron_cores=config["meltron"]["cores"],
    conda:
        "../env/meltron.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.intervals}")" "$(dirname "{log}")"

        python {input.script_prepare} \
            --enriched "{input.knn_enriched}" \
            --out "{output.intervals}" \
            --min-significant {params.min_significant} \
            --fdr-threshold 0.05 \
            --log-file "{log}" \
            >> "{log}" 2>&1

        python {input.script_run} \
            --is-a "{input.is_a}" \
            --is-b "{input.is_b}" \
            --genomic-intervals "{output.intervals}" \
            --outfile "{output.decond}" \
            --comparison-direction decondensation \
            --cutoff-na {params.cutoff_na} \
            --cores {params.meltron_cores} \
            --log-file "{log}" \
            >> "{log}" 2>&1

        python {input.script_run} \
            --is-a "{input.is_a}" \
            --is-b "{input.is_b}" \
            --genomic-intervals "{output.intervals}" \
            --outfile "{output.condens}" \
            --comparison-direction condensation \
            --cutoff-na {params.cutoff_na} \
            --cores {params.meltron_cores} \
            --log-file "{log}" \
            >> "{log}" 2>&1

        python {input.script_merge} \
            --comparison "{wildcards.comparison}" \
            --decondensation "{output.decond}" \
            --condensation "{output.condens}" \
            --interval-map "{output.interval_map}" \
            --out "{output.merged}" \
            --score-threshold {params.score_threshold} \
            --log-file "{log}" \
            >> "{log}" 2>&1

        touch "{output.manifest}"
        """
