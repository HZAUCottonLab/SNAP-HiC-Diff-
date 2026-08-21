# Meltron insulation-level validation per comparison pair.
# Runs decondensation AND condensation, then merges results.
#
# Inputs:
#   - Final per-sample IS tables produced by merge_IS (sample A and sample B)
#   - knn_enriched.tsv from Phase 03 (candidate windows)
# Outputs:
#   - {comparison}_intervals.tsv (4-col: chrom, start, end, ID)
#   - {comparison}_interval_map.tsv (window_id mapping)
#   - {comparison}_decondensation.tsv
#   - {comparison}_condensation.tsv
#   - {comparison}_meltron_merged.tsv (final output with best_score, best_direction)

rule meltron_per_comparison:
    input:
        is_a=lambda w: OUTDIR + "/IS/{}/{}_IS.tsv".format(
            w.comparison.split("_vs_")[0],
            w.comparison.split("_vs_")[0],
        ),
        is_b=lambda w: OUTDIR + "/IS/{}/{}_IS.tsv".format(
            w.comparison.split("_vs_")[1],
            w.comparison.split("_vs_")[1],
        ),
        knn_enriched = SCALE_ROOT + "/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
        script_prepare = SCRIPT_DIR + "/scripts/prepare_meltron_intervals.py",
        script_run = SCRIPT_DIR + "/scripts/run_meltron.py",
        script_merge = SCRIPT_DIR + "/scripts/merge_meltron_scores.py",
    output:
        manifest = SCALE_ROOT + "/logs/meltron/{comparison}.done",
        intervals = SCALE_ROOT + "/meltron/{comparison}/{comparison}_intervals.tsv",
        interval_map = SCALE_ROOT + "/meltron/{comparison}/{comparison}_interval_map.tsv",
        decond = SCALE_ROOT + "/meltron/{comparison}/{comparison}_decondensation.tsv",
        condens = SCALE_ROOT + "/meltron/{comparison}/{comparison}_condensation.tsv",
        merged = SCALE_ROOT + "/meltron/{comparison}/{comparison}_meltron_merged.tsv",
    log:
        SCALE_ROOT + "/logs/meltron/{comparison}.log"
    benchmark:
        SCALE_ROOT + "/benchmark/meltron/{comparison}.tsv"
    threads:
        config["meltron"]["cores"]
    params:
        outdir = OUTDIR,
        cutoff_na = config["meltron"]["cutoff_na"],
        scale_root = lambda w: SCALE_ROOT.format(window_bins=w.window_bins),
        score_threshold = config["meltron"]["score_threshold"],
        min_significant = config["meltron"]["min_significant"],
        meltron_cores = config["meltron"]["cores"],
        script_prepare = SCRIPT_DIR + "/scripts/prepare_meltron_intervals.py",
        script_run = SCRIPT_DIR + "/scripts/run_meltron.py",
        script_merge = SCRIPT_DIR + "/scripts/merge_meltron_scores.py",
    conda:
        "../env/meltron.yaml"
    shell:
        r"""
        set -euo pipefail

        COMPARISON="{wildcards.comparison}"
        COMPARISON_STR="{wildcards.comparison}"
        SAMPLE_A="${{COMPARISON_STR%%_vs_*}}"
        SAMPLE_B="${{COMPARISON_STR#*_vs_}}"

        IS_A="{params.outdir}/IS/${{SAMPLE_A}}/${{SAMPLE_A}}_IS.tsv"
        IS_B="{params.outdir}/IS/${{SAMPLE_B}}/${{SAMPLE_B}}_IS.tsv"

        mkdir -p "{params.scale_root}/meltron/$COMPARISON"

        # Step 1: prepare Meltron intervals from knn_enriched.tsv
        python {params.script_prepare} \
            --enriched {input.knn_enriched} \
            --out {output.intervals} \
            --min-significant {params.min_significant} \
            --raw-pvalue-threshold 0.05 \
            --log-file "{log}" \
            >> "{log}" 2>&1

        # Step 2: decondensation Meltron
        python {params.script_run} \
            --is-a "$IS_A" \
            --is-b "$IS_B" \
            --genomic-intervals {output.intervals} \
            --outfile {output.decond} \
            --comparison-direction decondensation \
            --cutoff-na {params.cutoff_na} \
            --cores {params.meltron_cores} \
            --log-file "{log}" \
            >> "{log}" 2>&1

        # Step 3: condensation Meltron
        python {params.script_run} \
            --is-a "$IS_A" \
            --is-b "$IS_B" \
            --genomic-intervals {output.intervals} \
            --outfile {output.condens} \
            --comparison-direction condensation \
            --cutoff-na {params.cutoff_na} \
            --cores {params.meltron_cores} \
            --log-file "{log}" \
            >> "{log}" 2>&1

        # Step 4: merge both directions
        python {params.script_merge} \
            --comparison "$COMPARISON" \
            --decondensation {output.decond} \
            --condensation {output.condens} \
            --interval-map {output.interval_map} \
            --out {output.merged} \
            --score-threshold {params.score_threshold} \
            --log-file "{log}" \
            >> "{log}" 2>&1

        touch {output.manifest}
        """
