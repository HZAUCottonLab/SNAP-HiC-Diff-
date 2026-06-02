# Compute gap detection and per-window QC for all chromosomes of one sample.
# Shell block loops over all chroms internally.
# Outputs produced by the Python script at:
#   results/windows/{sample}/{sample}__{chrom}_gaps.tsv
#   results/windows/{sample}/{sample}__{chrom}_window_qc.tsv
#
# checkpoint: depends on raw data, runs the computation, creates .done marker.
# Downstream rules (sliding_window) depend on real output files directly.

checkpoint compute_gap_qc:
    input:
        ice_done = lambda w: checkpoints.extract_and_balance_sample.get(sample=w.sample).output.manifest,
        script = SCRIPT_DIR + "/scripts/compute_gap_qc.py"
    output:
        touch = OUTDIR + "/logs/gap_qc/{sample}.done",
        gaps = expand(
            OUTDIR + "/windows/{{sample}}/{{sample}}__{chrom}_gaps.tsv",
            chrom=CHROMS
        ),
        qc = expand(
            OUTDIR + "/windows/{{sample}}/{{sample}}__{chrom}_window_qc.tsv",
            chrom=CHROMS
        ),
    log:
        OUTDIR + "/logs/gap_qc/{sample}.log"
    benchmark:
        "{}/benchmark/gap_qc/{{sample}}.tsv".format(OUTDIR)
    threads:
        config["parallel"]["threads_per_job"]
    params:
        chroms = CHROMS,
        resolution = RESOLUTION,
        window_size = WINDOW_SIZE,
        step_size = STEP_SIZE,
        gradient_threshold = config["gap"]["gradient_threshold"],
        min_low_ratio = config["gap"]["min_low_ratio"],
        dynamic_raw_sum_cpm_pct = config["qc"]["dynamic_raw_sum_cpm_pct"],
        max_bad_bin_frac = config["qc"]["max_bad_bin_frac"],
        nan_abnormal_threshold = config["nan_density_breakpoint"]["abnormal_threshold"],
        outdir = OUTDIR
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"

        mkdir -p "{params.outdir}/windows/${{SAMPLE}}" "{params.outdir}/logs/gap_qc"

        BALANCED=""
        for chrom in {params.chroms}; do
            BASE="${{SAMPLE}}__${{chrom}}"
            BALANCED="$BALANCED {params.outdir}/balanced/${{SAMPLE}}/${{BASE}}_iced.cool"
        done

        python {input.script} \
            --chroms {params.chroms} \
            --sample "$SAMPLE" \
            --balanced-cool-files $BALANCED \
            --outdir "{params.outdir}/windows/${{SAMPLE}}" \
            --resolution {params.resolution} \
            --window-size {params.window_size} \
            --step-size {params.step_size} \
            --gradient-threshold {params.gradient_threshold} \
            --min-low-ratio {params.min_low_ratio} \
            --dynamic-raw-sum-cpm-pct {params.dynamic_raw_sum_cpm_pct} \
            --max-bad-bin-frac {params.max_bad_bin_frac} \
            --nan-abnormal-threshold {params.nan_abnormal_threshold} \
            --log-file "{log}" \
            --n-jobs {threads} \
            > "{log}" 2>&1

        touch "{output.touch}"
        """
