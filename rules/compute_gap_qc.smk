# Compute gap detection and per-window QC for all chromosomes of one sample.
# Shell block loops over all chroms internally.
# Outputs produced by the Python script at:
#   results/window_scales/r<resolution>/w<window_bins>/windows/{sample}/{sample}__{chrom}_gaps.tsv
#   results/window_scales/r<resolution>/w<window_bins>/windows/{sample}/{sample}__{chrom}_window_qc.tsv
#
# checkpoint: depends on raw data, runs the computation, creates .done marker.
# Downstream rules (sliding_window) depend on real output files directly.

checkpoint compute_gap_qc:
    input:
        ice_done = lambda w: checkpoints.extract_and_balance_sample.get(sample=w.sample).output.manifest,
        script = SCRIPT_DIR + "/scripts/compute_gap_qc.py"
    output:
        touch = SCALE_ROOT + "/logs/gap_qc/{sample}.done",
        gaps = expand(
            SCALE_ROOT + "/windows/{{sample}}/{{sample}}__{chrom}_gaps.tsv",
            chrom=CHROMS, window_bins="{{window_bins}}"
        ),
        qc = expand(
            SCALE_ROOT + "/windows/{{sample}}/{{sample}}__{chrom}_window_qc.tsv",
            chrom=CHROMS, window_bins="{{window_bins}}"
        ),
    log:
        SCALE_ROOT + "/logs/gap_qc/{sample}.log"
    benchmark:
        SCALE_ROOT + "/benchmark/gap_qc/{sample}.tsv"
    threads:
        config["parallel"]["threads_per_job"]
    params:
        chroms = CHROMS,
        resolution = RESOLUTION,
        window_size = lambda w: get_window_scale(w.window_bins).window_size_bp,
        step_size = lambda w: get_window_scale(w.window_bins).step_size_bp,
        gradient_threshold = config["gap"]["gradient_threshold"],
        min_low_ratio = config["gap"]["min_low_ratio"],
        dynamic_raw_sum_cpm_pct = config["qc"]["dynamic_raw_sum_cpm_pct"],
        max_bad_bin_frac = config["qc"]["max_bad_bin_frac"],
        nan_abnormal_threshold = config["nan_density_breakpoint"]["abnormal_threshold"],
        scale_root = lambda w: SCALE_ROOT.format(window_bins=w.window_bins),
        shared_outdir = OUTDIR
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"

        mkdir -p "{params.scale_root}/windows/${{SAMPLE}}" "{params.scale_root}/logs/gap_qc" \
            "{params.scale_root}/benchmark/gap_qc"

        BALANCED=""
        for chrom in {params.chroms}; do
            BASE="${{SAMPLE}}__${{chrom}}"
            BALANCED="$BALANCED {params.shared_outdir}/balanced/${{SAMPLE}}/${{BASE}}_iced.cool"
        done

        python {input.script} \
            --chroms {params.chroms} \
            --sample "$SAMPLE" \
            --balanced-cool-files $BALANCED \
            --outdir "{params.scale_root}/windows/${{SAMPLE}}" \
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
