# Compute gap detection and per-window QC for one sample at one window scale.

checkpoint compute_gap_qc:
    input:
        ice_done=lambda w: checkpoints.extract_and_balance_sample.get(
            sample=w.sample
        ).output.manifest,
        script=SCRIPT_DIR + "/scripts/compute_gap_qc.py",
    output:
        touch=SCALE_ROOT + "/{scale}/logs/gap_qc/{sample}.done",
        gaps=expand(
            SCALE_ROOT + "/{{scale}}/windows/{{sample}}/{{sample}}__{chrom}_gaps.tsv",
            chrom=CHROMS,
        ),
        qc=expand(
            SCALE_ROOT + "/{{scale}}/windows/{{sample}}/{{sample}}__{chrom}_window_qc.tsv",
            chrom=CHROMS,
        ),
    log:
        SCALE_ROOT + "/{scale}/logs/gap_qc/{sample}.log",
    benchmark:
        SCALE_ROOT + "/{scale}/benchmark/gap_qc/{sample}.tsv",
    threads:
        config["parallel"]["threads_per_job"]
    params:
        chroms=CHROMS,
        resolution=lambda w: get_scale_value(w, "resolution_bp"),
        window_size=lambda w: get_scale_value(w, "window_size_bp"),
        step_size=lambda w: get_scale_value(w, "step_size_bp"),
        scale_outdir=get_scale_outdir,
        gradient_threshold=config["gap"]["gradient_threshold"],
        min_low_ratio=config["gap"]["min_low_ratio"],
        dynamic_raw_sum_cpm_pct=config["qc"]["dynamic_raw_sum_cpm_pct"],
        max_bad_bin_frac=config["qc"]["max_bad_bin_frac"],
        nan_abnormal_threshold=config["nan_density_breakpoint"]["abnormal_threshold"],
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"
        mkdir -p "{params.scale_outdir}/windows/${{SAMPLE}}" "{params.scale_outdir}/logs/gap_qc"

        BALANCED=""
        for chrom in {params.chroms}; do
            BASE="${{SAMPLE}}__${{chrom}}"
            BALANCED="$BALANCED {OUTDIR}/balanced/${{SAMPLE}}/${{BASE}}_iced.cool"
        done

        python {input.script} \
            --chroms {params.chroms} \
            --sample "$SAMPLE" \
            --balanced-cool-files $BALANCED \
            --outdir "{params.scale_outdir}/windows/${{SAMPLE}}" \
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
