# Sliding-window processing for one sample at one scale.

checkpoint sliding_window:
    input:
        gaps=expand(
            SCALE_ROOT + "/{{scale}}/windows/{{sample}}/{{sample}}__{chrom}_gaps.tsv",
            chrom=CHROMS,
        ),
        oe_files=lambda w: [
            "{}/oe/{}/{}__{}_oe.npz".format(OUTDIR, w.sample, w.sample, chrom)
            for chrom in CHROMS
        ],
        script=SCRIPT_DIR + "/scripts/sliding_window.py",
    output:
        manifest=SCALE_ROOT + "/{scale}/logs/windows/{sample}.done",
    log:
        SCALE_ROOT + "/{scale}/logs/windows/{sample}.log",
    benchmark:
        SCALE_ROOT + "/{scale}/benchmark/windows/{sample}.tsv",
    threads:
        config["parallel"]["threads_per_job"]
    params:
        chroms=CHROMS,
        resolution=lambda w: get_scale_value(w, "resolution_bp"),
        window_size=lambda w: get_scale_value(w, "window_size_bp"),
        step_size=lambda w: get_scale_value(w, "step_size_bp"),
        scale_outdir=get_scale_outdir,
        sigma=config["smoothing"]["sigma"],
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"
        mkdir -p "{params.scale_outdir}/windows/${{SAMPLE}}" "{params.scale_outdir}/logs/windows"

        GAPS=""
        OE=""
        for chrom in {params.chroms}; do
            BASE="${{SAMPLE}}__${{chrom}}"
            GAPS="$GAPS {params.scale_outdir}/windows/${{SAMPLE}}/${{BASE}}_gaps.tsv"
            OE="$OE {OUTDIR}/oe/${{SAMPLE}}/${{BASE}}_oe.npz"
        done

        python {input.script} \
            --chroms {params.chroms} \
            --sample "$SAMPLE" \
            --gap-tsv-files $GAPS \
            --ln-oe-files $OE \
            --outdir "{params.scale_outdir}/windows/${{SAMPLE}}" \
            --resolution {params.resolution} \
            --window-size {params.window_size} \
            --step-size {params.step_size} \
            --sigma {params.sigma} \
            --log-file "{log}" \
            --n-jobs {threads} \
            > "{log}" 2>&1

        touch "{output.manifest}"
        """


rule sliding_window_h5:
    """Declare scale-qualified H5 files from the sliding-window checkpoint."""
    input:
        lambda w: checkpoints.sliding_window.get(
            scale=w.scale, sample=w.sample
        ).output.manifest
    output:
        h5_files=expand(
            SCALE_ROOT + "/{{scale}}/windows/{{sample}}/{{sample}}__{{chrom}}_{{stage}}.h5",
            chrom=CHROMS,
            stage=["ln_oe", "gap_filled", "zscore", "gaussian"],
        )
    shell:
        "true"
