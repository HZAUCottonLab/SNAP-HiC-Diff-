# Sliding window processing: gap fill, z-score, and Gaussian smoothing.
# Two-stage architecture:
#   checkpoint: guarantees gap.tsv files exist, creates .done marker
#   rule sliding_window_h5: declares H5 output files, depends on checkpoint
#
# Inputs from compute_gap_qc at: results/windows/{sample}/{sample}__{chrom}_gaps.tsv
# Inputs from compute_oe_matrix at: results/oe/{sample}/{sample}__{chrom}_oe.npz
# Outputs produced at: results/windows/{sample}/{sample}__{chrom}_gap_filled.h5 etc.

checkpoint sliding_window:
    input:
        gaps = expand(
            OUTDIR + "/windows/{{sample}}/{{sample}}__{chrom}_gaps.tsv",
            chrom=CHROMS
        ),
        oe_files = lambda w: [
            "{}/oe/{}/{}__{}_oe.npz".format(OUTDIR, w.sample, w.sample, chrom)
            for chrom in CHROMS
        ],
        script = SCRIPT_DIR + "/scripts/sliding_window.py"
    output:
        manifest = OUTDIR + "/logs/windows/{sample}.done"
    log:
        OUTDIR + "/logs/windows/{sample}.log"
    benchmark:
        "{}/benchmark/windows/{{sample}}.tsv".format(OUTDIR)
    threads:
        config["parallel"]["threads_per_job"]
    params:
        chroms = CHROMS,
        resolution = RESOLUTION,
        window_size = WINDOW_SIZE,
        step_size = STEP_SIZE,
        sigma = config["smoothing"]["sigma"],
        outdir = OUTDIR
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"

        mkdir -p "{params.outdir}/windows/${{SAMPLE}}" "{params.outdir}/logs/windows"

        GAPS=""
        for chrom in {params.chroms}; do
            BASE="${{SAMPLE}}__${{chrom}}"
            GAPS="$GAPS {params.outdir}/windows/${{SAMPLE}}/${{BASE}}_gaps.tsv"
        done

        OE=""
        for chrom in {params.chroms}; do
            BASE="${{SAMPLE}}__${{chrom}}"
            OE="$OE {params.outdir}/oe/${{SAMPLE}}/${{BASE}}_oe.npz"
        done

        python {input.script} \
            --chroms {params.chroms} \
            --sample "$SAMPLE" \
            --gap-tsv-files $GAPS \
            --ln-oe-files $OE \
            --outdir "{params.outdir}/windows/${{SAMPLE}}" \
            --resolution {params.resolution} \
            --window-size {params.window_size} \
            --step-size {params.step_size} \
            --sigma {params.sigma} \
            --log-file "{log}" \
            --n-jobs {threads} \
            > "{log}" 2>&1

        touch "{output.manifest}"
        """


# Declares the H5 files produced by the sliding_window checkpoint.
# Downstream rules depend on this rule's output to get mtime tracking.
rule sliding_window_h5:
    """Declare H5 output files for Snakemake dependency tracking."""
    input:
        lambda w: checkpoints.sliding_window.get(sample=w.sample).output.manifest
    output:
        h5_files = expand(
            OUTDIR + "/windows/{{sample}}/{{sample}}__{{chrom}}_{{stage}}.h5",
            chrom=CHROMS,
            stage=["ln_oe", "gap_filled", "zscore", "gaussian"]
        )
