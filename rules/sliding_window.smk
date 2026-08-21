# Sliding window processing: gap fill, z-score, and Gaussian smoothing.
# Two-stage architecture:
#   checkpoint: guarantees gap.tsv files exist, creates .done marker
#   rule sliding_window_h5: declares H5 output files, depends on checkpoint
#
# Inputs from compute_gap_qc at: results/window_scales/r<resolution>/w<window_bins>/windows/{sample}/{sample}__{chrom}_gaps.tsv
# Inputs from compute_oe_matrix at: results/oe/{sample}/{sample}__{chrom}_oe.npz
# Outputs produced at: results/window_scales/r<resolution>/w<window_bins>/windows/{sample}/{sample}__{chrom}_gap_filled.h5 etc.

checkpoint sliding_window:
    input:
        gaps = expand(
            SCALE_ROOT + "/windows/{{sample}}/{{sample}}__{chrom}_gaps.tsv",
            chrom=CHROMS, window_bins="{{window_bins}}"
        ),
        oe_files = lambda w: [
            "{}/oe/{}/{}__{}_oe.npz".format(OUTDIR, w.sample, w.sample, chrom)
            for chrom in CHROMS
        ],
        script = SCRIPT_DIR + "/scripts/sliding_window.py"
    output:
        manifest = SCALE_ROOT + "/logs/windows/{sample}.done"
    log:
        SCALE_ROOT + "/logs/windows/{sample}.log"
    benchmark:
        SCALE_ROOT + "/benchmark/windows/{sample}.tsv"
    threads:
        config["parallel"]["threads_per_job"]
    params:
        chroms = CHROMS,
        resolution = RESOLUTION,
        window_size = lambda w: get_window_scale(w.window_bins).window_size_bp,
        step_size = lambda w: get_window_scale(w.window_bins).step_size_bp,
        sigma = config["smoothing"]["sigma"],
        scale_root = lambda w: SCALE_ROOT.format(window_bins=w.window_bins),
        shared_outdir = OUTDIR
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"

        mkdir -p "{params.scale_root}/windows/${{SAMPLE}}" "{params.scale_root}/logs/windows" \
            "{params.scale_root}/benchmark/windows"

        GAPS=""
        for chrom in {params.chroms}; do
            BASE="${{SAMPLE}}__${{chrom}}"
            GAPS="$GAPS {params.scale_root}/windows/${{SAMPLE}}/${{BASE}}_gaps.tsv"
        done

        OE=""
        for chrom in {params.chroms}; do
            BASE="${{SAMPLE}}__${{chrom}}"
            OE="$OE {params.shared_outdir}/oe/${{SAMPLE}}/${{BASE}}_oe.npz"
        done

        python {input.script} \
            --chroms {params.chroms} \
            --sample "$SAMPLE" \
            --gap-tsv-files $GAPS \
            --ln-oe-files $OE \
            --outdir "{params.scale_root}/windows/${{SAMPLE}}" \
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
        lambda w: checkpoints.sliding_window.get(
            sample=w.sample,
            window_bins=w.window_bins,
        ).output.manifest
    output:
        h5_files = expand(
            SCALE_ROOT + "/windows/{{sample}}/{{sample}}__{{chrom}}_{{stage}}.h5",
            chrom=CHROMS,
            stage=["ln_oe", "gap_filled", "zscore", "gaussian"],
            window_bins="{{window_bins}}"
        )
