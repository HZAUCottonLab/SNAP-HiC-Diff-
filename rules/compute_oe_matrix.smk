# Compute O/E matrix from balanced per-chromosome .cool using pre-computed expected values
# One job per sample, loops over all chroms.
#
# Inputs: results/balanced/{sample}/{sample}__{chrom}_iced.cool
# Inputs: results/expected/{sample}/{sample}__{chrom}_expected.npz
# Outputs: results/oe/{sample}/{sample}__{chrom}_oe.npz

rule compute_oe_matrix:
    input:
        lambda wildcards: [
            "{}/logs/expected/{}_to_npz.done".format(OUTDIR, wildcards.sample)
        ],
        script = SCRIPT_DIR + "/scripts/compute_oe_matrix_batch.py"
    output:
        manifest = "{}/logs/oe/{{sample}}.done".format(OUTDIR),
        npz_files = expand(
            OUTDIR + "/oe/{{sample}}/{{sample}}__{chrom}_oe.npz",
            chrom=CHROMS
        ),
    log:
        "{}/logs/oe/{{sample}}_oe.log".format(OUTDIR)
    benchmark:
        "{}/benchmark/oe/{{sample}}_oe.tsv".format(OUTDIR)
    threads:
        config["parallel"]["threads_per_job"]
    params:
        chroms = CHROMS,
        outdir = OUTDIR
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"

        mkdir -p "{params.outdir}/oe" "{params.outdir}/logs/oe"

        INPUTS=""
        EXPECTED=""
        for chrom in {params.chroms}; do
            BASE="${{SAMPLE}}__${{chrom}}"
            INPUTS="$INPUTS {params.outdir}/balanced/${{SAMPLE}}/${{BASE}}_iced.cool"
            EXPECTED="$EXPECTED {params.outdir}/expected/${{SAMPLE}}/${{BASE}}_expected.npz"
        done

        python {input.script} \
            --chroms {params.chroms} \
            --sample "$SAMPLE" \
            --outdir {params.outdir} \
            --inputs $INPUTS \
            --expected $EXPECTED \
            --log "{log}" \
            > "{log}" 2>&1

        touch "{output.manifest}"
        """
