# Compute insulation scores (cooltools insulation) and convert to Meltron IS format.
# One job per sample, loops over all chromosomes.
#
# Step 1: cooltools insulation per chromosome (all distances in one call per chrom)
# Step 2: convert cooltools output to Meltron IS format (4-col: chrom, viewpoint, IS_distance, IS)
#
# Inputs: results/balanced/{sample}/{sample}__{chrom}_iced.cool
# Outputs: results/IS/{sample}/{sample}_IS.tsv (Meltron format)

checkpoint compute_IS:
    input:
        lambda wildcards: checkpoints.extract_and_balance_sample.get(
            sample=wildcards.sample
        ).output.manifest
    output:
        touch = OUTDIR + "/logs/compute_IS/{sample}.done",
        is_file = OUTDIR + "/IS/{sample}/{sample}_IS.tsv"
    log:
        "{}/logs/compute_IS/{{sample}}_compute_IS.log".format(OUTDIR)
    benchmark:
        "{}/benchmark/compute_IS/{{sample}}_compute_IS.tsv".format(OUTDIR)
    threads:
        config["parallel"]["threads_per_job"]
    params:
        chroms = " ".join(CHROMS),
        is_distances_str = ",".join(str(d) for d in config["meltron"]["is_distances"]),
        outdir = OUTDIR,
        tmpdir = TMP,
        script_batch = SCRIPT_DIR + "/scripts/compute_IS_batch.py",
        script_convert = SCRIPT_DIR + "/scripts/convert_cooltools_IS_to_meltron_IS.py"
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"

        mkdir -p "{params.outdir}/insulation/${{SAMPLE}}" \
                 "{params.outdir}/IS/${{SAMPLE}}" \
                 "{params.outdir}/logs/compute_IS"

        # Step 1: cooltools insulation per chromosome (all distances in one command)
        python {params.script_batch} \
            --sample "$SAMPLE" \
            --chroms "{params.chroms}" \
            --outdir "{params.outdir}" \
            --is-distances "{params.is_distances_str}" \
            --tmpdir "{params.tmpdir}" \
            --n-jobs {threads} \
            --log-file "{log}" \
            >> "{log}" 2>&1

        # Step 2: convert cooltools TSVs to Meltron IS format, merge all chroms
        python {params.script_convert} \
            --sample "$SAMPLE" \
            --chroms "{params.chroms}" \
            --outdir "{params.outdir}" \
            --is-distances "{params.is_distances_str}" \
            --log-file "{log}" \
            >> "{log}" 2>&1

        touch "{output.touch}"
        """
