# Validate input mcool files and extract chromosome list

rule validate_inputs:
    input:
        manifest = "data/manifest_mcool.tsv",
        script = SCRIPT_DIR + "/scripts/validate_inputs.py",
    output:
        validated = OUTDIR + "/validated/manifest.validated.tsv"
    log:
        OUTDIR + "/logs/validate_inputs.log"
    benchmark:
        "{}/benchmark/validate_inputs/validate_inputs.tsv".format(OUTDIR)
    threads:
        config["parallel"]["threads_per_job"]
    params:
        outdir = OUTDIR
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.outdir}/validated" "{params.outdir}/logs"

        INPUT="{input.manifest}"
        OUTPUT="{output.validated}"
        LOG="{log}"

        python {input.script} \
            --input "$INPUT" \
            --output "$OUTPUT" \
            --log "$LOG" \
            > "$LOG" 2>&1
        """
