# Convert .cool files to .mcool via cooler zoomify.
#
# Architecture: checkpoint + per-sample wildcard rule.
#
# Step 1 — checkpoint detect_cool_samples:
#   Reads data/manifest.tsv, outputs detect_manifest.tsv with action per sample.
#   Wildcard {sample} is expanded by referencing the checkpoint output.
#
# Step 2 — rule convert_single_cool (per sample):
#   Wildcards: {sample}
#   input: references checkpoint output → triggers checkpoint first, gets detect_manifest.tsv
#   Reads action column from detect_manifest.tsv:
#     SKIP   -> touch done marker, do nothing
#     CONVERT -> cooler zoomify -> results/mcool/{sample}.mcool
#
# Step 3 — rule generate_mcool_manifest (no wildcards):
#   Depends on all convert_single_cool done markers.
#   Reads detect_manifest.tsv, writes data/manifest_mcool.tsv.
#
# Downstream rules read data/manifest_mcool.tsv via get_sample_mcool().

MCOOL_OUTDIR = OUTDIR + "/mcool"
COOL_DETECT_DIR = OUTDIR + "/logs/convert_to_mcool"


checkpoint detect_cool_samples:
    input:
        manifest = config["manifest"],
    output:
        manifest = COOL_DETECT_DIR + "/detect_manifest.tsv",
    log:
        OUTDIR + "/logs/convert_to_mcool_detect.log"
    benchmark:
        "{}/benchmark/convert_to_mcool/detect.tsv".format(OUTDIR)
    threads:
        config["parallel"]["threads_per_job"]
    params:
        detect_dir = COOL_DETECT_DIR,
        script = SCRIPT_DIR + "/scripts/detect_cool_samples.py",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        python "{params.script}" \
            --input "{input.manifest}" \
            --output "{output.manifest}" \
            --log "{log}"
        """


rule convert_single_cool:
    input:
        detect_manifest = lambda w: checkpoints.detect_cool_samples.get().output.manifest,
    output:
        touch = COOL_DETECT_DIR + "/{sample}.converted",
    log:
        COOL_DETECT_DIR + "/{sample}.log"
    benchmark:
        "{}/benchmark/convert_to_mcool/{{sample}}.tsv".format(OUTDIR)
    threads:
        config["convert"]["threads"]
    params:
        outdir = MCOOL_OUTDIR,
        resolutions = ",".join(map(str, config["convert"]["resolutions"])),
        script = SCRIPT_DIR + "/scripts/convert_cool_to_mcool.py",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        python "{params.script}" \
            --sample "{wildcards.sample}" \
            --detect-manifest "{input.detect_manifest}" \
            --outdir "{params.outdir}" \
            --resolutions "{params.resolutions}" \
            --threads {threads} \
            --log-file "{log}" \
            --touch-file "{output.touch}"
        """


rule generate_mcool_manifest:
    input:
        detect_manifest = COOL_DETECT_DIR + "/detect_manifest.tsv",
        converted = expand(COOL_DETECT_DIR + "/{sample}.converted", sample=SAMPLES),
    output:
        manifest = "data/manifest_mcool.tsv",
    log:
        OUTDIR + "/logs/generate_mcool_manifest.log"
    benchmark:
        "{}/benchmark/convert_to_mcool/generate_manifest.tsv".format(OUTDIR)
    threads:
        config["parallel"]["threads_per_job"]
    params:
        outdir = MCOOL_OUTDIR,
        script = SCRIPT_DIR + "/scripts/generate_mcool_manifest.py",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        python "{params.script}" \
            --detect-manifest "{input.detect_manifest}" \
            --outdir "{params.outdir}" \
            --output "{output.manifest}"
        """
