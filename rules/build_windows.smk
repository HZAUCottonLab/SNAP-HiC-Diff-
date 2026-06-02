# Build genomic windows for analysis

rule build_windows:
    input:
        validated = OUTDIR + "/validated/manifest.validated.tsv",
        script = SCRIPT_DIR + "/scripts/build_windows.py"
    output:
        windows = "{}/windows/windows_{}.parquet".format(OUTDIR, RESOLUTION)
    log:
        OUTDIR + "/logs/build_windows.log"
    benchmark:
        "{}/benchmark/build_windows/build_windows.tsv".format(OUTDIR)
    threads:
        config["parallel"]["threads_per_job"]
    params:
        resolution = RESOLUTION,
        window_size = config["window_size"],
        step_size = config["step_size"],
        outdir = OUTDIR
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.outdir}/windows" "{params.outdir}/logs"

        python {input.script} \
            --input "{input.validated}" \
            --output "{output.windows}" \
            --log "{log}" \
            --resolution {params.resolution} \
            --window-size {params.window_size} \
            --step-size {params.step_size} \
            > "{log}" 2>&1
        """
