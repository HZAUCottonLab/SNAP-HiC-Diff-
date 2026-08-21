# Build genomic windows for analysis

rule build_windows:
    input:
        validated = OUTDIR + "/validated/manifest.validated.tsv",
        script = SCRIPT_DIR + "/scripts/build_windows.py"
    output:
        windows = SCALE_ROOT + "/windows/windows.parquet"
    log:
        SCALE_ROOT + "/logs/build_windows.log"
    benchmark:
        SCALE_ROOT + "/benchmark/build_windows/build_windows.tsv"
    threads:
        config["parallel"]["threads_per_job"]
    params:
        resolution = RESOLUTION,
        window_bins = lambda w: get_window_scale(w.window_bins).window_bins,
        window_size = lambda w: get_window_scale(w.window_bins).window_size_bp,
        step_bins = lambda w: get_window_scale(w.window_bins).step_bins,
        step_size = lambda w: get_window_scale(w.window_bins).step_size_bp,
        scale_root = lambda w: SCALE_ROOT.format(window_bins=w.window_bins)
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.scale_root}/windows" "{params.scale_root}/logs" \
            "{params.scale_root}/benchmark/build_windows"

        python {input.script} \
            --input "{input.validated}" \
            --output "{output.windows}" \
            --log "{log}" \
            --resolution {params.resolution} \
            --window-bins {params.window_bins} \
            --window-size {params.window_size} \
            --step-bins {params.step_bins} \
            --step-size {params.step_size} \
            > "{log}" 2>&1
        """
