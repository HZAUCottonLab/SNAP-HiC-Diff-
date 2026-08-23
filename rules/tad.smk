wildcard_constraints:
    resolution_label="|".join(TAD_RESOLUTION_BY_LABEL),
    signal="|".join(TAD_SIGNALS)


rule prepare_tad:
    input:
        manifest="data/manifest_mcool.tsv",
        mcool=lambda wc: get_sample_mcool(wc.sample),
    output:
        ready=OUTDIR + "/tad_prepare/{sample}.{resolution_label}.{signal}.ready",
    log:
        OUTDIR + "/logs/tad_prepare/{sample}.{resolution_label}.{signal}.log"
    benchmark:
        OUTDIR + "/benchmark/tad_prepare/{sample}.{resolution_label}.{signal}.tsv"
    threads:
        TAD.threads
    params:
        resolution_map=TAD_RESOLUTION_MAP,
        run_root=TAD.run_root,
        d_config_json=TAD_D_CONFIG_JSON,
        excluded_chromosomes=",".join(TAD.excluded_chromosomes),
        bad_region_bed=TAD.bad_region_bed or "None",
        badregionfilter=str(TAD.badregionfilter).lower(),
    conda:
        "../env/tad.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"
        RESOLUTION_LABEL="{wildcards.resolution_label}"
        SIGNAL="{wildcards.signal}"
        MCOOL={input.mcool:q}
        READY_FLAG={output.ready:q}
        LOG_FILE={log:q}
        D_CONFIG_JSON={params.d_config_json:q}
        EXCLUDED_CHROMOSOMES={params.excluded_chromosomes:q}
        BAD_REGION_BED={params.bad_region_bed:q}
        BADREGIONFILTER={params.badregionfilter:q}
        TAD_RESOLUTION_MAP={params.resolution_map:q}
        RUN_ROOT={params.run_root:q}
        THREADS="{threads}"
        RESOLUTION="$(printf '%s\n' "$TAD_RESOLUTION_MAP" | tr ',' '\n' | grep -F "$RESOLUTION_LABEL=" | cut -d= -f2 || true)"
        if [ -z "$RESOLUTION" ]; then
            printf 'ERROR: unknown TAD resolution label %s\n' "$RESOLUTION_LABEL" >&2
            exit 1
        fi
        RUN_DIR="$RUN_ROOT/$SAMPLE/$SIGNAL/$RESOLUTION_LABEL/D"
        READY_DIR="$(dirname "$READY_FLAG")"
        LOG_DIR="$(dirname "$LOG_FILE")"
        mkdir -p "$READY_DIR" "$LOG_DIR" "$READY_DIR/../benchmark/tad_prepare"
        python scripts/tad_from_mcool.py \
            --mcool "$MCOOL" \
            --sample "$SAMPLE" \
            --signal "$SIGNAL" \
            --resolution "$RESOLUTION" \
            --run-dir "$RUN_DIR" \
            --ready-flag "$READY_FLAG" \
            --processors "$THREADS" \
            --badregionfilter "$BADREGIONFILTER" \
            --bad-region-bed "$BAD_REGION_BED" \
            --excluded-chromosomes "$EXCLUDED_CHROMOSOMES" \
            --d-config-json "$D_CONFIG_JSON" > "$LOG_FILE" 2>&1
        """


rule call_tad:
    input:
        ready=lambda wc: tad_ready_path(
            wc.sample, TAD_RESOLUTION_BY_LABEL[wc.resolution_label], wc.signal
        ),
        mcool=lambda wc: get_sample_mcool(wc.sample),
    output:
        bed=OUTDIR + "/tads/{sample}.{resolution_label}.{signal}.3dnetmod.bed",
        bedpe=OUTDIR + "/tads/{sample}.{resolution_label}.{signal}.3dnetmod.bedpe",
    log:
        OUTDIR + "/logs/tad_call/{sample}.{resolution_label}.{signal}.log"
    benchmark:
        OUTDIR + "/benchmark/tad_call/{sample}.{resolution_label}.{signal}.tsv"
    threads:
        TAD.threads
    params:
        resolution_map=TAD_RESOLUTION_MAP,
        run_root=TAD.run_root,
        code_dir=TAD.code_dir,
    conda:
        "../env/tad.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"
        RESOLUTION_LABEL="{wildcards.resolution_label}"
        SIGNAL="{wildcards.signal}"
        CODE_DIR={params.code_dir:q}
        MCOOL={input.mcool:q}
        BED_OUTPUT={output.bed:q}
        BEDPE_OUTPUT={output.bedpe:q}
        LOG_FILE={log:q}
        TAD_RESOLUTION_MAP={params.resolution_map:q}
        RUN_ROOT={params.run_root:q}
        THREADS="{threads}"
        RESOLUTION="$(printf '%s\n' "$TAD_RESOLUTION_MAP" | tr ',' '\n' | grep -F "$RESOLUTION_LABEL=" | cut -d= -f2 || true)"
        if [ -z "$RESOLUTION" ]; then
            printf 'ERROR: unknown TAD resolution label %s\n' "$RESOLUTION_LABEL" >&2
            exit 1
        fi
        RUN_DIR="$RUN_ROOT/$SAMPLE/$SIGNAL/$RESOLUTION_LABEL/D"
        SETTINGS="$RUN_DIR/settings_${{SAMPLE}}.${{SIGNAL}}.txt"
        OUTPUT_DIR="$(dirname "$BED_OUTPUT")"
        LOG_DIR="$(dirname "$LOG_FILE")"
        mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$OUTPUT_DIR/../benchmark/tad_call"
        python scripts/run_3dnetmod.py \
            --settings "$SETTINGS" \
            --code-dir "$CODE_DIR" \
            --mcool "$MCOOL" \
            --sample "$SAMPLE" \
            --signal "$SIGNAL" \
            --resolution "$RESOLUTION" \
            --bed-output "$BED_OUTPUT" \
            --bedpe-output "$BEDPE_OUTPUT" > "$LOG_FILE" 2>&1
        """


rule extract_tad_boundaries:
    input:
        tad=OUTDIR + "/tads/{sample}.{resolution_label}.{signal}.3dnetmod.bed",
        mcool=lambda wc: get_sample_mcool(wc.sample),
    output:
        bed=OUTDIR + "/tad_boundaries/{sample}.{resolution_label}.{signal}.boundaries.bed",
    log:
        OUTDIR + "/logs/tad_boundaries/{sample}.{resolution_label}.{signal}.log"
    benchmark:
        OUTDIR + "/benchmark/tad_boundaries/{sample}.{resolution_label}.{signal}.tsv"
    threads: 1
    params:
        resolution_map=TAD_RESOLUTION_MAP,
    conda:
        "../env/tad.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"
        RESOLUTION_LABEL="{wildcards.resolution_label}"
        SIGNAL="{wildcards.signal}"
        TAD_INPUT={input.tad:q}
        MCOOL={input.mcool:q}
        BED_OUTPUT={output.bed:q}
        LOG_FILE={log:q}
        TAD_RESOLUTION_MAP={params.resolution_map:q}
        RESOLUTION="$(printf '%s\n' "$TAD_RESOLUTION_MAP" | tr ',' '\n' | grep -F "$RESOLUTION_LABEL=" | cut -d= -f2 || true)"
        if [ -z "$RESOLUTION" ]; then
            printf 'ERROR: unknown TAD resolution label %s\n' "$RESOLUTION_LABEL" >&2
            exit 1
        fi
        OUTPUT_DIR="$(dirname "$BED_OUTPUT")"
        LOG_DIR="$(dirname "$LOG_FILE")"
        mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$OUTPUT_DIR/../benchmark/tad_boundaries"
        python scripts/tad_boundaries.py extract \
            --tad "$TAD_INPUT" --mcool "$MCOOL" \
            --resolution "$RESOLUTION" \
            --sample "$SAMPLE" --signal "$SIGNAL" \
            --output "$BED_OUTPUT" > "$LOG_FILE" 2>&1
        """


rule merge_tad_boundaries:
    input:
        boundaries=lambda wc: [
            tad_boundary_path(wc.sample, resolution, wc.signal)
            for resolution in TAD_RESOLUTIONS
        ],
        mcool=lambda wc: get_sample_mcool(wc.sample),
    output:
        bed=OUTDIR + "/tad_boundaries/{sample}.multires.{signal}.boundaries.bed",
    log:
        OUTDIR + "/logs/tad_boundaries/{sample}.multires.{signal}.log"
    benchmark:
        OUTDIR + "/benchmark/tad_boundaries/{sample}.multires.{signal}.tsv"
    threads: 1
    params:
        chrom_resolution=min(TAD_RESOLUTIONS),
    conda:
        "../env/tad.yaml"
    shell:
        r"""
        set -euo pipefail

        SAMPLE="{wildcards.sample}"
        SIGNAL="{wildcards.signal}"
        BOUNDARIES=({input.boundaries:q})
        MCOOL={input.mcool:q}
        BED_OUTPUT={output.bed:q}
        LOG_FILE={log:q}
        CHROM_RESOLUTION="{params.chrom_resolution}"
        OUTPUT_DIR="$(dirname "$BED_OUTPUT")"
        LOG_DIR="$(dirname "$LOG_FILE")"
        mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$OUTPUT_DIR/../benchmark/tad_boundaries"
        python scripts/tad_boundaries.py merge \
            --inputs "${{BOUNDARIES[@]}}" --mcool "$MCOOL" \
            --chrom-resolution "$CHROM_RESOLUTION" \
            --sample "$SAMPLE" --signal "$SIGNAL" \
            --output "$BED_OUTPUT" > "$LOG_FILE" 2>&1
        """
