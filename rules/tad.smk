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
        resolution=lambda wc: TAD_RESOLUTION_BY_LABEL[wc.resolution_label],
        run_dir=lambda wc: tad_run_dir(
            wc.sample, TAD_RESOLUTION_BY_LABEL[wc.resolution_label], wc.signal
        ),
        d_config_json=TAD_D_CONFIG_JSON,
        excluded_chromosomes=",".join(TAD.excluded_chromosomes),
        bad_region_bed=TAD.bad_region_bed or "None",
        badregionfilter=str(TAD.badregionfilter).lower(),
    conda:
        "../env/tad.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.ready})" "$(dirname {log})" "$(dirname {output.ready})/../benchmark/tad_prepare"
        python scripts/tad_from_mcool.py \
            --mcool {input.mcool:q} \
            --sample {wildcards.sample:q} \
            --signal {wildcards.signal:q} \
            --resolution {params.resolution} \
            --run-dir {params.run_dir:q} \
            --ready-flag {output.ready:q} \
            --processors {threads} \
            --badregionfilter {params.badregionfilter:q} \
            --bad-region-bed {params.bad_region_bed:q} \
            --excluded-chromosomes {params.excluded_chromosomes:q} \
            --d-config-json {params.d_config_json:q} > {log:q} 2>&1
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
        resolution=lambda wc: TAD_RESOLUTION_BY_LABEL[wc.resolution_label],
        settings=lambda wc: (
            f"{tad_run_dir(wc.sample, TAD_RESOLUTION_BY_LABEL[wc.resolution_label], wc.signal)}"
            f"/settings_{wc.sample}.{wc.signal}.txt"
        ),
        code_dir=TAD.code_dir,
    conda:
        "../env/tad.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bed})" "$(dirname {log})" "$(dirname {output.bed})/../benchmark/tad_call"
        python scripts/run_3dnetmod.py \
            --settings {params.settings:q} \
            --code-dir {params.code_dir:q} \
            --mcool {input.mcool:q} \
            --sample {wildcards.sample:q} \
            --signal {wildcards.signal:q} \
            --resolution {params.resolution} \
            --bed-output {output.bed:q} \
            --bedpe-output {output.bedpe:q} > {log:q} 2>&1
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
        resolution=lambda wc: TAD_RESOLUTION_BY_LABEL[wc.resolution_label],
    conda:
        "../env/tad.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bed})" "$(dirname {log})" "$(dirname {output.bed})/../benchmark/tad_boundaries"
        python scripts/tad_boundaries.py extract \
            --tad {input.tad:q} --mcool {input.mcool:q} \
            --resolution {params.resolution} \
            --sample {wildcards.sample:q} --signal {wildcards.signal:q} \
            --output {output.bed:q} > {log:q} 2>&1
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
    conda:
        "../env/tad.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bed})" "$(dirname {log})" "$(dirname {output.bed})/../benchmark/tad_boundaries"
        python scripts/tad_boundaries.py merge \
            --inputs {input.boundaries:q} --mcool {input.mcool:q} \
            --chrom-resolution 4000 \
            --sample {wildcards.sample:q} --signal {wildcards.signal:q} \
            --output {output.bed:q} > {log:q} 2>&1
        """
