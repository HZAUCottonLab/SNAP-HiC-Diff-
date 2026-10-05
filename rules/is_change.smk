# Aggregate all scale evidence once before calling and plotting IS-change loci.

rule build_is_change_candidates:
    input:
        knn_enriched=lambda w: expand(
            SCALE_ROOT + "/{scale}/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
            scale=SCALE_LABELS,
            comparison=w.comparison,
        ),
        meltron_merged=lambda w: expand(
            SCALE_ROOT + "/{scale}/meltron/{comparison}/{comparison}_meltron_merged.tsv",
            scale=SCALE_LABELS,
            comparison=w.comparison,
        ),
    output:
        candidates=MULTISCALE_ROOT + "/{comparison}/all_scale_candidates.tsv",
        union_tsv=MULTISCALE_ROOT + "/{comparison}/candidate_union.tsv",
        union_bed=MULTISCALE_ROOT + "/{comparison}/candidate_union.bed",
    log:
        MULTISCALE_ROOT + "/{comparison}/logs/is_change/candidates.log",
    params:
        script=SCRIPT_DIR + "/scripts/build_is_change_candidates.py",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.candidates}")" "$(dirname "{log}")"
        python "{params.script}" \
            --knn-enriched {input.knn_enriched} \
            --meltron-merged {input.meltron_merged} \
            --config config/config.yaml \
            --out-candidates "{output.candidates}" \
            --out-union-tsv "{output.union_tsv}" \
            --out-union-bed "{output.union_bed}" \
            --log "{log}" \
            >> "{log}" 2>&1
        """


rule call_is_change_loci:
    input:
        regions=MULTISCALE_ROOT + "/{comparison}/candidate_union.tsv",
        is_a=lambda w: get_comparison_is_file(w.comparison, 0),
        is_b=lambda w: get_comparison_is_file(w.comparison, 1),
        cool_a=lambda w: expand(
            OUTDIR + "/balanced/{sample}/{sample}__{chrom}_iced.cool",
            sample=get_pair_samples(w.comparison)[0],
            chrom=CHROMS,
        ),
        cool_b=lambda w: expand(
            OUTDIR + "/balanced/{sample}/{sample}__{chrom}_iced.cool",
            sample=get_pair_samples(w.comparison)[1],
            chrom=CHROMS,
        ),
    output:
        primary=MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci.tsv",
        preliminary=MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci_preliminary.tsv",
        masked=MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_masked.bed",
    log:
        MULTISCALE_ROOT + "/{comparison}/logs/is_change/caller.log",
    params:
        script=SCRIPT_DIR + "/scripts/call_is_change_loci.py",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.primary}")" "$(dirname "{log}")"
        python "{params.script}" \
            --candidate-regions "{input.regions}" \
            --is-a "{input.is_a}" \
            --is-b "{input.is_b}" \
            --cool-a {input.cool_a} \
            --cool-b {input.cool_b} \
            --config config/config.yaml \
            --out-primary "{output.primary}" \
            --out-preliminary "{output.preliminary}" \
            --out-masked "{output.masked}" \
            --log "{log}" \
            >> "{log}" 2>&1
        """


rule plot_is_change_loci:
    input:
        primary=lambda w: (
            MULTISCALE_ROOT + "/{}/is_change/{}_is_change_loci.boundary_annotated.tsv".format(
                w.comparison, w.comparison
            )
            if IS_BOUNDARY.enabled
            else MULTISCALE_ROOT + "/{}/is_change/{}_is_change_loci.tsv".format(
                w.comparison, w.comparison
            )
        ),
        regions=MULTISCALE_ROOT + "/{comparison}/candidate_union.tsv",
        boundary_events=lambda w: (
            boundary_event_path(w.comparison) if IS_BOUNDARY.enabled else []
        ),
        is_a=lambda w: get_comparison_is_file(w.comparison, 0),
        is_b=lambda w: get_comparison_is_file(w.comparison, 1),
        cool_a=lambda w: expand(
            OUTDIR + "/balanced/{sample}/{sample}__{chrom}_iced.cool",
            sample=get_pair_samples(w.comparison)[0],
            chrom=CHROMS,
        ),
        cool_b=lambda w: expand(
            OUTDIR + "/balanced/{sample}/{sample}__{chrom}_iced.cool",
            sample=get_pair_samples(w.comparison)[1],
            chrom=CHROMS,
        ),
    output:
        manifest=MULTISCALE_ROOT + "/{comparison}/is_change/figures/manifest.tsv",
    log:
        MULTISCALE_ROOT + "/{comparison}/logs/is_change/plot.log",
    params:
        script=SCRIPT_DIR + "/scripts/plot_is_change_loci.py",
        boundary_events=lambda w: (
            "--boundary-events \"{}\"".format(boundary_event_path(w.comparison))
            if IS_BOUNDARY.enabled
            else ""
        ),
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.manifest}")" "$(dirname "{log}")"
        python "{params.script}" \
            --primary "{input.primary}" \
            --candidate-regions "{input.regions}" \
            --is-a "{input.is_a}" \
            --is-b "{input.is_b}" \
            --config config/config.yaml \
            --outdir "$(dirname "{output.manifest}")" \
            --cool-a {input.cool_a} \
            --cool-b {input.cool_b} \
            --manifest "{output.manifest}" \
            {params.boundary_events} \
            >> "{log}" 2>&1
        """
