# Cross-scale candidate construction, locus calling, and diagnostic plots.

rule build_is_change_candidates:
    input:
        knn_enriched=lambda w: expand(
            SCALE_ROOT + "/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
            window_bins=WINDOW_BINS,
            comparison=w.comparison,
        ),
        meltron_merged=lambda w: expand(
            SCALE_ROOT + "/meltron/{comparison}/{comparison}_meltron_merged.tsv",
            window_bins=WINDOW_BINS,
            comparison=w.comparison,
        ),
    output:
        all_candidates=MULTISCALE_ROOT + "/{comparison}/all_scale_candidates.tsv",
        regions=MULTISCALE_ROOT + "/{comparison}/candidate_union.tsv",
        bed=MULTISCALE_ROOT + "/{comparison}/candidate_union.bed",
    log:
        MULTISCALE_ROOT + "/{comparison}/is_change/logs/candidates.log"
    params:
        expected_window_bins=" ".join(map(str, WINDOW_BINS)),
        script=SCRIPT_DIR + "/scripts/build_is_change_candidates.py",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.all_candidates}")" "$(dirname "{log}")"
        python "{params.script}" \
            --knn-enriched {input.knn_enriched} \
            --meltron-merged {input.meltron_merged} \
            --expected-window-bins {params.expected_window_bins} \
            --config config/config.yaml \
            --out-all "{output.all_candidates}" \
            --out-regions "{output.regions}" \
            --out-bed "{output.bed}" \
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
        MULTISCALE_ROOT + "/{comparison}/is_change/logs/caller.log"
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
        primary=MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci.tsv",
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
        manifest=MULTISCALE_ROOT + "/{comparison}/is_change/figures/manifest.tsv",
    log:
        MULTISCALE_ROOT + "/{comparison}/is_change/logs/plot.log"
    params:
        figure_dir=MULTISCALE_ROOT + "/{comparison}/is_change/figures",
        script=SCRIPT_DIR + "/scripts/plot_is_change_loci.py",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.figure_dir}" "$(dirname "{log}")"
        python "{params.script}" \
            --primary "{input.primary}" \
            --candidate-regions "{input.regions}" \
            --is-a "{input.is_a}" \
            --is-b "{input.is_b}" \
            --config config/config.yaml \
            --outdir "{params.figure_dir}" \
            --cool-a {input.cool_a} \
            --cool-b {input.cool_b} \
            --manifest "{output.manifest}" \
            >> "{log}" 2>&1
        """
