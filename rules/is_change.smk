# IS-change candidate construction, locus calling, and diagnostic plots.

rule build_is_change_candidates:
    input:
        knn_enriched=OUTDIR + "/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
        meltron_merged=OUTDIR + "/meltron/{comparison}/{comparison}_meltron_merged.tsv",
    output:
        regions=OUTDIR + "/is_change/{comparison}/{comparison}_candidate_regions.tsv",
    log:
        "{}/logs/is_change/{{comparison}}/candidates.log".format(OUTDIR)
    params:
        outdir=OUTDIR,
        script=SCRIPT_DIR + "/scripts/build_is_change_candidates.py",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}/is_change/{wildcards.comparison}" "{params.outdir}/logs/is_change/{wildcards.comparison}"
        python "{params.script}" \
            --knn-enriched "{input.knn_enriched}" \
            --meltron-merged "{input.meltron_merged}" \
            --config config/config.yaml \
            --out "{output.regions}" \
            --log "{log}" \
            >> "{log}" 2>&1
        """


rule call_is_change_loci:
    input:
        regions=OUTDIR + "/is_change/{comparison}/{comparison}_candidate_regions.tsv",
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
        primary=OUTDIR + "/is_change/{comparison}/{comparison}_is_change_loci.tsv",
        preliminary=OUTDIR + "/is_change/{comparison}/{comparison}_is_change_loci_preliminary.tsv",
        masked=OUTDIR + "/is_change/{comparison}/{comparison}_is_change_masked.bed",
    log:
        "{}/logs/is_change/{{comparison}}/caller.log".format(OUTDIR)
    params:
        outdir=OUTDIR,
        script=SCRIPT_DIR + "/scripts/call_is_change_loci.py",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}/is_change/{wildcards.comparison}" "{params.outdir}/logs/is_change/{wildcards.comparison}"
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
        primary=OUTDIR + "/is_change/{comparison}/{comparison}_is_change_loci.tsv",
        regions=OUTDIR + "/is_change/{comparison}/{comparison}_candidate_regions.tsv",
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
        manifest=OUTDIR + "/figures/is_change/{comparison}/manifest.tsv",
    log:
        "{}/logs/is_change/{{comparison}}/plot.log".format(OUTDIR)
    params:
        outdir=OUTDIR,
        script=SCRIPT_DIR + "/scripts/plot_is_change_loci.py",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}/figures/is_change/{wildcards.comparison}" "{params.outdir}/logs/is_change/{wildcards.comparison}"
        python "{params.script}" \
            --primary "{input.primary}" \
            --candidate-regions "{input.regions}" \
            --is-a "{input.is_a}" \
            --is-b "{input.is_b}" \
            --config config/config.yaml \
            --outdir "{params.outdir}/figures/is_change/{wildcards.comparison}" \
            --cool-a {input.cool_a} \
            --cool-b {input.cool_b} \
            --manifest "{output.manifest}" \
            >> "{log}" 2>&1
        """
