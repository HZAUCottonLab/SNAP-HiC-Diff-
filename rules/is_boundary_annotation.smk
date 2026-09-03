rule unify_tad_boundaries:
    input:
        raw_multires=lambda w: tad_merged_boundary_path(w.sample, "raw"),
        balanced_multires=lambda w: tad_merged_boundary_path(w.sample, "balanced"),
        raw_resolution=lambda w: [
            tad_boundary_path(w.sample, resolution, "raw")
            for resolution in TAD_RESOLUTIONS
        ],
        balanced_resolution=lambda w: [
            tad_boundary_path(w.sample, resolution, "balanced")
            for resolution in TAD_RESOLUTIONS
        ],
        mcool=lambda w: get_sample_mcool(w.sample),
    output:
        unified=OUTDIR + "/tad_boundaries/{sample}.multires.unified.boundaries.tsv",
        offsets=OUTDIR + "/tad_boundaries/{sample}.boundary_jitter_offsets.tsv",
        views=OUTDIR + "/tad_boundaries/{sample}.annotation_views.tsv",
    log:
        OUTDIR + "/logs/is_boundary/{sample}.unify.log",
    benchmark:
        OUTDIR + "/benchmark/is_boundary/{sample}.unify.tsv"
    threads: 1
    params:
        script=SCRIPT_DIR + "/scripts/unify_tad_boundaries.py",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.unified}")" "$(dirname "{log}")"
        python "{params.script}" \
            --sample "{wildcards.sample}" \
            --raw-multires "{input.raw_multires}" \
            --balanced-multires "{input.balanced_multires}" \
            --raw-resolution-boundaries {input.raw_resolution} \
            --balanced-resolution-boundaries {input.balanced_resolution} \
            --mcool "{input.mcool}" \
            --config config/config.yaml \
            --output-unified "{output.unified}" \
            --output-offsets "{output.offsets}" \
            --output-views "{output.views}" \
            --log "{log}"
        """


rule call_tad_boundary_events:
    input:
        boundaries_a=lambda w: unified_boundary_path(get_pair_samples(w.comparison)[0]),
        boundaries_b=lambda w: unified_boundary_path(get_pair_samples(w.comparison)[1]),
        offsets_a=lambda w: boundary_jitter_path(get_pair_samples(w.comparison)[0]),
        offsets_b=lambda w: boundary_jitter_path(get_pair_samples(w.comparison)[1]),
        views_a=lambda w: boundary_view_path(get_pair_samples(w.comparison)[0]),
        views_b=lambda w: boundary_view_path(get_pair_samples(w.comparison)[1]),
    output:
        events=OUTDIR + "/tad_boundary_events/{comparison}.boundary_events.tsv",
        stability=OUTDIR + "/tad_boundary_events/{comparison}.boundary_event_stability.tsv",
        qc=OUTDIR + "/tad_boundary_events/{comparison}.matching_qc.tsv",
    log:
        OUTDIR + "/logs/is_boundary/{comparison}.matching.log",
    benchmark:
        OUTDIR + "/benchmark/is_boundary/{comparison}.matching.tsv"
    threads: 1
    params:
        script=SCRIPT_DIR + "/scripts/call_tad_boundary_events.py",
        sample_a=lambda w: get_pair_samples(w.comparison)[0],
        sample_b=lambda w: get_pair_samples(w.comparison)[1],
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.events}")" "$(dirname "{log}")"
        python "{params.script}" \
            --comparison "{wildcards.comparison}" \
            --sample-a "{params.sample_a}" \
            --sample-b "{params.sample_b}" \
            --boundaries-a "{input.boundaries_a}" \
            --boundaries-b "{input.boundaries_b}" \
            --offsets-a "{input.offsets_a}" \
            --offsets-b "{input.offsets_b}" \
            --views-a "{input.views_a}" \
            --views-b "{input.views_b}" \
            --config config/config.yaml \
            --output-events "{output.events}" \
            --output-stability "{output.stability}" \
            --output-qc "{output.qc}" \
            --log "{log}"
        """


rule annotate_is_change_boundaries:
    input:
        events=OUTDIR + "/tad_boundary_events/{comparison}.boundary_events.tsv",
        primary=MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci.tsv",
        preliminary=MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci_preliminary.tsv",
        masked=MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_masked.bed",
        is_a=lambda w: get_comparison_is_file(w.comparison, 0),
        is_b=lambda w: get_comparison_is_file(w.comparison, 1),
    output:
        primary=MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci.boundary_annotated.tsv",
        preliminary=MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci_preliminary.boundary_annotated.tsv",
        strength=OUTDIR + "/tad_boundary_events/{comparison}.boundary_strength_by_scale.tsv",
        strength_summary=OUTDIR + "/tad_boundary_events/{comparison}.boundary_strength_summary.tsv",
        qc=OUTDIR + "/tad_boundary_events/{comparison}.annotation_qc.tsv",
        summary=OUTDIR + "/tad_boundary_events/{comparison}.annotation_summary.tsv",
    log:
        OUTDIR + "/logs/is_boundary/{comparison}.annotation.log",
    benchmark:
        OUTDIR + "/benchmark/is_boundary/{comparison}.annotation.tsv"
    threads: 1
    params:
        script=SCRIPT_DIR + "/scripts/annotate_is_change_boundaries.py",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.primary}")" "$(dirname "{log}")"
        python "{params.script}" \
            --comparison "{wildcards.comparison}" \
            --events "{input.events}" \
            --primary "{input.primary}" \
            --preliminary "{input.preliminary}" \
            --masked "{input.masked}" \
            --is-a "{input.is_a}" \
            --is-b "{input.is_b}" \
            --config config/config.yaml \
            --out-primary "{output.primary}" \
            --out-preliminary "{output.preliminary}" \
            --out-strength "{output.strength}" \
            --out-strength-summary "{output.strength_summary}" \
            --out-qc "{output.qc}" \
            --out-summary "{output.summary}" \
            --log "{log}"
        """


rule permute_is_boundary_association:
    input:
        loci=MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci.boundary_annotated.tsv",
        regions=MULTISCALE_ROOT + "/{comparison}/candidate_union.tsv",
        events=OUTDIR + "/tad_boundary_events/{comparison}.boundary_events.tsv",
        masked=MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_masked.bed",
    output:
        permutation=OUTDIR + "/tad_boundary_events/{comparison}.permutation.tsv",
    log:
        OUTDIR + "/logs/is_boundary/{comparison}.permutation.log",
    benchmark:
        OUTDIR + "/benchmark/is_boundary/{comparison}.permutation.tsv"
    threads: 1
    params:
        script=SCRIPT_DIR + "/scripts/permute_is_boundary_association.py",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.permutation}")" "$(dirname "{log}")"
        python "{params.script}" \
            --loci "{input.loci}" \
            --candidate-regions "{input.regions}" \
            --events "{input.events}" \
            --masked "{input.masked}" \
            --config config/config.yaml \
            --output "{output.permutation}" \
            --log "{log}"
        """
