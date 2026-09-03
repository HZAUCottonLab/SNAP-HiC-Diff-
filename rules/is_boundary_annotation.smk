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
        config="config/config.yaml",
        benchmark_dir=OUTDIR + "/benchmark/is_boundary",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.unified:q})" "$(dirname {log:q})" {params.benchmark_dir:q}
        python {params.script:q} \
            --sample {wildcards.sample:q} \
            --raw-multires {input.raw_multires:q} \
            --balanced-multires {input.balanced_multires:q} \
            --raw-resolution-boundaries {input.raw_resolution:q} \
            --balanced-resolution-boundaries {input.balanced_resolution:q} \
            --mcool {input.mcool:q} \
            --config {params.config:q} \
            --output-unified {output.unified:q} \
            --output-offsets {output.offsets:q} \
            --output-views {output.views:q} \
            --log {log:q} \
            >> {log:q} 2>&1
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
        config="config/config.yaml",
        benchmark_dir=OUTDIR + "/benchmark/is_boundary",
        sample_a=lambda w: get_pair_samples(w.comparison)[0],
        sample_b=lambda w: get_pair_samples(w.comparison)[1],
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.events:q})" "$(dirname {log:q})" {params.benchmark_dir:q}
        python {params.script:q} \
            --comparison {wildcards.comparison:q} \
            --sample-a {params.sample_a:q} \
            --sample-b {params.sample_b:q} \
            --boundaries-a {input.boundaries_a:q} \
            --boundaries-b {input.boundaries_b:q} \
            --offsets-a {input.offsets_a:q} \
            --offsets-b {input.offsets_b:q} \
            --views-a {input.views_a:q} \
            --views-b {input.views_b:q} \
            --config {params.config:q} \
            --output-events {output.events:q} \
            --output-stability {output.stability:q} \
            --output-qc {output.qc:q} \
            --log {log:q} \
            >> {log:q} 2>&1
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
        config="config/config.yaml",
        benchmark_dir=OUTDIR + "/benchmark/is_boundary",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.primary:q})" "$(dirname {log:q})" {params.benchmark_dir:q}
        python {params.script:q} \
            --comparison {wildcards.comparison:q} \
            --events {input.events:q} \
            --primary {input.primary:q} \
            --preliminary {input.preliminary:q} \
            --masked {input.masked:q} \
            --is-a {input.is_a:q} \
            --is-b {input.is_b:q} \
            --config {params.config:q} \
            --out-primary {output.primary:q} \
            --out-preliminary {output.preliminary:q} \
            --out-strength {output.strength:q} \
            --out-strength-summary {output.strength_summary:q} \
            --out-qc {output.qc:q} \
            --out-summary {output.summary:q} \
            --log {log:q} \
            >> {log:q} 2>&1
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
        config="config/config.yaml",
        benchmark_dir=OUTDIR + "/benchmark/is_boundary",
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.permutation:q})" "$(dirname {log:q})" {params.benchmark_dir:q}
        python {params.script:q} \
            --loci {input.loci:q} \
            --candidate-regions {input.regions:q} \
            --events {input.events:q} \
            --masked {input.masked:q} \
            --config {params.config:q} \
            --output {output.permutation:q} \
            --log {log:q} \
            >> {log:q} 2>&1
        """
