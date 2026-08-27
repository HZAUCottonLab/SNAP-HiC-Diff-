# Window-pair comparison for each independent scale branch.

rule compute_window_pair_covariates:
    """Compute EC, IB, and SP covariates for one comparison and scale."""
    input:
        script=SCRIPT_DIR + "/scripts/compute_window_pair_covariates.py",
        qc_a=lambda w: [
            "{}/windows/{}/{}__{}_window_qc.tsv".format(
                get_scale_outdir(w), get_pair_samples(w.comparison)[0],
                get_pair_samples(w.comparison)[0], chrom,
            )
            for chrom in CHROMS
        ],
        qc_b=lambda w: [
            "{}/windows/{}/{}__{}_window_qc.tsv".format(
                get_scale_outdir(w), get_pair_samples(w.comparison)[1],
                get_pair_samples(w.comparison)[1], chrom,
            )
            for chrom in CHROMS
        ],
    output:
        covariates=SCALE_ROOT + "/{scale}/window_compare/{comparison}/{comparison}_covariates.tsv",
    log:
        SCALE_ROOT + "/{scale}/logs/window_compare/{comparison}/compute_covariates.log",
    benchmark:
        SCALE_ROOT + "/{scale}/benchmark/window_compare/{comparison}/compute_covariates.tsv",
    threads:
        config["parallel"]["threads_per_job"]
    params:
        sample_a=lambda w: get_pair_samples(w.comparison)[0],
        sample_b=lambda w: get_pair_samples(w.comparison)[1],
        chroms=CHROMS,
        cpm_column=config["window_compare"]["covariates"]["cpm_column"],
        nz_column=config["window_compare"]["covariates"]["nz_column"],
        pseudocount=config["window_compare"]["covariates"]["pseudocount_cpm"],
        weight_ec=config["window_compare"]["covariates"]["weights"]["EC"],
        weight_ib=config["window_compare"]["covariates"]["weights"]["IB"],
        weight_sp=config["window_compare"]["covariates"]["weights"]["SP"],
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.covariates}")" "$(dirname "{log}")"
        python {input.script} \
            --comparison "{wildcards.comparison}" \
            --sample-a "{params.sample_a}" \
            --sample-b "{params.sample_b}" \
            --chroms {params.chroms} \
            --qc-a-files {input.qc_a} \
            --qc-b-files {input.qc_b} \
            --out "{output.covariates}" \
            --cpm-column "{params.cpm_column}" \
            --nz-column "{params.nz_column}" \
            --pseudocount {params.pseudocount} \
            --weight-ec {params.weight_ec} \
            --weight-ib {params.weight_ib} \
            --weight-sp {params.weight_sp} \
            --n-jobs {threads} \
            --log "{log}" \
            > "{log}" 2>&1
        """


rule compare_gaussian_smoothed_windows:
    """Compute Gaussian window metrics for one comparison and scale."""
    input:
        script=SCRIPT_DIR + "/scripts/compare_gaussian_smoothed_windows.py",
        covariates=rules.compute_window_pair_covariates.output.covariates,
        checkpoint_a=lambda w: checkpoints.sliding_window.get(
            scale=w.scale, sample=get_pair_samples(w.comparison)[0]
        ).output.manifest,
        checkpoint_b=lambda w: checkpoints.sliding_window.get(
            scale=w.scale, sample=get_pair_samples(w.comparison)[1]
        ).output.manifest,
    output:
        gaussian=SCALE_ROOT + "/{scale}/window_compare/{comparison}/{comparison}_gaussian_metrics.tsv",
    log:
        SCALE_ROOT + "/{scale}/logs/window_compare/{comparison}/compare_gaussian.log",
    benchmark:
        SCALE_ROOT + "/{scale}/benchmark/window_compare/{comparison}/compare_gaussian.tsv",
    threads:
        config["window_compare"]["gaussian_metrics"]["n_jobs"]
    params:
        sample_a=lambda w: get_pair_samples(w.comparison)[0],
        sample_b=lambda w: get_pair_samples(w.comparison)[1],
        chroms=CHROMS,
        scale_outdir=get_scale_outdir,
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.gaussian}")" "$(dirname "{log}")"
        H5_A=""
        H5_B=""
        for chrom in {params.chroms}; do
            H5_A="$H5_A {params.scale_outdir}/windows/{params.sample_a}/{params.sample_a}__${{chrom}}_gaussian.h5"
            H5_B="$H5_B {params.scale_outdir}/windows/{params.sample_b}/{params.sample_b}__${{chrom}}_gaussian.h5"
        done
        python {input.script} \
            --comparison "{wildcards.comparison}" \
            --sample-a "{params.sample_a}" \
            --sample-b "{params.sample_b}" \
            --chroms {params.chroms} \
            --covariates "{input.covariates}" \
            --h5-a-files $H5_A \
            --h5-b-files $H5_B \
            --out "{output.gaussian}" \
            --n-jobs {threads} \
            --log "{log}" \
            > "{log}" 2>&1
        """


rule compare_iced_windows_tw:
    """Compute Tracy-Widom statistics with shared balanced cool inputs."""
    input:
        script=SCRIPT_DIR + "/scripts/compare_iced_windows_tw.py",
        covariates=rules.compute_window_pair_covariates.output.covariates,
        checkpoint_ice_a=lambda w: checkpoints.extract_and_balance_sample.get(
            sample=get_pair_samples(w.comparison)[0]
        ).output.manifest,
        checkpoint_ice_b=lambda w: checkpoints.extract_and_balance_sample.get(
            sample=get_pair_samples(w.comparison)[1]
        ).output.manifest,
    output:
        tw=SCALE_ROOT + "/{scale}/window_compare/{comparison}/{comparison}_tw.tsv",
    log:
        SCALE_ROOT + "/{scale}/logs/window_compare/{comparison}/compare_tw.log",
    benchmark:
        SCALE_ROOT + "/{scale}/benchmark/window_compare/{comparison}/compare_tw.tsv",
    threads:
        config["window_compare"]["tw"]["n_jobs"]
    params:
        sample_a=lambda w: get_pair_samples(w.comparison)[0],
        sample_b=lambda w: get_pair_samples(w.comparison)[1],
        chroms=CHROMS,
        resolution=lambda w: get_scale_value(w, "resolution_bp"),
        min_nbin=config["tracy_widom"]["min_nbin"],
        filter_threshold=config["tracy_widom"]["na_fraction"],
        common_outdir=OUTDIR,
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tw}")" "$(dirname "{log}")"
        python {input.script} \
            --comparison "{wildcards.comparison}" \
            --sample-a "{params.sample_a}" \
            --sample-b "{params.sample_b}" \
            --chroms {params.chroms} \
            --covariates "{input.covariates}" \
            --out "{output.tw}" \
            --outdir "{params.common_outdir}" \
            --resolution {params.resolution} \
            --min-nbin {params.min_nbin} \
            --filter-threshold {params.filter_threshold} \
            --n-jobs {threads} \
            --log "{log}" \
            > "{log}" 2>&1
        """


rule merge_window_pair_metrics:
    """Merge scale-specific covariates, Gaussian, and TW tables."""
    input:
        script=SCRIPT_DIR + "/scripts/merge_window_pair_metrics.py",
        covariates=rules.compute_window_pair_covariates.output.covariates,
        gaussian=rules.compare_gaussian_smoothed_windows.output.gaussian,
        tw=rules.compare_iced_windows_tw.output.tw,
    output:
        merged=SCALE_ROOT + "/{scale}/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv",
    log:
        SCALE_ROOT + "/{scale}/logs/window_compare/{comparison}/merge_metrics.log",
    benchmark:
        SCALE_ROOT + "/{scale}/benchmark/window_compare/{comparison}/merge_metrics.tsv",
    threads:
        1
    params:
        sample_a=lambda w: get_pair_samples(w.comparison)[0],
        sample_b=lambda w: get_pair_samples(w.comparison)[1],
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.merged}")" "$(dirname "{log}")"
        python {input.script} \
            --comparison "{wildcards.comparison}" \
            --sample-a "{params.sample_a}" \
            --sample-b "{params.sample_b}" \
            --covariates "{input.covariates}" \
            --gaussian "{input.gaussian}" \
            --tw "{input.tw}" \
            --out "{output.merged}" \
            --log "{log}" \
            > "{log}" 2>&1
        """
