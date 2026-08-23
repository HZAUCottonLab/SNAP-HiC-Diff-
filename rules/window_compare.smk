# Window Compare Module
# Computes window pair metrics: covariates, Gaussian metrics, TW test, KNN matching
#
# Helper functions (get_comparison_pairs, get_pair_samples, get_window_compare_pairs)
# are defined in Snakefile and accessible here via include.
#
# Config keys accessed:
#   - window_compare.pairs_file
#   - window_compare.covariates.cpm_column
#   - window_compare.covariates.nz_column
#   - window_compare.covariates.pseudocount_cpm
#   - window_compare.covariates.weights.EC
#   - window_compare.covariates.weights.IB
#   - window_compare.covariates.weights.SP
#   - window_compare.gaussian_metrics.n_jobs
#   - window_compare.tw.n_jobs
#   - tracy_widom.min_nbin
#   - tracy_widom.na_fraction

# ==============================================================================
# Compute_window_pair_covariates
# ==============================================================================

def get_pair_gap_qc_files(wildcards, sample_index):
    """Resolve one comparison member's QC files through its producing checkpoint."""
    sample = get_pair_samples(wildcards.comparison)[sample_index]
    return checkpoints.compute_gap_qc.get(
        sample=sample,
        window_bins=wildcards.window_bins,
    ).output.qc

rule compute_window_pair_covariates:
    """Compute EC, IB, SP, weighted features for all window pairs of a comparison."""
    input:
        script = SCRIPT_DIR + "/scripts/compute_window_pair_covariates.py",
        qc_a = lambda w: get_pair_gap_qc_files(w, 0),
        qc_b = lambda w: get_pair_gap_qc_files(w, 1),
    output:
        covariates = SCALE_ROOT + "/window_compare/{comparison}/{comparison}_covariates.tsv",
    log:
        SCALE_ROOT + "/logs/window_compare/{comparison}/compute_covariates.log",
    benchmark:
        SCALE_ROOT + "/benchmark/window_compare/{comparison}/compute_covariates.tsv",
    threads: config["parallel"]["threads_per_job"]
    params:
        sample_a = lambda w: get_pair_samples(w.comparison)[0],
        sample_b = lambda w: get_pair_samples(w.comparison)[1],
        chroms = CHROMS,
        cpm_column = config["window_compare"]["covariates"]["cpm_column"],
        nz_column = config["window_compare"]["covariates"]["nz_column"],
        pseudocount = config["window_compare"]["covariates"]["pseudocount_cpm"],
        weight_ec = config["window_compare"]["covariates"]["weights"]["EC"],
        weight_ib = config["window_compare"]["covariates"]["weights"]["IB"],
        weight_sp = config["window_compare"]["covariates"]["weights"]["SP"],
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        COMPARISON="{wildcards.comparison}"
        SAMPLE_A="{params.sample_a}"
        SAMPLE_B="{params.sample_b}"

        mkdir -p "$(dirname "{output.covariates}")" "$(dirname "{log}")"

        QC_A=""
        for f in {input.qc_a}; do
            QC_A="$QC_A $f"
        done

        QC_B=""
        for f in {input.qc_b}; do
            QC_B="$QC_B $f"
        done

        python {input.script} \
            --comparison "$COMPARISON" \
            --sample-a "$SAMPLE_A" \
            --sample-b "$SAMPLE_B" \
            --chroms {params.chroms} \
            --qc-a-files $QC_A \
            --qc-b-files $QC_B \
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


# ==============================================================================
# Compare_gaussian_smoothed_windows
# ==============================================================================

rule compare_gaussian_smoothed_windows:
    """Compute MSE, SSIM, Spearman for Gaussian-smoothed window pairs."""
    input:
        script = SCRIPT_DIR + "/scripts/compare_gaussian_smoothed_windows.py",
        covariates = rules.compute_window_pair_covariates.output.covariates,
        checkpoint_a = lambda w: checkpoints.sliding_window.get(
            sample=get_pair_samples(w.comparison)[0],
            window_bins=w.window_bins,
        ).output.manifest,
        checkpoint_b = lambda w: checkpoints.sliding_window.get(
            sample=get_pair_samples(w.comparison)[1],
            window_bins=w.window_bins,
        ).output.manifest,
    output:
        gaussian = SCALE_ROOT + "/window_compare/{comparison}/{comparison}_gaussian_metrics.tsv",
    log:
        SCALE_ROOT + "/logs/window_compare/{comparison}/compare_gaussian.log",
    benchmark:
        SCALE_ROOT + "/benchmark/window_compare/{comparison}/compare_gaussian.tsv",
    threads: config["window_compare"]["gaussian_metrics"]["n_jobs"]
    params:
        sample_a = lambda w: get_pair_samples(w.comparison)[0],
        sample_b = lambda w: get_pair_samples(w.comparison)[1],
        chroms = CHROMS,
        scale_root = lambda w: SCALE_ROOT.format(window_bins=w.window_bins),
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        COMPARISON="{wildcards.comparison}"
        SAMPLE_A="{params.sample_a}"
        SAMPLE_B="{params.sample_b}"

        mkdir -p "$(dirname "{output.gaussian}")" "$(dirname "{log}")"

        # H5 files are generated by the sliding_window checkpoint
        # Use sample names to construct paths
        H5_A=""
        for chrom in {params.chroms}; do
            H5_A="$H5_A {params.scale_root}/windows/${{SAMPLE_A}}/${{SAMPLE_A}}__${{chrom}}_gaussian.h5"
        done

        H5_B=""
        for chrom in {params.chroms}; do
            H5_B="$H5_B {params.scale_root}/windows/${{SAMPLE_B}}/${{SAMPLE_B}}__${{chrom}}_gaussian.h5"
        done

        python {input.script} \
            --comparison "$COMPARISON" \
            --sample-a "$SAMPLE_A" \
            --sample-b "$SAMPLE_B" \
            --chroms {params.chroms} \
            --covariates "{input.covariates}" \
            --h5-a-files $H5_A \
            --h5-b-files $H5_B \
            --out "{output.gaussian}" \
            --n-jobs {threads} \
            --log "{log}" \
            > "{log}" 2>&1
        """


# ==============================================================================
# Compare_iced_windows_tw
# ==============================================================================

rule compare_iced_windows_tw:
    """Compute Tracy-Widom test for ICE window pairs."""
    input:
        script = SCRIPT_DIR + "/scripts/compare_iced_windows_tw.py",
        covariates = rules.compute_window_pair_covariates.output.covariates,
        checkpoint_ice_a = lambda w: checkpoints.extract_and_balance_sample.get(
            sample=get_pair_samples(w.comparison)[0]
        ).output.manifest,
        checkpoint_ice_b = lambda w: checkpoints.extract_and_balance_sample.get(
            sample=get_pair_samples(w.comparison)[1]
        ).output.manifest,
    output:
        tw = SCALE_ROOT + "/window_compare/{comparison}/{comparison}_tw.tsv",
    log:
        SCALE_ROOT + "/logs/window_compare/{comparison}/compare_tw.log",
    benchmark:
        SCALE_ROOT + "/benchmark/window_compare/{comparison}/compare_tw.tsv",
    threads: config["window_compare"]["tw"]["n_jobs"]
    params:
        sample_a = lambda w: get_pair_samples(w.comparison)[0],
        sample_b = lambda w: get_pair_samples(w.comparison)[1],
        chroms = CHROMS,
        resolution = RESOLUTION,
        min_nbin = config["tracy_widom"]["min_nbin"],
        filter_threshold = config["tracy_widom"]["na_fraction"],
        outdir = OUTDIR,
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        COMPARISON="{wildcards.comparison}"
        SAMPLE_A="{params.sample_a}"
        SAMPLE_B="{params.sample_b}"

        mkdir -p "$(dirname "{output.tw}")" "$(dirname "{log}")"

        python {input.script} \
            --comparison "$COMPARISON" \
            --sample-a "$SAMPLE_A" \
            --sample-b "$SAMPLE_B" \
            --chroms {params.chroms} \
            --covariates "{input.covariates}" \
            --out "{output.tw}" \
            --outdir "{params.outdir}" \
            --resolution {params.resolution} \
            --min-nbin {params.min_nbin} \
            --filter-threshold {params.filter_threshold} \
            --n-jobs {threads} \
            --log "{log}" \
            > "{log}" 2>&1
        """


# ==============================================================================
# Merge_window_pair_metrics
# ==============================================================================

rule merge_window_pair_metrics:
    """Merge covariates, Gaussian metrics, and TW results into a single file."""
    input:
        script = SCRIPT_DIR + "/scripts/merge_window_pair_metrics.py",
        covariates = rules.compute_window_pair_covariates.output.covariates,
        gaussian = rules.compare_gaussian_smoothed_windows.output.gaussian,
        tw = rules.compare_iced_windows_tw.output.tw,
    output:
        merged = SCALE_ROOT + "/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv",
    log:
        SCALE_ROOT + "/logs/window_compare/{comparison}/merge_metrics.log",
    benchmark:
        SCALE_ROOT + "/benchmark/window_compare/{comparison}/merge_metrics.tsv",
    threads: 1
    params:
        sample_a = lambda w: get_pair_samples(w.comparison)[0],
        sample_b = lambda w: get_pair_samples(w.comparison)[1],
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        COMPARISON="{wildcards.comparison}"
        SAMPLE_A="{params.sample_a}"
        SAMPLE_B="{params.sample_b}"

        mkdir -p "$(dirname "{output.merged}")" "$(dirname "{log}")"

        python {input.script} \
            --comparison "$COMPARISON" \
            --sample-a "$SAMPLE_A" \
            --sample-b "$SAMPLE_B" \
            --covariates "{input.covariates}" \
            --gaussian "{input.gaussian}" \
            --tw "{input.tw}" \
            --out "{output.merged}" \
            --log "{log}" \
            > "{log}" 2>&1
        """


# ==============================================================================
# Merge_all_window_pair_metrics
# ==============================================================================

rule merge_all_window_pair_metrics:
    """Merge all comparison-level window pair metrics into one background file."""
    input:
        merged_files = lambda w: expand(
            SCALE_ROOT + "/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv",
            comparison=get_comparison_pairs(), window_bins=w.window_bins,
        )
    output:
        merged_all = SCALE_ROOT + "/window_compare/merged_window_pair_metrics.tsv"
    log:
        SCALE_ROOT + "/logs/window_compare/merge_all_window_pair_metrics.log"
    benchmark:
        SCALE_ROOT + "/benchmark/window_compare/merge_all_window_pair_metrics.tsv"
    threads: 1
    conda:
        "../env/SNAP-HiC-Diff.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.merged_all}")" "$(dirname "{log}")"

        awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input.merged_files} \
            > "{output.merged_all}"

        echo "[INFO] Merged all window_pair_metrics into {output.merged_all}" > "{log}"
        """


