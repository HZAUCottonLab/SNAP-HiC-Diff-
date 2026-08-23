# SNAP-HiC-Diff Window Comparison Pipeline
# Main Snakefile entry point

import os
import glob
import pandas as pd
import yaml

# Load config
configfile: "config/config.yaml"
configfile: "config/chroms.yaml"

# Snakefile directory (resolved at parse time)
SCRIPT_DIR = os.path.dirname(os.path.abspath(workflow.snakefile))
import sys

sys.path.insert(0, os.path.join(SCRIPT_DIR, "scripts"))
from cworld_is_config import parse_cworld_is_settings
from window_analysis_config import parse_window_analysis_settings
CWORLD_IS = parse_cworld_is_settings(config)
# Scripts that read configuration must receive the CLI-selected configuration
# when one is supplied. Snakemake initializes ``configfiles`` with CLI paths
# before it evaluates ``configfile:`` statements in this Snakefile, so the
# default configuration files would otherwise become the last entries. Without
# an explicit CLI file, the first declared configuration is the pipeline's
# complete configuration (the following chroms file is only an overlay).
_CLI_CONFIGFILES = tuple(workflow.overwrite_configfiles)
ACTIVE_CONFIGFILE = str(
    _CLI_CONFIGFILES[-1] if _CLI_CONFIGFILES else workflow.configfiles[0]
)



# Derived paths
OUTDIR = config["outdir"]
# Mapping completeness is a global preflight check. File existence is left to
# the scale-specific validation rule so unrelated shared-IS targets remain runnable.
WINDOW_ANALYSIS = parse_window_analysis_settings(config, check_files=False)
RESOLUTION = WINDOW_ANALYSIS.resolution_bp
WINDOW_BINS = [scale.window_bins for scale in WINDOW_ANALYSIS.scales]
SCALE_ROOT = "{}/window_scales/r{}/w{{window_bins}}".format(OUTDIR, RESOLUTION)
MULTISCALE_ROOT = "{}/multiscale/r{}".format(OUTDIR, RESOLUTION)
TMP = config["tmpdir"]

# ==============================================================================
# Helper Functions
# ==============================================================================

def get_chroms():
    """Get chromosome list."""
    return CHROMS


def get_sample_mcool(sample):
    """Get mcool path for a sample.
    
    Reads data/manifest_mcool.tsv if it exists (output of convert_to_mcool),
    otherwise falls back to the original config manifest.
    Raises KeyError if sample not found in the active manifest.
    """
    if os.path.exists("data/manifest_mcool.tsv"):
        manifest = pd.read_csv("data/manifest_mcool.tsv", sep="\t")
    else:
        manifest = pd.read_csv(config["manifest"], sep="\t")
    row = manifest[manifest["sample_id"] == sample].iloc[0]
    return row["cool_uri"]



def get_sample_pairs():
    """Get all sample pairs to compare."""
    pairs = []
    for a in SAMPLES:
        for b in SAMPLES:
            if a < b:
                pairs.append((a, b))
    return pairs

def get_window_scale(window_bins):
    return WINDOW_ANALYSIS.scale(int(window_bins))


def get_background_preset(window_bins):
    return WINDOW_ANALYSIS.backgrounds_by_window_bins[int(window_bins)]



# Chromosome list (loaded from config/chroms.yaml)
CHROMS = config.get("chroms", [])

# Samples loaded dynamically from manifest.tsv (avoid redundant config)
MANIFEST_PATH = config.get("manifest", "data/manifest.tsv")
if os.path.exists(MANIFEST_PATH):
    SAMPLES = pd.read_csv(MANIFEST_PATH, sep="\t")["sample_id"].tolist()
else:
    SAMPLES = config.get("samples", [])


# ==============================================================================
# Window Compare Helper Functions
# ==============================================================================

def get_comparison_pairs():
    """Load comparison pairs from config."""
    pairs_file = config.get("window_compare", {}).get("pairs_file", "data/window_compare_pairs.yaml")
    if os.path.exists(pairs_file):
        with open(pairs_file, "r") as f:
            data = yaml.safe_load(f)
            return data.get("comparisons", [])
    return []


def get_window_compare_pairs():
    """Get list of comparison names for window_compare."""
    return [p["name"] for p in get_comparison_pairs()]


def get_pair_samples(comparison):
    """Get sample_a and sample_b for a given comparison name."""
    pairs = get_comparison_pairs()
    for p in pairs:
        if p["name"] == comparison:
            return p["sample_a"], p["sample_b"]
    return None, None


def get_comparison_is_file(comparison, member_index):
    sample = get_pair_samples(comparison)[member_index]
    return OUTDIR + "/IS/{0}/{0}_IS.tsv".format(sample)


def knn_heatmap_chrom_done_by_comp(wildcards):
    """Return per-chrom done markers for a comparison.

    Uses checkpoint output so dry-run schedules checkpoint job.
    """
    ck = checkpoints.plot_knn_significant_heatmap.get(
        comparison=wildcards.comparison,
        window_bins=wildcards.window_bins,
    )
    manifest = ck.output.manifest
    df = pd.read_csv(manifest, sep="\t")
    chroms = sorted(df["chrom"].unique())
    scale_root = SCALE_ROOT.format(window_bins=wildcards.window_bins)
    return [
        f"{scale_root}/logs/knn_significant_heatmap/{wildcards.comparison}/{chrom}.done"
        for chrom in chroms
    ]


# ==============================================================================
# Include Rule Files
# ==============================================================================

include: "rules/convert_to_mcool.smk"
include: "rules/validate_inputs.smk"
include: "rules/build_windows.smk"
include: "rules/extract_chrom.smk"
include: "rules/extract_expected.smk"
include: "rules/compute_oe_matrix.smk"
include: "rules/compute_gap_qc.smk"
include: "rules/sliding_window.smk"
include: "rules/window_compare.smk"
include: "rules/knn_empirical_pvalue.smk"
include: "rules/knn_significant_heatmap.smk"
include: "rules/compute_IS.smk"
include: "rules/meltron.smk"
include: "rules/is_change.smk"

# ==============================================================================
# Default Rules
# ==============================================================================

rule all:
    input:
        # 0. Convert .cool to .mcool (if needed) → data/manifest_mcool.tsv
        "data/manifest_mcool.tsv",
        # 1. Validated manifest
        OUTDIR + "/validated/manifest.validated.tsv",
        # 2. Windows
        expand(SCALE_ROOT + "/windows/windows.parquet", window_bins=WINDOW_BINS),
        # 3. Extract + Balance (one job per sample, for-loop over chroms in shell)
        expand(OUTDIR + "/logs/extract_balance/{sample}.done", sample=SAMPLES),
        # 4. Expected files (checkpoint manifests)
        expand(OUTDIR + "/logs/expected/{sample}_cooltools.done", sample=SAMPLES),
        expand(OUTDIR + "/logs/expected/{sample}_to_npz.done", sample=SAMPLES),
        # 5. O/E matrix (checkpoint manifest)
        expand(OUTDIR + "/logs/oe/{sample}.done", sample=SAMPLES),
        # 6. Gap QC + Sliding window (checkpoint manifests)
        expand(
            SCALE_ROOT + "/logs/gap_qc/{sample}.done",
            sample=SAMPLES,
            window_bins=WINDOW_BINS,
        ),
        expand(
            SCALE_ROOT + "/logs/windows/{sample}.done",
            sample=SAMPLES,
            window_bins=WINDOW_BINS,
        ),
        # 7. Window compare — every intermediate and final metric is scale-qualified.
        expand(
            SCALE_ROOT + "/window_compare/{comparison}/{comparison}_covariates.tsv",
            comparison=get_window_compare_pairs(),
            window_bins=WINDOW_BINS,
        ),
        expand(
            SCALE_ROOT + "/window_compare/{comparison}/{comparison}_gaussian_metrics.tsv",
            comparison=get_window_compare_pairs(),
            window_bins=WINDOW_BINS,
        ),
        expand(
            SCALE_ROOT + "/window_compare/{comparison}/{comparison}_tw.tsv",
            comparison=get_window_compare_pairs(),
            window_bins=WINDOW_BINS,
        ),
        expand(
            SCALE_ROOT + "/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv",
            comparison=get_window_compare_pairs(),
            window_bins=WINDOW_BINS,
        ),
        # 8. KNN calibration and its diagnostic heatmaps run once per scale.
        expand(
            SCALE_ROOT + "/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
            comparison=get_window_compare_pairs(),
            window_bins=WINDOW_BINS,
        ),
        expand(
            SCALE_ROOT + "/logs/knn_significant_heatmap/{comparison}/{chrom}.done",
            comparison=get_window_compare_pairs(),
            chrom=get_chroms(),
            window_bins=WINDOW_BINS,
        ),
        # 9. CWorld IS is shared by all five genomic-window branches.
        expand(OUTDIR + "/IS/{sample}/{sample}_IS.tsv", sample=SAMPLES),
        # 10. Meltron is evaluated independently for every genomic window scale.
        expand(
            SCALE_ROOT + "/meltron/{comparison}/{comparison}_meltron_merged.tsv",
            comparison=get_window_compare_pairs(),
            window_bins=WINDOW_BINS,
        ),
        # 11. Publish only aggregate multiscale candidate and IS-change products.
        expand(
            MULTISCALE_ROOT + "/{comparison}/candidate_union.tsv",
            comparison=get_window_compare_pairs(),
        ),
        expand(
            MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci.tsv",
            comparison=get_window_compare_pairs(),
        ),
        expand(
            MULTISCALE_ROOT + "/{comparison}/is_change/figures/manifest.tsv",
            comparison=get_window_compare_pairs(),
        ),


rule clean:
    """Remove all output files."""
    shell:
        "rm -rf {OUTDIR}/* logs/*"
