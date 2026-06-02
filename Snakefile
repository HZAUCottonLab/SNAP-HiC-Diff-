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

# Derived paths
OUTDIR = config["outdir"]
RESOLUTION = config["resolution"]
WINDOW_SIZE = config["window_size"]
TMP = config["tmpdir"]
STEP_SIZE = config.get("step_size", config["window_size"])  # Default: non-overlapping

# ==============================================================================
# Helper Functions
# ==============================================================================

def get_chroms():
    """Get chromosome list."""
    return CHROMS


def get_sample_mcool(sample):
    """Get mcool path for a sample."""
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


def knn_heatmap_chrom_done_by_comp(wildcards):
    """Return per-chrom done markers for a comparison.

    Uses checkpoint output so dry-run schedules checkpoint job.
    """
    ck = checkpoints.plot_knn_significant_heatmap.get(comparison=wildcards.comparison)
    manifest = ck.output.manifest
    df = pd.read_csv(manifest, sep="\t")
    chroms = sorted(df["chrom"].unique())
    return [
        f"{OUTDIR}/logs/knn_significant_heatmap/{wildcards.comparison}/{chrom}.done"
        for chrom in chroms
    ]


# ==============================================================================
# Include Rule Files
# ==============================================================================

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

# ==============================================================================
# Default Rules
# ==============================================================================

rule all:
    input:
        # 1. Validated manifest
        OUTDIR + "/validated/manifest.validated.tsv",
        # 2. Windows
        OUTDIR + "/windows/windows_" + str(RESOLUTION) + ".parquet",
        # 3. Extract + Balance (one job per sample, for-loop over chroms in shell)
        expand(OUTDIR + "/logs/extract_balance/{sample}.done", sample=SAMPLES),
        # 4. Expected files (checkpoint manifests)
        expand(OUTDIR + "/logs/expected/{sample}_cooltools.done", sample=SAMPLES),
        expand(OUTDIR + "/logs/expected/{sample}_to_npz.done", sample=SAMPLES),
        # 5. O/E matrix (checkpoint manifest)
        expand(OUTDIR + "/logs/oe/{sample}.done", sample=SAMPLES),
        # 6. Gap QC + Sliding window (checkpoint manifests)
        expand(OUTDIR + "/logs/gap_qc/{sample}.done", sample=SAMPLES),
        expand(OUTDIR + "/logs/windows/{sample}.done", sample=SAMPLES),
        # 7. Window compare — intermediate stage outputs (explicit, so deletion triggers rebuild)
        expand(OUTDIR + "/window_compare/{comparison}/{comparison}_covariates.tsv",
               comparison=get_window_compare_pairs()),
        expand(OUTDIR + "/window_compare/{comparison}/{comparison}_gaussian_metrics.tsv",
               comparison=get_window_compare_pairs()),
        expand(OUTDIR + "/window_compare/{comparison}/{comparison}_tw.tsv",
               comparison=get_window_compare_pairs()),
        # 8. Window compare — final merged output
        expand(
            OUTDIR + "/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv",
            comparison=get_window_compare_pairs()
        ),
        # 9. KNN empirical P-values (per comparison pair)
        expand(
            OUTDIR + "/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
            comparison=get_window_compare_pairs()
        ),
        # 10. KNN significant window heatmaps (per-chrom done markers)
        expand(
            OUTDIR + "/logs/knn_significant_heatmap/{comparison}/{chrom}.done",
            comparison=get_window_compare_pairs(),
            chrom=get_chroms(),
        ),


rule clean:
    """Remove all output files."""
    shell:
        "rm -rf {OUTDIR}/* logs/*"
