# SNAP-HiC-Diff Window Comparison Pipeline
# Main Snakefile entry point

import json
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
from tad_config import parse_tad_settings

CWORLD_IS = parse_cworld_is_settings(config)
TAD = parse_tad_settings(config)
TAD_RESOLUTIONS = TAD.resolutions
TAD_SIGNALS = TAD.signals
TAD_LABELS = {resolution: f"{resolution // 1000}kb" for resolution in TAD_RESOLUTIONS}
TAD_RESOLUTION_BY_LABEL = {label: resolution for resolution, label in TAD_LABELS.items()}
TAD_D_CONFIG_JSON = json.dumps(TAD.d, sort_keys=True)


# Derived paths
OUTDIR = config["outdir"]
RESOLUTION = config["resolution"]
WINDOW_SIZE = config["window_size"]
TMP = config["tmpdir"]
STEP_SIZE = config.get("step_size", config["window_size"])  # Default: non-overlapping

def tad_run_dir(sample, resolution, signal):
    return f"{TAD.run_root}/{sample}/{signal}/{TAD_LABELS[resolution]}/D"


def tad_ready_path(sample, resolution, signal):
    return f"{OUTDIR}/tad_prepare/{sample}.{TAD_LABELS[resolution]}.{signal}.ready"


def tad_bed_path(sample, resolution, signal):
    return f"{OUTDIR}/tads/{sample}.{TAD_LABELS[resolution]}.{signal}.3dnetmod.bed"


def tad_bedpe_path(sample, resolution, signal):
    return f"{OUTDIR}/tads/{sample}.{TAD_LABELS[resolution]}.{signal}.3dnetmod.bedpe"


def tad_boundary_path(sample, resolution, signal):
    return f"{OUTDIR}/tad_boundaries/{sample}.{TAD_LABELS[resolution]}.{signal}.boundaries.bed"


def tad_merged_boundary_path(sample, signal):
    return f"{OUTDIR}/tad_boundaries/{sample}.multires.{signal}.boundaries.bed"

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
include: "rules/tad.smk"

# ==============================================================================
# Default Rules
# ==============================================================================

TAD_ALL_OUTPUTS = [
    output
    for sample in SAMPLES
    for resolution in TAD_RESOLUTIONS
    for signal in TAD_SIGNALS
    for output in (
        tad_bed_path(sample, resolution, signal),
        tad_bedpe_path(sample, resolution, signal),
        tad_boundary_path(sample, resolution, signal),
    )
] + [
    tad_merged_boundary_path(sample, signal)
    for sample in SAMPLES
    for signal in TAD_SIGNALS
]

rule all:
    input:
        # 0. Convert .cool to .mcool (if needed) → data/manifest_mcool.tsv
        "data/manifest_mcool.tsv",
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
        # 11. Compute IS per sample (checkpoint → per-chrom jobs → merge)
        expand(OUTDIR + "/IS/{sample}/{sample}_IS.tsv", sample=SAMPLES),
        # 12. Meltron per comparison (both directions, merged)
        expand(
            OUTDIR + "/meltron/{comparison}/{comparison}_meltron_merged.tsv",
            comparison=get_window_compare_pairs()
        ),
        # 13. IS-change loci and diagnostic plot manifests
        expand(
            OUTDIR + "/is_change/{comparison}/{comparison}_is_change_loci.tsv",
            comparison=get_window_compare_pairs()
        ),
        expand(
            OUTDIR + "/figures/is_change/{comparison}/manifest.tsv",
            comparison=get_window_compare_pairs()
        ),
        TAD_ALL_OUTPUTS,


rule clean:
    """Remove all output files."""
    shell:
        "rm -rf {OUTDIR}/* logs/*"
