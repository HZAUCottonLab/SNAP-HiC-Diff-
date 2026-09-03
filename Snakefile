# SNAP-HiC-Diff Window Comparison Pipeline
# Main Snakefile entry point

import json
import os
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
from is_boundary_config import parse_is_boundary_settings
from tad_config import parse_tad_settings
from window_scale_config import build_window_scales

CWORLD_IS = parse_cworld_is_settings(config)
TAD = parse_tad_settings(config)
IS_BOUNDARY = parse_is_boundary_settings(config, CWORLD_IS, TAD)
TAD_RESOLUTIONS = TAD.resolutions
TAD_SIGNALS = TAD.signals
TAD_LABELS = {resolution: "{}kb".format(resolution // 1000) for resolution in TAD_RESOLUTIONS}
TAD_RESOLUTION_BY_LABEL = {label: resolution for resolution, label in TAD_LABELS.items()}
TAD_RESOLUTION_MAP = ",".join(
    "{}={}".format(label, resolution)
    for label, resolution in sorted(TAD_RESOLUTION_BY_LABEL.items())
)
TAD_D_CONFIG_JSON = json.dumps(TAD.d, sort_keys=True)


# Derived paths
OUTDIR = config["outdir"]
RESOLUTION = config["resolution"]
TMP = config["tmpdir"]
WINDOW_SCALES = build_window_scales(config)
SCALE_LABELS = tuple(WINDOW_SCALES)
SCALE_ROOT = "{}/window_scales/r{}".format(config["outdir"], config["resolution"])
MULTISCALE_ROOT = "{}/multiscale/r{}".format(config["outdir"], config["resolution"])


def get_scale_value(wildcards, key):
    return WINDOW_SCALES[wildcards.scale][key]


def get_scale_outdir(wildcards):
    return "{}/{}".format(SCALE_ROOT, wildcards.scale)


wildcard_constraints:
    scale="w25|w50|w75|w100|w125"

def tad_run_dir(sample, resolution, signal):
    return "{}/{}/{}/{}/D".format(TAD.run_root, sample, signal, TAD_LABELS[resolution])


def tad_ready_path(sample, resolution, signal):
    return "{}/tad_prepare/{}.{}.{}.ready".format(OUTDIR, sample, TAD_LABELS[resolution], signal)


def tad_bed_path(sample, resolution, signal):
    return "{}/tads/{}.{}.{}.3dnetmod.bed".format(OUTDIR, sample, TAD_LABELS[resolution], signal)


def tad_bedpe_path(sample, resolution, signal):
    return "{}/tads/{}.{}.{}.3dnetmod.bedpe".format(OUTDIR, sample, TAD_LABELS[resolution], signal)


def tad_boundary_path(sample, resolution, signal):
    return "{}/tad_boundaries/{}.{}.{}.boundaries.bed".format(OUTDIR, sample, TAD_LABELS[resolution], signal)


def tad_merged_boundary_path(sample, signal):
    return "{}/tad_boundaries/{}.multires.{}.boundaries.bed".format(OUTDIR, sample, signal)


def unified_boundary_path(sample):
    return "{}/tad_boundaries/{}.multires.unified.boundaries.tsv".format(OUTDIR, sample)


def boundary_jitter_path(sample):
    return "{}/tad_boundaries/{}.boundary_jitter_offsets.tsv".format(OUTDIR, sample)


def boundary_view_path(sample):
    return "{}/tad_boundaries/{}.annotation_views.tsv".format(OUTDIR, sample)


def boundary_tad_unevaluable_path(sample):
    return "{}/tad_boundaries/{}.tad_unevaluable.tsv".format(OUTDIR, sample)


def boundary_event_path(comparison):
    return "{}/tad_boundary_events/{}.boundary_events.tsv".format(OUTDIR, comparison)


def boundary_stability_path(comparison):
    return "{}/tad_boundary_events/{}.boundary_event_stability.tsv".format(OUTDIR, comparison)


def annotated_locus_path(comparison, tier):
    suffix = "_preliminary" if tier == "preliminary" else ""
    return MULTISCALE_ROOT + "/{}/is_change/{}_is_change_loci{}.boundary_annotated.tsv".format(
        comparison, comparison, suffix
    )

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


# ==============================================================================
# Include Rule Files
# ==============================================================================

include: "rules/convert_to_mcool.smk"
include: "rules/validate_inputs.smk"
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
include: "rules/is_boundary_annotation.smk"

# ==============================================================================
# Default Rules
# ==============================================================================

COMMON_ALL_OUTPUTS = [
    "data/manifest_mcool.tsv",
    OUTDIR + "/validated/manifest.validated.tsv",
] + expand(
    OUTDIR + "/logs/extract_balance/{sample}.done", sample=SAMPLES
) + expand(
    OUTDIR + "/logs/expected/{sample}_{stage}.done",
    sample=SAMPLES, stage=["cooltools", "to_npz"],
) + expand(
    OUTDIR + "/logs/oe/{sample}.done", sample=SAMPLES
) + expand(
    OUTDIR + "/IS/{sample}/{sample}_IS.tsv", sample=SAMPLES
)

SCALE_ALL_OUTPUTS = (
    expand(
        SCALE_ROOT + "/{scale}/audits/knn_background.json",
        scale=SCALE_LABELS,
    )
    + expand(
        SCALE_ROOT + "/{scale}/logs/{stage}/{sample}.done",
        scale=SCALE_LABELS,
        stage=["gap_qc", "windows"],
        sample=SAMPLES,
    )
    + expand(
        SCALE_ROOT + "/{scale}/window_compare/{comparison}/{comparison}_covariates.tsv",
        scale=SCALE_LABELS,
        comparison=get_window_compare_pairs(),
    )
    + expand(
        SCALE_ROOT + "/{scale}/window_compare/{comparison}/{comparison}_gaussian_metrics.tsv",
        scale=SCALE_LABELS,
        comparison=get_window_compare_pairs(),
    )
    + expand(
        SCALE_ROOT + "/{scale}/window_compare/{comparison}/{comparison}_tw.tsv",
        scale=SCALE_LABELS,
        comparison=get_window_compare_pairs(),
    )
    + expand(
        SCALE_ROOT + "/{scale}/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv",
        scale=SCALE_LABELS,
        comparison=get_window_compare_pairs(),
    )
    + expand(
        SCALE_ROOT + "/{scale}/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv",
        scale=SCALE_LABELS,
        comparison=get_window_compare_pairs(),
    )
    + expand(
        SCALE_ROOT + "/{scale}/logs/knn_significant_heatmap/{comparison}/{chrom}.done",
        scale=SCALE_LABELS,
        comparison=get_window_compare_pairs(),
        chrom=get_chroms(),
    )
    + expand(
        SCALE_ROOT + "/{scale}/meltron/{comparison}/{comparison}_meltron_merged.tsv",
        scale=SCALE_LABELS,
        comparison=get_window_compare_pairs(),
    )
)

MULTISCALE_ALL_OUTPUTS = (
    expand(
        MULTISCALE_ROOT + "/{comparison}/all_scale_candidates.tsv",
        comparison=get_window_compare_pairs(),
    )
    + expand(
        MULTISCALE_ROOT + "/{comparison}/candidate_union.tsv",
        comparison=get_window_compare_pairs(),
    )
    + expand(
        MULTISCALE_ROOT + "/{comparison}/candidate_union.bed",
        comparison=get_window_compare_pairs(),
    )
    + expand(
        MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci.tsv",
        comparison=get_window_compare_pairs(),
    )
    + expand(
        MULTISCALE_ROOT + "/{comparison}/is_change/figures/manifest.tsv",
        comparison=get_window_compare_pairs(),
    )
)

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

IS_BOUNDARY_ALL_OUTPUTS = (
    expand(
        OUTDIR + "/tad_boundaries/{sample}.multires.unified.boundaries.tsv",
        sample=SAMPLES,
    )
    + expand(
        OUTDIR + "/tad_boundaries/{sample}.boundary_jitter_offsets.tsv",
        sample=SAMPLES,
    )
    + expand(
        OUTDIR + "/tad_boundaries/{sample}.annotation_views.tsv",
        sample=SAMPLES,
    )
    + expand(
        OUTDIR + "/tad_boundaries/{sample}.tad_unevaluable.tsv",
        sample=SAMPLES,
    )
    + expand(
        OUTDIR + "/tad_boundary_events/{comparison}.boundary_events.tsv",
        comparison=get_window_compare_pairs(),
    )
    + expand(
        OUTDIR + "/tad_boundary_events/{comparison}.boundary_event_stability.tsv",
        comparison=get_window_compare_pairs(),
    )
    + expand(
        OUTDIR + "/tad_boundary_events/{comparison}.matching_qc.tsv",
        comparison=get_window_compare_pairs(),
    )
    + expand(
        OUTDIR + "/tad_boundary_events/{comparison}.boundary_strength_by_scale.tsv",
        comparison=get_window_compare_pairs(),
    )
    + expand(
        OUTDIR + "/tad_boundary_events/{comparison}.boundary_strength_summary.tsv",
        comparison=get_window_compare_pairs(),
    )
    + expand(
        OUTDIR + "/tad_boundary_events/{comparison}.annotation_qc.tsv",
        comparison=get_window_compare_pairs(),
    )
    + expand(
        OUTDIR + "/tad_boundary_events/{comparison}.annotation_summary.tsv",
        comparison=get_window_compare_pairs(),
    )
    + expand(
        MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci.boundary_annotated.tsv",
        comparison=get_window_compare_pairs(),
    )
    + expand(
        MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci_preliminary.boundary_annotated.tsv",
        comparison=get_window_compare_pairs(),
    )
    + expand(
        MULTISCALE_ROOT + "/{comparison}/is_change/figures/manifest.tsv",
        comparison=get_window_compare_pairs(),
    )
    if IS_BOUNDARY.enabled
    else []
)

PERMUTATION_ALL_OUTPUTS = (
    expand(
        OUTDIR + "/tad_boundary_events/{comparison}.{tier}.permutation.tsv",
        comparison=get_window_compare_pairs(), tier=("primary", "preliminary"),
    )
    if IS_BOUNDARY.permutation.enabled
    else []
)

rule all:
    input:
        COMMON_ALL_OUTPUTS,
        SCALE_ALL_OUTPUTS,
        MULTISCALE_ALL_OUTPUTS,
        TAD_ALL_OUTPUTS,
        IS_BOUNDARY_ALL_OUTPUTS,
        PERMUTATION_ALL_OUTPUTS,
    default_target: True


rule clean:
    """Remove all output files."""
    shell:
        "rm -rf {OUTDIR}/* logs/*"
