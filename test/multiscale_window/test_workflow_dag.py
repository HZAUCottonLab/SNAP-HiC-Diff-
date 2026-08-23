from math import ceil
from pathlib import Path
import re
import subprocess

import yaml


REPO = Path(__file__).resolve().parents[2]
COMPARISON = "H1398_vs_L134"
WINDOW_BINS = [25, 50, 75, 100, 125]
BACKGROUND_COLUMNS = [
    "comparison", "chrom", "window_start", "window_end",
    "EC", "IB", "SP", "MSE", "SSIM", "Spearman", "tw_pvalue", "NCC",
]


def write_background(path, resolution_bp, window_bins):
    """Write the smallest schema-valid external background preset fixture."""
    window_size_bp = resolution_bp * window_bins
    row = {
        "comparison": "background_a_vs_background_b",
        "chrom": "chr1",
        "window_start": 0,
        "window_end": window_size_bp,
        "EC": 1.0,
        "IB": 1.0,
        "SP": 1.0,
        "MSE": 1.0,
        "SSIM": 0.9,
        "Spearman": 0.8,
        "tw_pvalue": 0.1,
        "NCC": 0.7,
    }
    path.write_text(
        "\t".join(BACKGROUND_COLUMNS)
        + "\n"
        + "\t".join(str(row[column]) for column in BACKGROUND_COLUMNS)
        + "\n"
    )


def make_config(tmp_path, resolution_bp, missing_window_bins=()):
    """Build a portable full config with one external preset per active scale."""
    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    config["outdir"] = str(tmp_path / "results")
    config["tmpdir"] = str(tmp_path / "tmp")
    config["resolution"] = resolution_bp
    config["convert"]["resolutions"] = [resolution_bp]
    config["window_analysis"] = {
        "window_bins": WINDOW_BINS,
        "step_fraction": 0.20,
    }
    config["knn_empirical_pvalue"]["knn_k"] = 1

    presets = []
    for window_bins in WINDOW_BINS:
        if window_bins in missing_window_bins:
            continue
        background = tmp_path / f"background_r{resolution_bp}_w{window_bins}.tsv"
        write_background(background, resolution_bp, window_bins)
        presets.append(
            {
                "resolution_bp": resolution_bp,
                "window_bins": window_bins,
                "step_bins": ceil(window_bins * 0.20),
                "file": str(background),
            }
        )
    config["knn_empirical_pvalue"]["background_presets"] = presets

    # The dry-run must also remain valid when its source resolution changes.
    is_distances_bp = [resolution_bp * factor for factor in (8, 15, 23, 30, 38, 45, 53, 60, 68, 75)]
    cworld_is = config["meltron"]["cworld_is"]
    cworld_is["source_resolution_bp"] = resolution_bp
    cworld_is["is_distances_bp"] = is_distances_bp
    is_change = config["is_change"]
    is_change["context_bp"] = resolution_bp * 80
    is_change["spatial_cluster_tolerance_bp"] = resolution_bp * 3
    is_change["scale_bands"] = {
        "small": is_distances_bp[:3],
        "medium": is_distances_bp[3:6],
        "large": is_distances_bp[6:],
    }

    config_path = tmp_path / f"config-r{resolution_bp}.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    return config_path, Path(config["outdir"])


def dry_run(config_path, target):
    return subprocess.run(
        [
            "snakemake",
            "-n",
            "-p",
            "--configfile",
            str(config_path),
            "--",
            str(target),
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )


def job_count(plan, rule_name):
    match = re.search(rf"^\s*{rule_name}\s+(\d+)\s*$", plan, re.MULTILINE)
    assert match, plan
    return int(match.group(1))


def test_final_candidate_union_dag_expands_five_scales_at_4kb_and_shares_is(tmp_path):
    config_path, outdir = make_config(tmp_path, 4000)
    target = outdir / "multiscale" / "r4000" / COMPARISON / "candidate_union.tsv"

    result = dry_run(config_path, target)
    plan = result.stdout + result.stderr

    assert result.returncode == 0, plan
    for window_bins in WINDOW_BINS:
        assert f"window_scales/r4000/w{window_bins}" in plan
    assert job_count(plan, "meltron_per_comparison") == len(WINDOW_BINS)
    assert job_count(plan, "compute_IS_chrom") == 52
    assert job_count(plan, "merge_IS") == 2


def test_final_candidate_union_dag_keeps_bin_tags_at_10kb_and_shares_is(tmp_path):
    config_path, outdir = make_config(tmp_path, 10000)
    target = outdir / "multiscale" / "r10000" / COMPARISON / "candidate_union.tsv"

    result = dry_run(config_path, target)
    plan = result.stdout + result.stderr

    assert result.returncode == 0, plan
    for window_bins in WINDOW_BINS:
        assert f"window_scales/r10000/w{window_bins}" in plan
    assert job_count(plan, "meltron_per_comparison") == len(WINDOW_BINS)
    assert job_count(plan, "compute_IS_chrom") == 52
    assert job_count(plan, "merge_IS") == 2


def test_missing_scale_preset_fails_with_the_exact_expected_tuple(tmp_path):
    config_path, outdir = make_config(tmp_path, 4000, missing_window_bins=(125,))
    target = outdir / "multiscale" / "r4000" / COMPARISON / "candidate_union.tsv"

    result = dry_run(config_path, target)
    plan = result.stdout + result.stderr

    assert result.returncode != 0
    assert "missing background preset for r4000/w125 step_bins=25" in plan


def test_rule_all_only_publishes_scale_qualified_and_multiscale_final_targets():
    snakefile = (REPO / "Snakefile").read_text()
    rule_all = snakefile[snakefile.index("rule all:"):snakefile.index("\n\nrule clean:")]

    for legacy_root in (
        'OUTDIR + "/window_compare/',
        'OUTDIR + "/knn_empirical_pvalue/',
        'OUTDIR + "/meltron/',
        'OUTDIR + "/is_change/',
        'OUTDIR + "/figures/is_change/',
    ):
        assert legacy_root not in rule_all

    for scale_qualified_target in (
        'SCALE_ROOT + "/windows/windows.parquet"',
        'SCALE_ROOT + "/logs/gap_qc/{sample}.done"',
        'SCALE_ROOT + "/logs/windows/{sample}.done"',
        'SCALE_ROOT + "/window_compare/{comparison}/{comparison}_covariates.tsv"',
        'SCALE_ROOT + "/window_compare/{comparison}/{comparison}_gaussian_metrics.tsv"',
        'SCALE_ROOT + "/window_compare/{comparison}/{comparison}_tw.tsv"',
        'SCALE_ROOT + "/window_compare/{comparison}/{comparison}_window_pair_metrics.tsv"',
        'SCALE_ROOT + "/knn_empirical_pvalue/{comparison}/{comparison}_knn_enriched.tsv"',
        'SCALE_ROOT + "/logs/knn_significant_heatmap/{comparison}/{chrom}.done"',
        'SCALE_ROOT + "/meltron/{comparison}/{comparison}_meltron_merged.tsv"',
    ):
        assert scale_qualified_target in rule_all
    assert rule_all.count("window_bins=WINDOW_BINS") == 10
    assert 'MULTISCALE_ROOT + "/{comparison}/candidate_union.tsv"' in rule_all
    assert 'MULTISCALE_ROOT + "/{comparison}/is_change/{comparison}_is_change_loci.tsv"' in rule_all
    assert 'MULTISCALE_ROOT + "/{comparison}/is_change/figures/manifest.tsv"' in rule_all


def test_scripts_receive_the_active_configfile_not_the_repository_default():
    knn_rules = (REPO / "rules" / "knn_empirical_pvalue.smk").read_text()
    is_change_rules = (REPO / "rules" / "is_change.smk").read_text()

    assert 'config = ACTIVE_CONFIGFILE,' in knn_rules
    assert '--config "{params.config}"' in knn_rules
    assert is_change_rules.count("config=ACTIVE_CONFIGFILE,") == 3
    assert is_change_rules.count('--config "{params.config}"') == 3
    assert "config/config.yaml" not in knn_rules
    assert "config/config.yaml" not in is_change_rules

def test_knn_execution_honors_active_10kb_configfile(tmp_path):
    """A real KNN rule execution must pass the CLI-selected config to its script."""
    config_path, outdir = make_config(tmp_path, 10000)
    config = yaml.safe_load(config_path.read_text())
    config["knn_empirical_pvalue"]["knn_k"] = 2
    config["knn_empirical_pvalue"]["n_jobs"] = 1
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    background = tmp_path / "background_r10000_w25.tsv"
    with background.open("a") as handle:
        handle.write(
            "background_a_vs_background_b\tchr1\t50000\t300000\t"
            "1.1\t1.1\t1.1\t1.1\t0.8\t0.7\t0.2\t0.6\n"
        )

    merged = outdir / (
        "window_scales/r10000/w25/window_compare/H1398_vs_L134/"
        "H1398_vs_L134_window_pair_metrics.tsv"
    )
    merged.parent.mkdir(parents=True)
    merged.write_text(
        "chrom\twindow_id\twindow_start\twindow_end\tcomparison\tEC\tIB\tSP\t"
        "MSE\tSSIM\tSpearman\ttw_pvalue\tNCC\n"
        "chr1\t1\t0\t250000\tH1398_vs_L134\t1.0\t1.0\t1.0\t1.0\t0.9\t0.8\t0.1\t0.7\n"
        "chr1\t2\t50000\t300000\tH1398_vs_L134\t1.2\t1.2\t1.2\t1.2\t0.8\t0.7\t0.2\t0.6\n"
    )
    target = outdir / (
        "window_scales/r10000/w25/knn_empirical_pvalue/H1398_vs_L134/"
        "H1398_vs_L134_knn_enriched.tsv"
    )


    snakefile = tmp_path / "knn_execution.smk"
    snakefile.write_text(
        f"""import sys
configfile: {str(REPO / "config" / "config.yaml")!r}
configfile: {str(REPO / "config" / "chroms.yaml")!r}

SCRIPT_DIR = {str(REPO)!r}
sys.path.insert(0, SCRIPT_DIR + "/scripts")
from window_analysis_config import parse_window_analysis_settings

WINDOW_ANALYSIS = parse_window_analysis_settings(config, check_files=False)
OUTDIR = config["outdir"]
RESOLUTION = WINDOW_ANALYSIS.resolution_bp
SCALE_ROOT = OUTDIR + "/window_scales/r{{}}".format(RESOLUTION) + "/w{{window_bins}}"
_CLI_CONFIGFILES = tuple(workflow.overwrite_configfiles)
ACTIVE_CONFIGFILE = str(
    _CLI_CONFIGFILES[-1] if _CLI_CONFIGFILES else workflow.configfiles[0]
)

def get_window_scale(window_bins):
    return WINDOW_ANALYSIS.scale(int(window_bins))

def get_background_preset(window_bins):
    return WINDOW_ANALYSIS.backgrounds_by_window_bins[int(window_bins)]

include: {str(REPO / "rules" / "knn_empirical_pvalue.smk")!r}

rule all:
    input:
        {str(target)!r},
"""
    )

    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(snakefile),
            "--cores",
            "1",
            "--configfile",
            str(config_path),
            "--",
            str(target),
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )
    report = result.stdout + result.stderr

    assert result.returncode == 0, report
    assert target.exists()
    log = outdir / (
        "window_scales/r10000/w25/logs/knn_empirical_pvalue/H1398_vs_L134/"
        "knn_empirical_pvalue.log"
    )
    assert "KNN k: 2" in log.read_text()
    assert "resolution_bp\twindow_bins" in target.read_text()
    assert "10000\t25" in target.read_text()
