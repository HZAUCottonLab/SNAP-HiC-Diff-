from pathlib import Path

import yaml


REPO = Path(__file__).resolve().parents[2]
SCALES = ["w25", "w50", "w75", "w100", "w125"]
BACKGROUND_ROOT = (
    "/share/home/xhhuang24/project/PanTE/background_pool/results/"
    "multiscale_4k/window_scales"
)


def read(relative):
    return (REPO / relative).read_text()


def test_config_has_exact_bin_derived_presets():
    config = yaml.safe_load(read("config/config.yaml"))
    assert "window_size" not in config
    assert "step_size" not in config
    assert config["resolution"] == 4000
    assert config["window_analysis"] == {
        "window_bins": [25, 50, 75, 100, 125],
        "step_fraction": 0.2,
    }
    root = BACKGROUND_ROOT
    assert [
        (p["resolution_bp"], p["window_bins"], p["step_bins"], p["file"])
        for p in config["knn_empirical_pvalue"]["background_presets"]
    ] == [
        (4000, 25, 5, f"{root}/100kb/window_compare/HiC_window_background.tsv"),
        (4000, 50, 10, f"{root}/200kb/window_compare/HiC_window_background.tsv"),
        (4000, 75, 15, f"{root}/300kb/window_compare/HiC_window_background.tsv"),
        (4000, 100, 20, f"{root}/400kb/window_compare/HiC_window_background.tsv"),
        (4000, 125, 25, f"{root}/500kb/window_compare/HiC_window_background.tsv"),
    ]


def test_one_scale_wildcard_dag_keeps_common_is_and_tad_targets():
    snakefile = read("Snakefile")
    assert "WINDOW_SCALES = build_window_scales(config)" in snakefile
    assert "SCALE_LABELS = tuple(WINDOW_SCALES)" in snakefile
    assert 'scale="w25|w50|w75|w100|w125"' in snakefile
    assert 'include: "rules/build_windows.smk"' not in snakefile
    assert "WINDOW_SIZE" not in snakefile
    assert "STEP_SIZE" not in snakefile
    for value in (
        "validated/manifest.validated.tsv", "logs/extract_balance",
        "logs/expected", "logs/oe", "/IS/", "SCALE_ALL_OUTPUTS",
        "MULTISCALE_ALL_OUTPUTS", "TAD_ALL_OUTPUTS",
    ):
        assert value in snakefile
    for scale in SCALES:
        assert scale in snakefile


def test_window_dependent_rules_are_scale_qualified():
    for relative in (
        "rules/compute_gap_qc.smk", "rules/sliding_window.smk",
        "rules/window_compare.smk", "rules/knn_empirical_pvalue.smk",
        "rules/knn_significant_heatmap.smk", "rules/meltron.smk",
    ):
        text = read(relative)
        assert "SCALE_ROOT" in text, relative
        assert "{scale}" in text, relative
    is_change = read("rules/is_change.smk")
    assert "MULTISCALE_ROOT" in is_change
    assert "SCALE_LABELS" in is_change
    for option in (
        "--knn-enriched {input.knn_enriched}",
        "--meltron-merged {input.meltron_merged}",
        "--out-candidates", "--out-union-tsv", "--out-union-bed",
    ):
        assert option in is_change


def test_runner_has_no_outer_scale_submission_loop():
    runner = read("run_snakemake.sh")
    assert "snakemake \\" in runner
    assert "-j 1000" in runner
    assert "for scale in" not in runner
    assert "dsub -s" not in runner


def test_minimal_fixture_expands_five_sibling_jobs_in_one_dag(tmp_path):
    import os
    import subprocess

    fixture = tmp_path / "Snakefile"
    fixture.write_text(
        """
SCALES = ["w25", "w50", "w75", "w100", "w125"]

rule all:
    input:
        expand("results/{scale}.done", scale=SCALES)

rule branch:
    output:
        "results/{scale}.done"
    shell:
        "touch {output}"
"""
    )
    result = subprocess.run(
        ["snakemake", "-n", "--printshellcmds", "--cores", "1", "-s", str(fixture)],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        env={**os.environ, "XDG_CACHE_HOME": str(tmp_path / "cache")},
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert output.count("rule branch:") == 5
    for scale in SCALES:
        assert f"scale={scale}" in output
