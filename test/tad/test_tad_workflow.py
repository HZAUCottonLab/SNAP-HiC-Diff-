from pathlib import Path
import subprocess


REPO = Path(__file__).resolve().parents[2]


def test_rule_contract_contains_all_signal_resolution_stages():
    rules = (REPO / "rules" / "tad.smk").read_text()
    snakefile = (REPO / "Snakefile").read_text()
    for rule in (
        "rule prepare_tad:",
        "rule call_tad:",
        "rule extract_tad_boundaries:",
        "rule merge_tad_boundaries:",
    ):
        assert rule in rules
    assert rules.count('../env/tad.yaml') == 4
    assert 'include: "rules/tad.smk"' in snakefile
    for token in (
        "tad_bed_path(sample, resolution, signal)",
        "tad_bedpe_path(sample, resolution, signal)",
        "tad_boundary_path(sample, resolution, signal)",
        "tad_merged_boundary_path(sample, signal)",
    ):
        assert token in snakefile


def test_tad_rules_follow_wildcard_shell_standard():
    rules = (REPO / "rules" / "tad.smk").read_text()
    snakefile = (REPO / "Snakefile").read_text()

    assert 'f"' not in rules
    assert "f'" not in rules
    assert 'f"' not in snakefile
    assert "f'" not in snakefile
    for parameter in ("resolution", "run_dir", "settings"):
        assert f"{parameter}=lambda wc:" not in rules

    assert rules.count('SAMPLE="{wildcards.sample}"') == 4
    assert rules.count('SIGNAL="{wildcards.signal}"') == 4
    assert rules.count('RESOLUTION_LABEL="{wildcards.resolution_label}"') == 3


def test_two_signal_three_resolution_dag_builds():
    completed = subprocess.run(
        [
            "snakemake", "-n", "-p",
            "--configfile", "config/config.yaml",
            "--forcerun",
            "prepare_tad",
            "call_tad",
            "extract_tad_boundaries",
            "merge_tad_boundaries",
            "--",
            "results/tad_boundaries/H1394.multires.raw.boundaries.bed",
            "results/tad_boundaries/H1394.multires.balanced.boundaries.bed",
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
    )
    plan = completed.stdout + completed.stderr
    assert completed.returncode == 0, plan
    for rule in (
        "prepare_tad", "call_tad",
        "extract_tad_boundaries", "merge_tad_boundaries",
    ):
        assert rule in plan
    for signal in ("raw", "balanced"):
        for resolution in ("4kb", "8kb", "20kb"):
            assert f"H1394.{resolution}.{signal}.3dnetmod.bed" in plan
