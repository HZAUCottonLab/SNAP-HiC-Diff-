import os
from copy import deepcopy
from pathlib import Path
import subprocess

import yaml


REPO = Path(__file__).resolve().parents[2]


def test_boundary_annotation_rules_and_all_outputs_are_wired():
    snakefile = (REPO / "Snakefile").read_text()
    rules = (REPO / "rules" / "is_boundary_annotation.smk").read_text()
    for name in (
        "rule unify_tad_boundaries:",
        "rule call_tad_boundary_events:",
        "rule annotate_is_change_boundaries:",
        "rule permute_is_boundary_association:",
    ):
        assert name in rules
    assert 'include: "rules/tad.smk"' in snakefile
    assert 'include: "rules/is_boundary_annotation.smk"' in snakefile
    assert snakefile.index('include: "rules/tad.smk"') < snakefile.index(
        'include: "rules/is_boundary_annotation.smk"'
    )
    assert "IS_BOUNDARY_ALL_OUTPUTS" in snakefile
    assert "PERMUTATION_ALL_OUTPUTS" in snakefile


def test_boundary_rule_shells_quote_paths_and_prepare_benchmarks():
    rules = (REPO / "rules" / "is_boundary_annotation.smk").read_text()

    assert rules.count('benchmark_dir=OUTDIR + "/benchmark/is_boundary",') == 4
    assert rules.count('>> {log:q} 2>&1') == 4
    for argument in (
        "--sample {wildcards.sample:q}",
        "--raw-multires {input.raw_multires:q}",
        "--balanced-multires {input.balanced_multires:q}",
        "--raw-resolution-boundaries {input.raw_resolution:q}",
        "--balanced-resolution-boundaries {input.balanced_resolution:q}",
        "--boundaries-a {input.boundaries_a:q}",
        "--boundaries-b {input.boundaries_b:q}",
        "--events {input.events:q}",
        "--primary {input.primary:q}",
        "--preliminary {input.preliminary:q}",
        "--loci {input.loci:q}",
        "--output {output.permutation:q}",
    ):
        assert argument in rules


def _dry_run(target, cache_dir, configfile="config/config.yaml"):
    return subprocess.run(
        ["snakemake", "-n", "-p", "--configfile", configfile, "--", target],
        cwd=REPO,
        text=True,
        capture_output=True,
        env={
            **os.environ,
            "PYTHONDONTWRITEBYTECODE": "1",
            "XDG_CACHE_HOME": str(cache_dir),
        },
    )


def test_annotated_primary_dag_has_boundary_rules_but_not_permutation(tmp_path):
    target = (
        "results/multiscale/r4000/H1398_vs_L134/is_change/"
        "H1398_vs_L134_is_change_loci.boundary_annotated.tsv"
    )
    completed = _dry_run(target, tmp_path / "cache")
    plan = completed.stdout + completed.stderr

    assert completed.returncode == 0, plan
    for rule in (
        "unify_tad_boundaries",
        "call_tad_boundary_events",
        "annotate_is_change_boundaries",
    ):
        assert rule in plan
    assert "H1398" in plan
    assert "L134" in plan
    assert "permute_is_boundary_association" not in plan


def test_plot_rule_reverts_to_unannotated_primary_when_disabled(tmp_path):
    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    config = deepcopy(config)
    config["is_boundary_annotation"]["enabled"] = False
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))
    target = (
        "results/multiscale/r4000/H1398_vs_L134/is_change/figures/manifest.tsv"
    )
    completed = _dry_run(target, tmp_path / "cache", str(config_path))
    plan = completed.stdout + completed.stderr

    assert completed.returncode == 0, plan
    assert "annotate_is_change_boundaries" not in plan
    assert "--boundary-events" not in plan


def test_permutation_dag_is_targeted_when_enabled(tmp_path):
    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    config = deepcopy(config)
    config["is_boundary_annotation"]["permutation"]["enabled"] = True
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))
    target = "results/tad_boundary_events/H1398_vs_L134.permutation.tsv"
    completed = _dry_run(target, tmp_path / "cache", str(config_path))
    plan = completed.stdout + completed.stderr

    assert completed.returncode == 0, plan
    for rule in (
        "unify_tad_boundaries",
        "call_tad_boundary_events",
        "annotate_is_change_boundaries",
        "permute_is_boundary_association",
    ):
        assert rule in plan
