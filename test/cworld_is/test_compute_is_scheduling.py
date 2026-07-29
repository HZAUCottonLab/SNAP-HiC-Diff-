from pathlib import Path
import re
import subprocess


REPO = Path(__file__).resolve().parents[2]
RULE = REPO / "rules" / "compute_IS.smk"
SAMPLES = [
    "H1394",
    "H1395",
    "H1396",
    "H1397",
    "H1398",
    "H1406",
    "H1409",
    "H1410",
    "H1436",
    "L021",
]
TARGETS = [f"results/IS/{sample}/{sample}_IS.tsv" for sample in SAMPLES]


def test_compute_is_schedules_one_cworld_job_per_sample_chromosome():
    result = subprocess.run(
        [
            "snakemake",
            "-n",
            "-p",
            "--configfile",
            "config/config.yaml",
            "--",
            *TARGETS,
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )
    plan = result.stdout + result.stderr

    assert result.returncode == 0, plan
    assert re.search(r"^\s*compute_IS_chrom\s+260\s*$", plan, re.MULTILINE), plan
    assert re.search(r"^\s*merge_IS\s+10\s*$", plan, re.MULTILINE), plan
    assert "checkpoint compute_IS" not in RULE.read_text()
    assert 'config["meltron"]["is_distances"]' not in RULE.read_text()
    assert plan.count("matrix2insulation.pl") == 260, plan


def test_compute_is_rule_uses_chromosome_outputs_and_configured_cpus():
    rule = RULE.read_text()

    assert "rule compute_IS_chrom:" in rule
    assert "rule merge_IS:" in rule
    assert '"--expected-bin-size"' in rule
    assert "tmp/cworld_matrix" in rule
    assert "CWORLD_IS.chrom_job_threads" in rule
    assert "CWORLD_IS.merge_job_threads" in rule
    assert 'MANIFEST_TMP="$MANIFEST.tmp.$$"' in rule
    assert "find" in rule and "-newer" in rule
