from pathlib import Path
import re
import subprocess


REPO = Path(__file__).resolve().parents[2]
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


def test_compute_is_schedules_one_cworld_job_per_requested_sample():
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
    assert re.search(r"^\s*compute_IS\s+10\s*$", plan, re.MULTILINE), plan
    assert "localrule compute_IS_bedgraphs" not in plan
    assert "merge_IS" not in plan
    assert plan.count("matrix2insulation.pl") == len(SAMPLES) == 10, plan
