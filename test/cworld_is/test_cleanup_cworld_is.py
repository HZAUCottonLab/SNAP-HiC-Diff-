from pathlib import Path
import subprocess
import sys

import pytest


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "cleanup_cworld_is.py"


def populate(results, sample):
    (results / "cworld_insulation" / sample).mkdir(parents=True)
    (results / "IS" / sample).mkdir(parents=True)
    (results / "logs" / "compute_IS").mkdir(parents=True, exist_ok=True)
    (results / "benchmark" / "compute_IS").mkdir(parents=True, exist_ok=True)
    (results / "cworld_insulation" / sample / "a.bedGraph").write_text("x")
    (results / "IS" / sample / f"{sample}_IS.tsv").write_text("x")
    (results / "logs" / "compute_IS" / f"{sample}__A01.log").write_text("x")
    (results / "benchmark" / "compute_IS" / f"{sample}__A01.tsv").write_text("x")


def run_cleanup(results, *extra):
    return run_cleanup_for_sample(results, "S1", *extra)


def run_cleanup_for_sample(results, sample, *extra):
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--outdir",
            str(results),
            "--samples",
            sample,
            *extra,
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def test_cleanup_defaults_to_dry_run(tmp_path):
    results = tmp_path / "results"
    populate(results, "S1")

    result = run_cleanup(results)

    assert result.returncode == 0, result.stderr
    assert "[DRY-RUN]" in result.stdout
    assert (results / "cworld_insulation" / "S1").exists()
    assert (results / "IS" / "S1" / "S1_IS.tsv").exists()


def test_execute_removes_only_selected_scope(tmp_path):
    results = tmp_path / "results"
    populate(results, "S1")
    populate(results, "S2")
    protected = results / "balanced" / "S1" / "S1__HC04_A01_iced.cool"
    protected.parent.mkdir(parents=True)
    protected.write_text("must remain")

    result = run_cleanup(results, "--execute")

    assert result.returncode == 0, result.stderr
    assert not (results / "cworld_insulation" / "S1").exists()
    assert not (results / "IS" / "S1").exists()
    assert not list((results / "logs" / "compute_IS").glob("S1*"))
    assert not list((results / "benchmark" / "compute_IS").glob("S1*"))
    assert (results / "cworld_insulation" / "S2").exists()
    assert (results / "IS" / "S2").exists()
    assert protected.read_text() == "must remain"


def test_execute_keeps_shared_prefix_sample_artifacts(tmp_path):
    results = tmp_path / "results"
    populate(results, "S1")
    populate(results, "S10")

    result = run_cleanup(results, "--execute")

    assert result.returncode == 0, result.stderr
    assert not (results / "cworld_insulation" / "S1").exists()
    assert not (results / "IS" / "S1").exists()
    assert not list((results / "logs" / "compute_IS").glob("S1__*"))
    assert not list((results / "benchmark" / "compute_IS").glob("S1__*"))
    assert (results / "cworld_insulation" / "S10" / "a.bedGraph").is_file()
    assert (results / "IS" / "S10" / "S10_IS.tsv").is_file()
    assert (results / "logs" / "compute_IS" / "S10__A01.log").is_file()
    assert (results / "benchmark" / "compute_IS" / "S10__A01.tsv").is_file()


@pytest.mark.parametrize(
    "sample",
    ["*", "?", "[abc]", "S 1", "S1/", "S1\\", ".", "..", "S.1"],
)
def test_execute_rejects_unsafe_sample_before_deletion(tmp_path, sample):
    results = tmp_path / "results"
    populate(results, "S1")
    populate(results, "S2")

    result = run_cleanup_for_sample(results, sample, "--execute")

    assert result.returncode != 0
    for safe_sample in ("S1", "S2"):
        assert (results / "cworld_insulation" / safe_sample / "a.bedGraph").is_file()
        assert (results / "IS" / safe_sample / f"{safe_sample}_IS.tsv").is_file()
        assert (
            results / "logs" / "compute_IS" / f"{safe_sample}__A01.log"
        ).is_file()
        assert (
            results / "benchmark" / "compute_IS" / f"{safe_sample}__A01.tsv"
        ).is_file()
