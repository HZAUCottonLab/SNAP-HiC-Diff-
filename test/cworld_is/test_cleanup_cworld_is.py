from pathlib import Path
import subprocess
import sys


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
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--outdir", str(results), "--samples", "S1", *extra],
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
