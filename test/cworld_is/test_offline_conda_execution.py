from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
CWORLD_ENV = REPO / "env" / "cworld.yaml"
CONDA_HELPER = REPO / "scripts" / "conda_envs.sh"
RUN_SNAKEMAKE = REPO / "run_snakemake.sh"
DSUB_SUBMIT = REPO / "config" / "dsub_submit.sh"


def dependency_lines(path):
    return {
        line.strip()[2:]
        for line in path.read_text().splitlines()
        if line.strip().startswith("- ")
    }


def test_cworld_environment_uses_only_the_required_arm_runtime():
    dependencies = dependency_lines(CWORLD_ENV)

    assert {
        "python>=3.10",
        "cooler>=0.10",
        "numpy",
        "pandas",
        "perl=5.32",
        "perl-gd>=2.57",
        "zlib",
    } <= dependencies
    assert not {
        dependency
        for dependency in dependencies
        if dependency.startswith("perl-file-slurp")
        or dependency.startswith("perl-statistics-basic")
        or dependency.startswith("perl-math-statistics")
        or dependency.startswith("perl-set-intervaltree")
        or dependency.startswith("perl-list-moreutils")
        or dependency.startswith("perl-list-util")
        or dependency.startswith("perl-data-dumper")
        or dependency.startswith("perl-file-copy-recursive")
        or dependency.startswith("perl-text-csv")
    }


def test_conda_helper_declares_all_specs_marker_and_runtime_checks():
    script = CONDA_HELPER.read_text()

    for spec in (
        "env/SNAP-HiC-Diff.yaml",
        "env/cworld.yaml",
        "env/meltron.yaml",
    ):
        assert spec in script
    assert 'ACTION="${1:-}"' in script
    assert "offline.ready.tsv" in script
    assert "--conda-create-envs-only" not in script
    assert 'MAMBA_BIN="${MAMBA_BIN:-mamba}"' in script
    assert "environment_hash()" in script
    assert "md5sum" in script
    assert '"$MAMBA_BIN" env create --yes --prefix "$prefix"' in script
    assert "conda-meta/history" in script
    assert '-I "$ROOT/env/cworld-dekker/lib" -Mcworld::dekker' in script
    assert '"$snap_prefix/bin/hictk" --version' in script
    assert "Rscript" in script
    assert "source \"$READY_MARKER\"" not in script


def test_launcher_verifies_marker_and_forces_offline_conda_before_snakemake():
    script = RUN_SNAKEMAKE.read_text()

    verify_pos = script.index('scripts/conda_envs.sh" verify')
    snakemake_pos = script.index("    --executor cluster-generic")
    assert verify_pos < snakemake_pos
    assert "export CONDA_OFFLINE=true" in script
    assert 'export CONDA_PREFIX_DIR="${WORKDIR}/.snakemake/conda"' in script
    assert '--conda-prefix "${CONDA_PREFIX_DIR}"' in script
    assert "--conda-create-envs-only" not in script


def test_dsub_wrapper_injects_offline_mode_into_generated_worker_script():
    script = DSUB_SUBMIT.read_text()

    assert "printf '#!/usr/bin/env bash\\n'" in script
    assert "printf 'export CONDA_OFFLINE=true\\n'" in script
    assert 'cat "$_jobscript"' in script
    assert script.index("export CONDA_OFFLINE=true") < script.index('cat "$_jobscript"')
