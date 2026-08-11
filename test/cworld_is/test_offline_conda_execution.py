from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
CWORLD_ENV = REPO / "env" / "cworld.yaml"
MELTRON_ENV = REPO / "env" / "meltron.yaml"
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
        "r-base",
        "perl=5.32",
        "perl-gd>=2.57",
        "bedtools",
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


def test_meltron_environment_declares_python_dataframe_runtime_and_validates_it():
    dependencies = dependency_lines(MELTRON_ENV)
    script = CONDA_HELPER.read_text()
    meltron_case_start = script.index("        env/meltron.yaml)")
    meltron_case = script[
        meltron_case_start : script.index("            ;;", meltron_case_start)
    ]

    assert {"python>=3.10", "numpy", "pandas"} <= dependencies
    for package in ("numpy", "pandas"):
        assert (
            f'validate_python_package "$meltron_prefix/bin/python" {package}'
            in meltron_case
        )


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
    assert "write_snakemake_cache_metadata()" in script
    assert '"$CONDA_PREFIX_DIR/${hash}_.yaml"' in script
    assert '"$CONDA_PREFIX_DIR/${hash}_.env_setup_done"' in script
    assert 'cp "$ROOT/$spec" "$cache_yaml"' in script
    assert '-I "$ROOT/env/cworld-dekker/lib" -Mcworld::dekker' in script
    assert '[[ -x "$cworld_prefix/bin/bedtools" ]]' in script
    assert '"$cworld_prefix/bin/bedtools" --version' in script
    assert '[[ -x "$cworld_prefix/bin/Rscript" ]]' in script
    assert '"$cworld_prefix/bin/Rscript" --version' in script
    assert '"$snap_prefix/bin/hictk" --version' in script
    assert "Rscript" in script
    assert "source \"$READY_MARKER\"" not in script


def test_python_package_validation_reads_metadata_without_importing_modules():
    script = CONDA_HELPER.read_text()

    assert "from importlib.metadata import version" in script
    for module in ("cooler", "cooltools", "numpy", "pandas"):
        assert f"'import {module}'" not in script

    for package in ("cooler", "cooltools", "numpy", "pandas"):
        assert f'validate_python_package "$snap_prefix/bin/python" {package}' in script
    for package in ("cooler", "numpy", "pandas"):
        assert f'validate_python_package "$cworld_prefix/bin/python" {package}' in script


def test_conda_helper_requires_nonempty_conda_history():
    script = CONDA_HELPER.read_text()

    assert '[[ -s "$prefix/conda-meta/history" ]]' in script
    assert '[[ -f "$prefix/conda-meta/history" ]]' not in script


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

    wrapper_start = script.index("    printf '#!/usr/bin/env bash\\n'")
    wrapper_end = script.index('} > "$_wrapper"', wrapper_start)
    wrapper = script[wrapper_start:wrapper_end]

    assert "printf 'set -euo pipefail\\n'" in wrapper
    assert "printf 'export CONDA_OFFLINE=true\\n'" in wrapper
    assert 'cat "$_jobscript"' in wrapper
    assert wrapper.index("#!/usr/bin/env bash") < wrapper.index("set -euo pipefail") < wrapper.index("export CONDA_OFFLINE=true") < wrapper.index('cat "$_jobscript"')


def test_cworld_rscript_checks_follow_cworld_prefix_declaration():
    script = CONDA_HELPER.read_text()

    cworld_case = script.index("env/cworld.yaml)")
    cworld_prefix = script.index('local cworld_prefix="$prefix"', cworld_case)
    cworld_case_end = script.index("            ;;", cworld_case)

    for statement in (
        '[[ -x "$cworld_prefix/bin/Rscript" ]]',
        '"$cworld_prefix/bin/Rscript" --version',
    ):
        assert cworld_prefix < script.index(statement) < cworld_case_end


def test_conda_helper_checks_meltron_r_package_metadata_without_loading_packages():
    script = CONDA_HELPER.read_text()

    for package in ("argparser", "data.table", "DescTools"):
        assert f'"$meltron_prefix/lib/R/library/{package}/DESCRIPTION"' in script

    assert 'library(argparser)' not in script
    assert 'library(data.table)' not in script
    assert 'library(DescTools)' not in script


def test_meltron_gitlink_declares_its_official_source_mapping():
    gitmodules = REPO / ".gitmodules"

    assert gitmodules.is_file(), "Meltron gitlink requires a tracked .gitmodules mapping"
    contents = gitmodules.read_text()
    assert '[submodule "env/Meltron"]' in contents
    assert "\tpath = env/Meltron" in contents
    assert "\turl = https://github.com/pombo-lab/Meltron" in contents


def test_launcher_preflights_meltron_source_before_calling_snakemake():
    script = RUN_SNAKEMAKE.read_text()
    source_file = "${WORKDIR}/env/Meltron/command_line_apps/Meltron.R"

    assert source_file in script
    assert "git submodule update --init --recursive" in script
    assert script.index(source_file) < script.index("snakemake " + "\\")


def test_launcher_uses_requested_scheduler_limit_and_all_target():
    script = RUN_SNAKEMAKE.read_text()

    assert "    -j 1000" in script
    assert script.rstrip().endswith("    all")
