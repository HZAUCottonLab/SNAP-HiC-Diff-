from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "compute_is_oom_20260609"
    / "run_test.sh"
)
VALIDATOR = SCRIPT.parents[2] / "scripts" / "validate_cworld_is.py"


def test_smoke_script_allows_data_outside_code_checkout():
    script = SCRIPT.read_text()

    assert 'DATA_PROJ="${DATA_PROJ:-$PROJ}"' in script
    assert 'ENV_PROJ="${ENV_PROJ:-$PROJ}"' in script
    assert 'export CONDA_ENVS_PATH="${ENV_PROJ}/.snakemake/conda"' in script
    assert script.count('"$DATA_PROJ/results/balanced/L021/L021__HC04_D08_iced.cool"') == 2


def test_smoke_uses_the_prebuilt_cworld_prefix_for_perl_and_bedtools():
    script = SCRIPT.read_text()

    assert 'CONDA_OFFLINE=true CONDA_PREFIX_DIR="$CONDA_ENVS_PATH" bash "$PROJ/scripts/conda_envs.sh" verify' in script
    assert 'CWORLD_PREFIX="$(awk -F \'\\t\' \'$1 == "env/cworld.yaml" {print $3}\' "$CONDA_ENVS_PATH/offline.ready.tsv")"' in script
    assert 'export PATH="$CWORLD_PREFIX/bin:/usr/bin:/bin"' in script
    assert '"$CWORLD_PREFIX/bin/bedtools" --version' in script
    assert '"$CWORLD_PREFIX/bin/Rscript" --version' in script


def test_smoke_retains_cworld_diagnostics():
    script = SCRIPT.read_text()
    validator = VALIDATOR.read_text()

    assert 'BASE_DIR="$PROJ/tmp/cworld_smoke/L021_HC04_D08_is32000"' in script
    assert 'RUN_DIR="$(mktemp -d "$BASE_DIR/run_XXXXXX")"' in script
    assert '--outdir "$RUN_DIR"' in script
    assert 'MELTRON_TSV="$RUN_DIR/IS/L021/L021_IS.tsv"' in script
    assert '"-v",' in validator
    assert '"cworld_stdout.log"' in validator
    assert '"cworld_stderr.log"' in validator
    assert '"cworld_artifacts.txt"' in validator
    assert validator.index('"cworld_artifacts.txt"') < validator.index('if len(matched_bedgraph) != 1:')


def test_smoke_runs_rscript_after_setting_the_controlled_path():
    script = SCRIPT.read_text()

    path_setup = script.index('export PATH="$CWORLD_PREFIX/bin:/usr/bin:/bin"')
    rscript_check = script.index('"$CWORLD_PREFIX/bin/Rscript" --version')
    bedtools_check = script.index('"$CWORLD_PREFIX/bin/bedtools" --version')

    assert path_setup < rscript_check < bedtools_check
