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
    assert validator.index('"cworld_artifacts.txt"') < validator.index('if not matched_bedgraph:')
