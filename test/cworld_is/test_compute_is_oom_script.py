from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "compute_is_oom_20260609"
    / "run_test.sh"
)


def test_smoke_script_allows_data_outside_code_checkout():
    script = SCRIPT.read_text()

    assert 'DATA_PROJ="${DATA_PROJ:-$PROJ}"' in script
    assert 'ENV_PROJ="${ENV_PROJ:-$PROJ}"' in script
    assert 'export CONDA_ENVS_PATH="${ENV_PROJ}/.snakemake/conda"' in script
    assert script.count('"$DATA_PROJ/results/balanced/L021/L021__HC04_D08_iced.cool"') == 2
