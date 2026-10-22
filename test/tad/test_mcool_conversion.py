import importlib.util
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import yaml


REPO = Path(__file__).resolve().parents[2]


def load_converter():
    path = REPO / "scripts" / "convert_cool_to_mcool.py"
    spec = importlib.util.spec_from_file_location("convert_cool_to_mcool", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_conversion_builds_balanced_4_8_10_20_100kb_mcool():
    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    assert config["convert"]["resolutions"] == [4000, 8000, 10000, 20000, 100000]

    converter = load_converter()
    with patch.object(converter.subprocess, "run") as run:
        converter.run_convert(
            "input.cool",
            "output.mcool",
            config["convert"]["resolutions"],
            20,
            StringIO(),
        )

    run.assert_called_once_with(
        [
            "cooler", "zoomify", "input.cool",
            "-o", "output.mcool",
            "-r", "4000,8000,10000,20000,100000",
            "--balance",
            "--nproc", "20",
        ],
        check=True,
    )
