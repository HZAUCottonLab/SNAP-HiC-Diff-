from pathlib import Path
import sys

import numpy as np


REPO = Path(__file__).resolve().parents[2]
CODE_DIR = REPO / "scripts" / "python3_3dnetmod"
sys.path.insert(0, str(CODE_DIR))
import genlouvain


def test_genlouvain_tolerates_another_worker_creating_output(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    original_makedirs = genlouvain.os.makedirs

    monkeypatch.setattr(
        genlouvain.os.path,
        "isdir",
        lambda path: False if path == "output" else Path(path).is_dir(),
    )

    def create_only_when_race_safe(path, exist_ok=False):
        if not exist_ok:
            raise FileExistsError(path)
        original_makedirs(path, exist_ok=True)

    monkeypatch.setattr(genlouvain.os, "makedirs", create_only_when_race_safe)

    genlouvain.genlouvain(np.zeros((2, 2)), seed=1)

    assert (tmp_path / "output" / "gen_louvain_seed.txt").is_file()
