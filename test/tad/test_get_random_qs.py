from pathlib import Path
import sys

import numpy as np


REPO = Path(__file__).resolve().parents[2]
CODE_DIR = REPO / "scripts" / "python3_3dnetmod"
sys.path.insert(0, str(CODE_DIR))
import get_random_Qs


def test_get_random_qs_reloads_existing_csv_in_python3(tmp_path, monkeypatch):
    monkeypatch.setattr(get_random_Qs, "OUTDIR", str(tmp_path / "run"))
    network_dir = tmp_path / "run" / "GPS" / "random_networks"
    network_dir.mkdir(parents=True)
    filename = network_dir / "sample_chrregion_random_network_0_tag.csv"
    np.savetxt(filename, np.array([[1.0, 2.0], [3.0, 4.0]]), delimiter=",")

    networks = get_random_Qs.get_random_Qs(
        {"num_rand": 1}, {}, "sample", "region", "tag", "chr"
    )

    assert np.array_equal(networks[0], np.array([[1.0, 2.0], [3.0, 4.0]]))
