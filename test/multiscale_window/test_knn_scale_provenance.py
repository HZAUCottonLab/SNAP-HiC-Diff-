import importlib.util
from pathlib import Path
import subprocess
import sys
import types
import warnings

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]
KNN_SCRIPT = REPO / "scripts" / "compute_knn_empirical_pvalue.py"
PLOT_SCRIPT = REPO / "scripts" / "plot_knn_significant_heatmap.py"

SCALE_COLUMNS = [
    "resolution_bp",
    "window_bins",
    "window_size_bp",
    "step_bins",
    "step_size_bp",
]


def metric_row(window_id, window_start, window_end):
    return {
        "chrom": "chr1",
        "window_id": window_id,
        "window_start": window_start,
        "window_end": window_end,
        "comparison": "control_vs_treatment",
        "EC": 1.0 + window_id,
        "IB": 2.0 + window_id,
        "SP": 3.0 + window_id,
        "MSE": 4.0 + window_id,
        "SSIM": 0.9 - window_id / 100,
        "Spearman": 0.8 - window_id / 100,
        "tw_pvalue": 0.1 + window_id / 100,
        "NCC": 0.7 - window_id / 100,
    }


def write_knn_config(path):
    path.write_text(
        "knn_empirical_pvalue:\n"
        "  knn_k: 2\n"
        "  scaler: robust\n"
        "  n_jobs: 1\n"
        "  score_columns: [MSE, d_SSIM, d_Spearman, tw_score, d_NCC]\n"
        "  derive_scores:\n"
        "    d_SSIM: '1 - SSIM'\n"
        "    d_Spearman: '1 - Spearman'\n"
        "    tw_score: '-log10(max(tw_pvalue, 1e-30))'\n"
        "    d_NCC: '1 - NCC'\n"
    )


def run_knn(query, background, config, output, log):
    return subprocess.run(
        [
            sys.executable,
            str(KNN_SCRIPT),
            "--comparison",
            "control_vs_treatment",
            "--input",
            str(query),
            "--background",
            str(background),
            "--config",
            str(config),
            "--resolution-bp",
            "4000",
            "--window-bins",
            "25",
            "--window-size-bp",
            "100000",
            "--step-bins",
            "5",
            "--step-size-bp",
            "20000",
            "--out",
            str(output),
            "--log",
            str(log),
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def test_knn_cli_emits_constant_scale_provenance(tmp_path):
    query = tmp_path / "query.tsv"
    background = tmp_path / "background.tsv"
    config = tmp_path / "knn.yaml"
    output = tmp_path / "knn_enriched.tsv"
    log = tmp_path / "knn.log"
    pd.DataFrame(
        [metric_row(1, 0, 100000), metric_row(2, 20000, 120000)]
    ).to_csv(query, sep="\t", index=False)
    pd.DataFrame(
        [metric_row(3, 40000, 140000), metric_row(4, 60000, 160000)]
    ).to_csv(background, sep="\t", index=False)
    write_knn_config(config)

    result = run_knn(query, background, config, output, log)

    assert result.returncode == 0, result.stderr
    enriched = pd.read_csv(output, sep="\t")
    assert enriched.columns[:10].tolist() == SCALE_COLUMNS + [
        "chrom",
        "window_id",
        "window_start",
        "window_end",
        "comparison",
    ]
    assert enriched[SCALE_COLUMNS].drop_duplicates().to_dict("records") == [
        {
            "resolution_bp": 4000,
            "window_bins": 25,
            "window_size_bp": 100000,
            "step_bins": 5,
            "step_size_bp": 20000,
        }
    ]


def test_knn_cli_rejects_duplicate_scale_coordinate_identity(tmp_path):
    query = tmp_path / "query.tsv"
    background = tmp_path / "background.tsv"
    config = tmp_path / "knn.yaml"
    output = tmp_path / "knn_enriched.tsv"
    log = tmp_path / "knn.log"
    duplicate = metric_row(1, 0, 100000)
    pd.DataFrame([duplicate, duplicate]).to_csv(query, sep="\t", index=False)
    pd.DataFrame(
        [metric_row(3, 40000, 140000), metric_row(4, 60000, 160000)]
    ).to_csv(background, sep="\t", index=False)
    write_knn_config(config)

    result = run_knn(query, background, config, output, log)

    assert result.returncode != 0
    assert "duplicate query resolution/window/chrom/window coordinates" in result.stderr


def test_knn_rule_uses_scale_specific_background_audit_without_legacy_key():
    rule_text = (REPO / "rules" / "knn_empirical_pvalue.smk").read_text()

    assert "rule validate_knn_background:" in rule_text
    assert "background_validation.json" in rule_text
    assert "audit = rules.validate_knn_background.output.audit" in rule_text
    assert "get_background_preset(w.window_bins).file" in rule_text
    assert "background_file" not in rule_text


def load_plotter(monkeypatch):
    fake_cooler = types.ModuleType("cooler")
    fake_cooler.Cooler = None
    fake_h5py = types.ModuleType("h5py")
    fake_h5py.File = None
    monkeypatch.setitem(sys.modules, "cooler", fake_cooler)
    monkeypatch.setitem(sys.modules, "h5py", fake_h5py)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        spec = importlib.util.spec_from_file_location("knn_heatmap_under_test", PLOT_SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


def test_heatmap_cli_uses_explicit_shared_and_window_results_roots(monkeypatch, tmp_path):
    plotter = load_plotter(monkeypatch)
    shared_root = tmp_path / "shared-results"
    window_root = tmp_path / "window-results"
    manifest = tmp_path / "significant_windows.tsv"
    enriched = tmp_path / "knn_enriched.tsv"
    output = tmp_path / "plots"
    manifest.write_text("chrom\twindow_id\twindow_start\twindow_end\nchr1\t1\t0\t100000\n")
    enriched.write_text("chrom\twindow_id\twindow_start\twindow_end\tcomparison\n")

    cooler_paths = []
    h5_paths = []

    class FakeCooler:
        chromnames = ["chr1"]

        def __init__(self, path):
            cooler_paths.append(path)

        def matrix(self, **_kwargs):
            return self

        def fetch(self, _region):
            return self

        def toarray(self):
            return np.ones((2, 2))

    class FakeH5File:
        def __init__(self, path, _mode):
            h5_paths.append(path)

        def __enter__(self):
            return {"w000001": {"matrix": np.ones((2, 2))}}

        def __exit__(self, *_args):
            return False

    monkeypatch.setattr(plotter.cooler, "Cooler", FakeCooler)
    monkeypatch.setattr(plotter.h5py, "File", FakeH5File)
    monkeypatch.setattr(plotter, "plot_2x3_heatmap", lambda *_args: None)
    monkeypatch.setattr(plotter.os.path, "exists", lambda _path: True)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(PLOT_SCRIPT),
            "--chrom-mode",
            "--comparison",
            "control_vs_treatment",
            "--chrom",
            "chr1",
            "--enriched",
            str(enriched),
            "--manifest",
            str(manifest),
            "--out-dir",
            str(output),
            "--shared-results-root",
            str(shared_root),
            "--window-results-root",
            str(window_root),
        ],
    )

    plotter.main()

    assert cooler_paths == [
        str(shared_root / "balanced" / "control" / "control__chr1_iced.cool"),
        str(shared_root / "balanced" / "treatment" / "treatment__chr1_iced.cool"),
    ]
    assert h5_paths == [
        str(window_root / "windows" / "control" / "control__chr1_ln_oe.h5"),
        str(window_root / "windows" / "treatment" / "treatment__chr1_ln_oe.h5"),
        str(window_root / "windows" / "control" / "control__chr1_gaussian.h5"),
        str(window_root / "windows" / "treatment" / "treatment__chr1_gaussian.h5"),
    ]
