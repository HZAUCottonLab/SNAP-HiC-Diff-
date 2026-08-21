from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from build_windows import build_windows


def test_build_windows_records_the_active_scale_on_every_window():
    windows = build_windows(
        {"chr1": 140_000},
        window_size=100_000,
        step_size=20_000,
        resolution=4_000,
        window_bins=25,
        step_bins=5,
    )

    assert not windows.empty
    assert windows[
        [
            "resolution_bp",
            "window_bins",
            "window_size_bp",
            "step_bins",
            "step_size_bp",
        ]
    ].drop_duplicates().to_dict("records") == [
        {
            "resolution_bp": 4_000,
            "window_bins": 25,
            "window_size_bp": 100_000,
            "step_bins": 5,
            "step_size_bp": 20_000,
        }
    ]


def test_window_preprocessing_paths_are_scale_qualified_without_moving_shared_inputs():
    snakefile = (REPO / "Snakefile").read_text()
    assert 'SCALE_ROOT = "{}/window_scales/r{}/w{{window_bins}}".format(OUTDIR, RESOLUTION)' in snakefile

    build_windows_rule = (REPO / "rules" / "build_windows.smk").read_text()
    gap_qc_rule = (REPO / "rules" / "compute_gap_qc.smk").read_text()
    sliding_rule = (REPO / "rules" / "sliding_window.smk").read_text()

    for rule in (build_windows_rule, gap_qc_rule, sliding_rule):
        assert "SCALE_ROOT" in rule
        assert "window_bins" in rule

    assert gap_qc_rule.count('window_bins="{{window_bins}}"') == 2
    assert sliding_rule.count('window_bins="{{window_bins}}"') == 2
    assert "shared_outdir = OUTDIR" in gap_qc_rule
    assert "{params.shared_outdir}/balanced/" in gap_qc_rule
    assert "shared_outdir = OUTDIR" in sliding_rule
    assert "{params.shared_outdir}/oe/" in sliding_rule
    assert "window_bins=w.window_bins" in sliding_rule
