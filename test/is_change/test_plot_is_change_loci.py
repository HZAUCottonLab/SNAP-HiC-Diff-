from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml
import pytest
import cooler
import matplotlib.pyplot as plt


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "plot_is_change_loci.py"
SCALES = (32000, 60000, 92000, 120000, 152000, 180000, 212000, 240000, 272000, 300000)
sys.path.insert(0, str(REPO / "scripts"))
import plot_is_change_loci as plotter


def write_config(root):
    config = {
        "meltron": {"cworld_is": {
            "source_resolution_bp": 4000, "is_distances_bp": list(SCALES),
            "ids_divisor": 2, "chrom_job_threads": 1, "merge_job_threads": 1,
        }},
        "is_change": {
            "raw_knn_pvalue_threshold": 0.05, "meltron_score_threshold": 5.0,
            "merge_gap_bp": 0, "context_bp": 300000, "local_radius_bins": 1,
            "min_local_valid_fraction": 2 / 3, "flank_fraction_of_scale": 0.25,
            "min_flank_bins": 2, "min_flank_valid_fraction": 0.80,
            "max_missing_run_bins": 1, "high_percentile": 0.95,
            "broad_percentile": 0.90, "spatial_cluster_tolerance_bp": 12000,
            "min_primary_scale_bands": 2,
            "scale_bands": {
                "small": [32000, 60000, 92000],
                "medium": [120000, 152000, 180000],
                "large": [212000, 240000, 272000, 300000],
            }, "diagnostic_plot_top_n": 1,
        },
    }
    path = root / "config.yaml"
    path.write_text(yaml.safe_dump(config))
    return path


def settings_from_config(path):
    config = yaml.safe_load(path.read_text())
    return plotter.parse_is_change_settings(
        config, plotter.parse_cworld_is_settings(config)
    )


def region_row(start=320_000, end=324_000):
    return pd.Series(
        {"region_id": 1, "chrom": "chr1", "region_start": start, "region_end": end}
    )


def write_tracks(root, nan_weight_bin=None):
    starts = np.arange(160) * 4000
    rows_a, rows_b = [], []
    for scale in SCALES:
        values = np.zeros(len(starts))
        values[80:] = 1.0
        rows_a.extend(("chr1", start, scale, 0.0) for start in starts)
        rows_b.extend(("chr1", start, scale, value) for start, value in zip(starts, values))
    columns = ["chrom", "viewpoint", "IS_distance", "IS"]
    for name, rows in (("is_a.tsv", rows_a), ("is_b.tsv", rows_b)):
        pd.DataFrame(rows, columns=columns).to_csv(root / name, sep="\t", index=False)
    bins = pd.DataFrame({"chrom": "chr1", "start": starts, "end": starts + 4000, "weight": 1.0})
    if nan_weight_bin is not None:
        bins.loc[nan_weight_bin, "weight"] = np.nan
    diagonal = pd.DataFrame({"bin1_id": np.arange(len(bins)), "bin2_id": np.arange(len(bins)), "count": 20})
    adjacent = pd.DataFrame({"bin1_id": np.arange(len(bins) - 1), "bin2_id": np.arange(1, len(bins)), "count": 10})
    pixels_a = pd.concat((diagonal, adjacent), ignore_index=True)
    pixels_b = pixels_a.copy()
    pixels_b.loc[pixels_b.bin1_id + 1 == pixels_b.bin2_id, "count"] = 5
    cooler.create_cooler(str(root / "a_iced.cool"), bins, pixels_a, ordered=True)
    cooler.create_cooler(str(root / "b_iced.cool"), bins, pixels_b, ordered=True)
    return {
        "is_a": root / "is_a.tsv", "is_b": root / "is_b.tsv",
        "cool_a": root / "a_iced.cool", "cool_b": root / "b_iced.cool",
    }



def plot_inputs(tmp_path, events):
    config = write_config(tmp_path)
    paths = write_tracks(tmp_path, nan_weight_bin=10)
    settings = settings_from_config(config)
    return {
        "region": region_row(),
        "events": events,
        "track_a": plotter.read_track(paths["is_a"], settings),
        "track_b": plotter.read_track(paths["is_b"], settings),
        "bins_a": plotter.read_balanced_bins([str(paths["cool_a"])], settings.resolution_bp),
        "bins_b": plotter.read_balanced_bins([str(paths["cool_b"])], settings.resolution_bp),
        "coolers_a": plotter.open_coolers_by_chrom([str(paths["cool_a"])], "H1398"),
        "coolers_b": plotter.open_coolers_by_chrom([str(paths["cool_b"])], "L134"),
        "settings": settings,
    }


def test_region_figure_places_two_aligned_contact_heatmaps_above_is_tracks(tmp_path):
    events = pd.DataFrame([{
        "event_start": 320_000, "event_end": 324_000, "summit": 320_000,
        "is_change_type": "transition",
    }])
    inputs = plot_inputs(tmp_path, events)

    figure, heatmap_axes, is_axes = plotter.build_region_figure(
        **inputs, sample_a="H1398", sample_b="L134"
    )

    start, end = plotter.plot_interval(inputs["region"], inputs["settings"])
    assert len(heatmap_axes) == 2
    assert len(is_axes) == 10
    assert all(axis.get_xlim() == (start, end) for axis in (*heatmap_axes, *is_axes))
    assert all(tuple(axis.images[0].get_extent()) == (start, end, end, start) for axis in heatmap_axes)
    assert heatmap_axes[0].images[0].norm.vmax == heatmap_axes[1].images[0].norm.vmax
    assert "H1398" in heatmap_axes[0].get_ylabel()
    assert "L134" in heatmap_axes[1].get_ylabel()
    assert any(np.allclose(line.get_xdata(), [320_000, 320_000]) for line in heatmap_axes[0].lines)
    assert any(np.allclose(line.get_ydata(), [320_000, 320_000]) for line in heatmap_axes[1].lines)
    plt.close(figure)


def test_broad_event_has_vertical_and_horizontal_heatmap_bands(tmp_path):
    events = pd.DataFrame([{
        "event_start": 320_000, "event_end": 324_000, "summit": np.nan,
        "is_change_type": "broad_shift",
    }])
    inputs = plot_inputs(tmp_path, events)

    figure, heatmap_axes, _ = plotter.build_region_figure(
        **inputs, sample_a="H1398", sample_b="L134"
    )

    assert any(patch.get_alpha() == pytest.approx(0.22) for patch in heatmap_axes[0].patches)
    assert any(patch.get_alpha() == pytest.approx(0.22) for patch in heatmap_axes[1].patches)
    plt.close(figure)


def test_plotter_keeps_existing_png_when_target_chromosome_cooler_is_missing(tmp_path):
    events = pd.DataFrame([{
        "event_start": 320_000, "event_end": 324_000, "summit": 320_000,
        "is_change_type": "transition",
    }])
    inputs = plot_inputs(tmp_path, events)
    output = tmp_path / "plots" / "chr1_region.png"
    output.parent.mkdir()
    output.write_bytes(b"previous-valid-png")

    with pytest.raises(ValueError, match="no iced cooler is available for chr1"):
        plotter.plot_region(
            inputs["region"], inputs["events"], inputs["track_a"], inputs["track_b"],
            inputs["bins_a"], inputs["bins_b"], {}, inputs["coolers_b"], "H1398", "L134",
            inputs["settings"], output,
        )

    assert output.read_bytes() == b"previous-valid-png"

def test_plot_interval_is_4kb_aligned_and_shared_by_all_panels(tmp_path):
    settings = settings_from_config(write_config(tmp_path))
    region = region_row()

    start, end = plotter.plot_interval(region, settings)

    assert start == 20_000
    assert end == 624_000
    assert (start % 4_000, end % 4_000) == (0, 0)


def test_open_coolers_by_chrom_rejects_non_iced_path(tmp_path):
    paths = write_tracks(tmp_path)
    uniced_path = paths["cool_a"].with_name("a.cool")
    paths["cool_a"].rename(uniced_path)

    with pytest.raises(ValueError, match=r"requires _iced.cool"):
        plotter.open_coolers_by_chrom([str(uniced_path)], "sample A")

def test_balanced_contact_matrix_masks_unavailable_bins_and_uses_log1p(tmp_path):
    paths = write_tracks(tmp_path, nan_weight_bin=10)
    coolers = plotter.open_coolers_by_chrom([str(paths["cool_a"])], "sample A")

    matrix = plotter.load_log_contact_matrix(coolers, "chr1", 0, 640_000, "sample A")

    assert np.ma.isMaskedArray(matrix)
    assert np.ma.getmaskarray(matrix)[10].all()
    assert float(matrix.max()) == pytest.approx(np.log1p(20.0))


def test_shared_contact_norm_uses_pooled_finite_995_percentile():
    matrix_a = np.ma.masked_array(
        [[0.0, 1.0], [2.0, np.nan]], mask=[[False, False], [False, True]]
    )
    matrix_b = np.ma.masked_array([[3.0, 4.0], [5.0, 6.0]])

    norm = plotter.shared_contact_norm(matrix_a, matrix_b)

    assert norm.vmin == 0.0
    assert norm.vmax == pytest.approx(np.percentile([0, 1, 2, 3, 4, 5, 6], 99.5))


def test_balanced_contact_matrix_rejects_interval_beyond_chromosome(tmp_path):
    paths = write_tracks(tmp_path)
    coolers = plotter.open_coolers_by_chrom([str(paths["cool_a"])], "sample A")

    with pytest.raises(ValueError, match="invalid plotting interval"):
        plotter.load_log_contact_matrix(coolers, "chr1", 0, 644_000, "sample A")



def test_plot_manifest_contains_existing_pngs(tmp_path):
    config = write_config(tmp_path)
    write_tracks(tmp_path)
    primary = pd.DataFrame([{
        "comparison": "A__vs__B", "region_id": 1, "chrom": "chr1",
        "event_start": 320000, "event_end": 324000, "summit": 320000,
        "is_change_type": "transition", "confidence_tier": "primary",
        "support_scale_bands": "small,medium", "support_scales": "32000,120000",
        "max_amplitude_percentile": 0.99, "max_transition_percentile": 0.99,
        "max_pulse_percentile": 0.1, "broad_shift_length_bp": 0,
        "mask_status": "clean", "coverage_state_summary": "clean", "rank_within_region": 1,
    }])
    primary.to_csv(tmp_path / "primary.tsv", sep="\t", index=False)
    pd.DataFrame([{"region_id": 1, "chrom": "chr1", "region_start": 320000, "region_end": 324000}]).to_csv(
        tmp_path / "candidate_regions.tsv", sep="\t", index=False
    )

    result = subprocess.run(
        [
            sys.executable, str(SCRIPT), "--primary", str(tmp_path / "primary.tsv"),
            "--candidate-regions", str(tmp_path / "candidate_regions.tsv"),
            "--is-a", str(tmp_path / "is_a.tsv"), "--is-b", str(tmp_path / "is_b.tsv"),
            "--cool-a", str(tmp_path / "a_iced.cool"), "--cool-b", str(tmp_path / "b_iced.cool"),
            "--config", str(config), "--outdir", str(tmp_path / "plots"),
            "--manifest", str(tmp_path / "plots" / "manifest.tsv"),
        ], text=True, capture_output=True, check=False,
    )

    assert result.returncode == 0, result.stderr
    paths = pd.read_csv(tmp_path / "plots" / "manifest.tsv", sep="\t")["plot_path"]
    assert all(Path(path).is_file() for path in paths)



def test_plotter_prefers_three_band_primary_event_across_regions(tmp_path):
    config = write_config(tmp_path)
    write_tracks(tmp_path)
    pd.DataFrame([
        {
            "comparison": "A__vs__B", "region_id": 1, "chrom": "chr1",
            "event_start": 320000, "event_end": 324000, "summit": 320000,
            "is_change_type": "transition", "confidence_tier": "primary",
            "support_scale_bands": "small,medium", "support_scales": "32000,120000",
            "max_amplitude_percentile": 0.99, "max_transition_percentile": 0.99,
            "max_pulse_percentile": 0.99, "broad_shift_length_bp": 0,
            "mask_status": "clean", "coverage_state_summary": "clean", "rank_within_region": 1,
        },
        {
            "comparison": "A__vs__B", "region_id": 2, "chrom": "chr1",
            "event_start": 320000, "event_end": 324000, "summit": 320000,
            "is_change_type": "transition", "confidence_tier": "primary",
            "support_scale_bands": "small,medium,large", "support_scales": "32000,120000,212000",
            "max_amplitude_percentile": 0.90, "max_transition_percentile": 0.90,
            "max_pulse_percentile": 0.90, "broad_shift_length_bp": 0,
            "mask_status": "clean", "coverage_state_summary": "clean", "rank_within_region": 1,
        },
    ]).to_csv(tmp_path / "primary.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"region_id": 1, "chrom": "chr1", "region_start": 320000, "region_end": 324000},
        {"region_id": 2, "chrom": "chr1", "region_start": 320000, "region_end": 324000},
    ]).to_csv(tmp_path / "candidate_regions.tsv", sep="\t", index=False)

    result = subprocess.run(
        [
            sys.executable, str(SCRIPT), "--primary", str(tmp_path / "primary.tsv"),
            "--candidate-regions", str(tmp_path / "candidate_regions.tsv"),
            "--is-a", str(tmp_path / "is_a.tsv"), "--is-b", str(tmp_path / "is_b.tsv"),
            "--cool-a", str(tmp_path / "a_iced.cool"), "--cool-b", str(tmp_path / "b_iced.cool"),
            "--config", str(config), "--outdir", str(tmp_path / "plots"),
            "--manifest", str(tmp_path / "plots" / "manifest.tsv"),
        ], text=True, capture_output=True, check=False,
    )

    assert result.returncode == 0, result.stderr
    manifest = pd.read_csv(tmp_path / "plots" / "manifest.tsv", sep="\t")
    assert manifest["region_id"].tolist() == [2]

def test_plotter_rejects_preliminary_rows(tmp_path):
    config = write_config(tmp_path)
    write_tracks(tmp_path)
    pd.DataFrame([{
        "region_id": 1, "chrom": "chr1", "event_start": 320000,
        "event_end": 324000, "summit": 320000, "is_change_type": "transition",
        "confidence_tier": "single_band_preliminary",
        "comparison": "A__vs__B", "support_scale_bands": "small", "support_scales": "32000",
        "max_amplitude_percentile": 0.99, "max_transition_percentile": 0.99,
        "max_pulse_percentile": 0.1, "broad_shift_length_bp": 0,
        "mask_status": "clean", "coverage_state_summary": "clean", "rank_within_region": 1,
    }]).to_csv(tmp_path / "primary.tsv", sep="\t", index=False)
    pd.DataFrame([{"region_id": 1, "chrom": "chr1", "region_start": 320000, "region_end": 324000}]).to_csv(
        tmp_path / "candidate_regions.tsv", sep="\t", index=False
    )

    result = subprocess.run(
        [
            sys.executable, str(SCRIPT), "--primary", str(tmp_path / "primary.tsv"),
            "--candidate-regions", str(tmp_path / "candidate_regions.tsv"),
            "--is-a", str(tmp_path / "is_a.tsv"), "--is-b", str(tmp_path / "is_b.tsv"),
            "--cool-a", str(tmp_path / "a_iced.cool"), "--cool-b", str(tmp_path / "b_iced.cool"),
            "--config", str(config), "--outdir", str(tmp_path / "plots"),
            "--manifest", str(tmp_path / "plots" / "manifest.tsv"),
        ], text=True, capture_output=True, check=False,
    )

    assert result.returncode != 0
    assert "preliminary" in result.stderr


def test_plotter_rejects_off_grid_track_coordinates(tmp_path):
    config = write_config(tmp_path)
    write_tracks(tmp_path)
    track = pd.read_csv(tmp_path / "is_a.tsv", sep="\t")
    track["viewpoint"] = track["viewpoint"].astype(float)
    track.loc[0, "viewpoint"] = 0.5
    track.to_csv(tmp_path / "is_a.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "comparison": "A__vs__B", "region_id": 1, "chrom": "chr1",
        "event_start": 320000, "event_end": 324000, "summit": 320000,
        "is_change_type": "transition", "confidence_tier": "primary",
        "support_scale_bands": "small,medium", "support_scales": "32000,120000",
        "max_amplitude_percentile": 0.99, "max_transition_percentile": 0.99,
        "max_pulse_percentile": 0.1, "broad_shift_length_bp": 0,
        "mask_status": "clean", "coverage_state_summary": "clean", "rank_within_region": 1,
    }]).to_csv(tmp_path / "primary.tsv", sep="\t", index=False)
    pd.DataFrame([{"region_id": 1, "chrom": "chr1", "region_start": 320000, "region_end": 324000}]).to_csv(
        tmp_path / "candidate_regions.tsv", sep="\t", index=False
    )

    result = subprocess.run(
        [
            sys.executable, str(SCRIPT), "--primary", str(tmp_path / "primary.tsv"),
            "--candidate-regions", str(tmp_path / "candidate_regions.tsv"),
            "--is-a", str(tmp_path / "is_a.tsv"), "--is-b", str(tmp_path / "is_b.tsv"),
            "--cool-a", str(tmp_path / "a_iced.cool"), "--cool-b", str(tmp_path / "b_iced.cool"),
            "--config", str(config), "--outdir", str(tmp_path / "plots"),
            "--manifest", str(tmp_path / "plots" / "manifest.tsv"),
        ], text=True, capture_output=True, check=False,
    )

    assert result.returncode != 0
    assert "integers" in result.stderr


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        ("viewpoint", np.inf, "finite integers"),
        ("viewpoint", -4000, "non-negative"),
        ("IS_distance", 32000.5, "finite integers"),
        ("IS_distance", 40000, "unconfigured scales"),
    ],
)
def test_plotter_strictly_validates_track_coordinates_and_scales(tmp_path, column, value, message):
    sys.path.insert(0, str(REPO / "scripts"))
    import plot_is_change_loci as plotter
    from types import SimpleNamespace

    write_tracks(tmp_path)
    path = tmp_path / "is_a.tsv"
    track = pd.read_csv(path, sep="\t")
    track[column] = track[column].astype(float)
    track.loc[0, column] = value
    track.to_csv(path, sep="\t", index=False)

    settings = SimpleNamespace(resolution_bp=4000, scales_bp=SCALES)
    with pytest.raises(ValueError, match=message):
        plotter.read_track(path, settings)



def test_plotter_removes_partial_png_after_render_failure(tmp_path, monkeypatch):
    def fail_savefig(self, path, *args, **kwargs):
        Path(path).write_bytes(b"partial")
        raise OSError("render failed")

    monkeypatch.setattr(plotter.plt.Figure, "savefig", fail_savefig)
    inputs = plot_inputs(
        tmp_path,
        pd.DataFrame(columns=["event_start", "event_end", "summit", "is_change_type"]),
    )
    output = tmp_path / "plot.png"

    with pytest.raises(OSError, match="render failed"):
        plotter.plot_region(
            inputs["region"], inputs["events"], inputs["track_a"], inputs["track_b"],
            inputs["bins_a"], inputs["bins_b"], inputs["coolers_a"], inputs["coolers_b"],
            "H1398", "L134", inputs["settings"], output,
        )

    assert not output.exists()
    assert not list(tmp_path.glob(".plot.png.*.png"))

def test_plot_cli_requires_cooler_lists(monkeypatch):
    sys.path.insert(0, str(REPO / "scripts"))
    import plot_is_change_loci as plotter

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "plot_is_change_loci.py", "--primary", "primary.tsv",
            "--candidate-regions", "regions.tsv", "--is-a", "a.tsv", "--is-b", "b.tsv",
            "--config", "config.yaml", "--outdir", "plots", "--manifest", "manifest.tsv",
        ],
    )

    with pytest.raises(SystemExit):
        plotter.parse_args()


def test_plot_data_grays_invalid_cooler_weight_but_not_zero_is():
    sys.path.insert(0, str(REPO / "scripts"))
    import plot_is_change_loci as plotter
    from types import SimpleNamespace

    tracks = pd.DataFrame(
        {"chrom": ["chr1", "chr1"], "start": [0, 4000], "IS_distance": [32000, 32000], "IS": [0.0, 1.0]}
    )
    region = SimpleNamespace(chrom="chr1", region_start=0, region_end=4000)
    settings = SimpleNamespace(scales_bp=(32000,), resolution_bp=4000, context_bp=4000)
    bins_a = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 4000], "available": [True, False]})
    bins_b = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 4000], "available": [True, True]})

    frame = plotter.scale_plot_frame(region, tracks, tracks, bins_a, bins_b, 32000, settings)
