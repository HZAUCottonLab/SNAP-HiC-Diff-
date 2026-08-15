from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "plot_is_change_loci.py"
SCALES = (32000, 60000, 92000, 120000, 152000, 180000, 212000, 240000, 272000, 300000)


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


def write_tracks(root):
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
            "--config", str(config), "--outdir", str(tmp_path / "plots"),
            "--manifest", str(tmp_path / "plots" / "manifest.tsv"),
        ], text=True, capture_output=True, check=False,
    )

    assert result.returncode == 0, result.stderr
    paths = pd.read_csv(tmp_path / "plots" / "manifest.tsv", sep="\t")["plot_path"]
    assert all(Path(path).is_file() for path in paths)
