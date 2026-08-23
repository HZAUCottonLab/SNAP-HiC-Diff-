#!/usr/bin/env python3
"""Build an evidence ledger and genomic union from aligned KNN/Meltron scales."""

import argparse
import os
from pathlib import Path
import tempfile

import pandas as pd
import yaml

from cworld_is_config import parse_cworld_is_settings
from is_change_candidates import (
    SORT_COLUMNS,
    merge_candidate_windows,
    select_final_windows,
)
from is_change_config import parse_is_change_settings


BED_COLUMNS = ["chrom", "region_start", "region_end", "region_id"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--knn-enriched", required=True, nargs="+")
    parser.add_argument("--meltron-merged", required=True, nargs="+")
    parser.add_argument("--expected-window-bins", required=True, nargs="+", type=int)
    parser.add_argument("--config", required=True)
    parser.add_argument("--out-all", required=True)
    parser.add_argument("--out-regions", required=True)
    parser.add_argument("--out-bed", required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def log(message, path):
    print(message)
    log_path = Path(path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as handle:
        handle.write(f"{message}\n")


def _observed_scale(frame, label, path):
    required = ["resolution_bp", "window_bins"]
    missing = set(required).difference(frame.columns)
    if missing:
        raise ValueError(f"{label} {path} missing columns: {', '.join(sorted(missing))}")
    if frame.empty:
        return None

    values = frame.loc[:, required].apply(pd.to_numeric, errors="coerce")
    if values.isna().any().any():
        raise ValueError(f"{label} {path} has non-numeric scale provenance")
    unique = values.drop_duplicates()
    if len(unique) != 1:
        raise ValueError(f"{label} {path} must contain exactly one observed scale")
    resolution, window_bins = unique.iloc[0]
    if not float(resolution).is_integer() or not float(window_bins).is_integer():
        raise ValueError(f"{label} {path} has non-integer scale provenance")
    return int(resolution), int(window_bins)


def validate_aligned_scale_inputs(
    knn_frames, meltron_frames, knn_paths, meltron_paths, expected_window_bins, resolution_bp
):
    expected = tuple(expected_window_bins)
    if not expected or any(value <= 0 for value in expected):
        raise ValueError("--expected-window-bins must contain positive integers")
    if len(set(expected)) != len(expected):
        raise ValueError("--expected-window-bins contains duplicates")
    if not (
        len(knn_frames)
        == len(meltron_frames)
        == len(expected)
        == len(knn_paths)
        == len(meltron_paths)
    ):
        raise ValueError(
            "KNN inputs, Meltron inputs, and --expected-window-bins must have equal lengths"
        )

    observed = set()
    for expected_bins, knn, meltron, knn_path, meltron_path in zip(
        expected, knn_frames, meltron_frames, knn_paths, meltron_paths
    ):
        for frame, path, label in (
            (knn, knn_path, "KNN"),
            (meltron, meltron_path, "Meltron"),
        ):
            scale = _observed_scale(frame, label, path)
            if scale is None:
                continue
            if scale != (resolution_bp, expected_bins):
                raise ValueError(
                    f"{label} {path} reports r{scale[0]}/w{scale[1]}, "
                    f"expected r{resolution_bp}/w{expected_bins}"
                )
            observed.add(scale[1])

    unexpected = observed.difference(expected)
    if unexpected:
        raise ValueError(f"observed unexpected window bins: {sorted(unexpected)}")


def build_candidates(knn_frames, meltron_frames, settings):
    selected_frames = []
    audits = []
    for knn, meltron in zip(knn_frames, meltron_frames):
        selected, audit = select_final_windows(
            knn, meltron, settings, return_audit=True
        )
        selected_frames.append(selected)
        audits.append(audit)

    all_candidates = pd.concat(selected_frames, ignore_index=True, sort=False)
    if not all_candidates.empty:
        all_candidates = all_candidates.sort_values(
            SORT_COLUMNS, kind="stable"
        ).reset_index(drop=True)
    regions = merge_candidate_windows(all_candidates, settings)
    return all_candidates, regions, audits


def _temporary_sibling(path):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
    )
    os.close(descriptor)
    return Path(temporary)


def write_outputs_atomically(all_candidates, regions, out_all, out_regions, out_bed):
    destinations = [Path(out_all), Path(out_regions), Path(out_bed)]
    if len({destination.resolve() for destination in destinations}) != len(destinations):
        raise ValueError("all three output paths must be distinct")

    temporaries = [_temporary_sibling(destination) for destination in destinations]
    try:
        all_candidates.to_csv(temporaries[0], sep="\t", index=False)
        regions.to_csv(temporaries[1], sep="\t", index=False)
        regions.loc[:, BED_COLUMNS].to_csv(
            temporaries[2], sep="\t", index=False, header=False
        )
        for temporary, destination in zip(temporaries, destinations):
            os.replace(temporary, destination)
    except BaseException:
        for temporary in temporaries:
            temporary.unlink(missing_ok=True)
        raise


def main():
    args = parse_args()
    with open(args.config) as handle:
        config = yaml.safe_load(handle)
    settings = parse_is_change_settings(config, parse_cworld_is_settings(config))
    knn_frames = [pd.read_csv(path, sep="\t") for path in args.knn_enriched]
    meltron_frames = [pd.read_csv(path, sep="\t") for path in args.meltron_merged]
    validate_aligned_scale_inputs(
        knn_frames,
        meltron_frames,
        args.knn_enriched,
        args.meltron_merged,
        args.expected_window_bins,
        settings.resolution_bp,
    )
    all_candidates, regions, audits = build_candidates(
        knn_frames, meltron_frames, settings
    )
    write_outputs_atomically(
        all_candidates, regions, args.out_all, args.out_regions, args.out_bed
    )

    log(
        f"Matched {sum(audit['matched_windows'] for audit in audits)} "
        "KNN/Meltron windows by scale-qualified source identity",
        args.log,
    )
    log(
        f"Selected {len(all_candidates)} windows after KNN and Meltron thresholds",
        args.log,
    )
    log(
        f"Wrote {len(all_candidates)} source windows to {args.out_all}; "
        f"{len(regions)} union regions to {args.out_regions} and {args.out_bed}",
        args.log,
    )


if __name__ == "__main__":
    main()
