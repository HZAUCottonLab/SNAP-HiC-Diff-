#!/usr/bin/env python3
"""Build merged IS-change candidate regions from KNN and Meltron TSVs."""

import argparse
import os
from pathlib import Path
import tempfile

import pandas as pd
import yaml

from cworld_is_config import parse_cworld_is_settings
from is_change_candidates import merge_candidate_windows, select_final_windows
from is_change_config import parse_is_change_settings


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--knn-enriched", nargs="+", required=True)
    parser.add_argument("--meltron-merged", nargs="+", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--out-candidates", required=True)
    parser.add_argument("--out-union-tsv", required=True)
    parser.add_argument("--out-union-bed", required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def log(message, path):
    print(message)
    log_path = Path(path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as handle:
        handle.write(f"{message}\n")


def atomic_table(frame, path, header=True):
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp"
    )
    os.close(descriptor)
    try:
        frame.to_csv(temporary, sep="\t", index=False, header=header)
        os.replace(temporary, output)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def main():
    args = parse_args()
    with open(args.config) as handle:
        config = yaml.safe_load(handle)
    settings = parse_is_change_settings(config, parse_cworld_is_settings(config))
    if len(args.knn_enriched) != len(args.meltron_merged):
        raise ValueError("KNN and Meltron input lists must have equal lengths")
    selected_scales = []
    for knn_path, meltron_path in zip(args.knn_enriched, args.meltron_merged):
        selected, audit = select_final_windows(
            pd.read_csv(knn_path, sep="\t"),
            pd.read_csv(meltron_path, sep="\t"),
            settings,
            return_audit=True,
        )
        selected_scales.append(selected)
        log(
            f"{knn_path}: matched={audit['matched_windows']} "
            f"selected={audit['selected_windows']}",
            args.log,
        )
    if all(selected.empty for selected in selected_scales):
        candidates = selected_scales[0].copy()
    else:
        candidates = pd.concat(selected_scales, ignore_index=True)
    candidates = candidates.sort_values(
        ["chrom", "window_start", "window_end", "window_bins"], kind="stable"
    ).reset_index(drop=True)
    regions = merge_candidate_windows(candidates, settings)
    bed = regions[["chrom", "region_start", "region_end", "region_id"]]
    atomic_table(candidates, args.out_candidates)
    atomic_table(regions, args.out_union_tsv)
    atomic_table(bed, args.out_union_bed, header=False)


if __name__ == "__main__":
    main()
