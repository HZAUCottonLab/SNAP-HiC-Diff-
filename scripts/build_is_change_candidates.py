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
    parser.add_argument("--knn-enriched", required=True)
    parser.add_argument("--meltron-merged", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def log(message, path):
    print(message)
    log_path = Path(path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as handle:
        handle.write(f"{message}\n")


def main():
    args = parse_args()
    with open(args.config) as handle:
        config = yaml.safe_load(handle)
    settings = parse_is_change_settings(config, parse_cworld_is_settings(config))
    selected = select_final_windows(
        pd.read_csv(args.knn_enriched, sep="\t"),
        pd.read_csv(args.meltron_merged, sep="\t"),
        settings,
    )
    regions = merge_candidate_windows(selected, settings)

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp"
    )
    os.close(descriptor)
    try:
        regions.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, output)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise
    log(f"Wrote {len(regions)} candidate regions to {output}", args.log)


if __name__ == "__main__":
    main()
