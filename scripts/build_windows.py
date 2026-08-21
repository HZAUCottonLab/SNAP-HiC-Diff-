#!/usr/bin/env python3
"""
Build genomic windows for analysis.
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Build genomic windows.")
    p.add_argument("--input", required=True, help="Validated manifest TSV")
    p.add_argument("--output", required=True, help="Output parquet path")
    p.add_argument("--log", required=True, help="Log file path")
    p.add_argument("--resolution", type=int, default=32000)
    p.add_argument("--window-bins", type=int, required=True)
    p.add_argument("--window-size", type=int, default=300000)
    p.add_argument("--step-bins", type=int, required=True)
    p.add_argument("--step-size", type=int, default=300000)
    return p.parse_args()


def log(msg, fh):
    fh.write(f"{msg}\n")
    fh.flush()


def get_chrom_sizes_from_mcool(cool_uri, resolution):
    """
    Get chromosome sizes from mcool file.
    """
    import subprocess

    uri = f"{cool_uri}::/resolutions/{resolution}"
    cmd = ["hictk", "dump", uri, "--table", "chroms"]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            return {}

        chrom_sizes = {}
        lines = result.stdout.strip().split("\n")
        for line in lines:
            parts = line.split()
            if len(parts) >= 2:
                chrom = parts[0]
                size = int(parts[1])  # hictk outputs: chrom<sep>size (2 columns)
                chrom_sizes[chrom] = size

        return chrom_sizes

    except Exception as e:
        log(f"Error getting chrom sizes: {e}")
        return {}


def build_windows(chrom_sizes, window_size, step_size, resolution, window_bins, step_bins):
    """
    Generate sliding windows across genome.

    Args:
        chrom_sizes: dict of {chrom: length}
        window_size: window size in bp (default 300000)
        step_size: sliding step in bp (default 32000)
        resolution: bin resolution in bp (default 4000)
        window_bins: window width in bins
        step_bins: sliding step in bins

    Returns:
        DataFrame with window information
    """
    window_data = []

    for chrom, chrom_len in sorted(chrom_sizes.items()):
        for start in range(0, chrom_len, step_size):
            end = start + window_size
            if end > chrom_len:
                if start >= chrom_len - window_size * 0.5:
                    end = chrom_len
                else:
                    break

            if end - start < window_size * 0.5:
                continue

            window_data.append({
                "window_id": f"{chrom}_{start}_{end}",
                "chrom": chrom,
                "start": start,
                "end": end,
                "n_bins": (end - start) // resolution,
                "window_size": end - start,
                "resolution_bp": resolution,
                "window_bins": window_bins,
                "window_size_bp": window_size,
                "step_bins": step_bins,
                "step_size_bp": step_size,
            })

    return pd.DataFrame(window_data)


def main():
    args = parse_args()
    resolution = args.resolution
    window_bins = args.window_bins
    window_size = args.window_size
    step_bins = args.step_bins
    step_size = args.step_size

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.log).parent.mkdir(parents=True, exist_ok=True)

    with open(args.log, "w") as log_fh:
        log(f"Starting window building", log_fh)
        log(f"  Input: {args.input}", log_fh)
        log(f"  Output: {args.output}", log_fh)
        log(f"  Resolution: {resolution}", log_fh)
        log(f"  Window bins: {window_bins}", log_fh)
        log(f"  Window size: {window_size}", log_fh)
        log(f"  Step bins: {step_bins}", log_fh)
        log(f"  Step size: {step_size}", log_fh)

        manifest = pd.read_csv(args.input, sep="\t")
        log(f"Found {len(manifest)} validated samples", log_fh)

        all_chrom_sizes = {}
        for _, row in manifest.iterrows():
            sample_id = row["sample_id"]
            cool_uri = row["cool_uri"]

            chrom_sizes = get_chrom_sizes_from_mcool(cool_uri, resolution)
            log(f"Sample {sample_id}: {len(chrom_sizes)} chromosomes", log_fh)

            for chrom, size in chrom_sizes.items():
                if chrom not in all_chrom_sizes or size > all_chrom_sizes[chrom]:
                    all_chrom_sizes[chrom] = size

        log(f"Total unique chromosomes: {len(all_chrom_sizes)}", log_fh)
        log(f"Chromosomes: {list(all_chrom_sizes.keys())}", log_fh)

        windows_df = build_windows(
            all_chrom_sizes,
            window_size,
            step_size,
            resolution,
            window_bins,
            step_bins,
        )
        log(f"Generated {len(windows_df)} windows", log_fh)

        chrom_counts = windows_df["chrom"].value_counts().to_dict()
        log(f"Windows per chromosome:", log_fh)
        for chrom, count in sorted(chrom_counts.items()):
            log(f"  {chrom}: {count}", log_fh)

        windows_df.to_parquet(args.output, index=False)
        log(f"Saved to {args.output}", log_fh)

    return len(windows_df) > 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
