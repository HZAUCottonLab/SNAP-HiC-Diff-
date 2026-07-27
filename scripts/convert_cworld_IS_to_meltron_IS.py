#!/usr/bin/env python3
"""
Convert cworld-dekker insulation bedGraph output to Meltron IS TSV format.

cworld bedGraph format:
    Line 1: track type=bedGraph name='...' description='...' maxHeightPixels=...
    Data lines: chrom\tstart\tend\tnormalized_insulation

Meltron IS format:
    chrom\tviewpoint\tIS_distance\tIS
    chrom\t(start+end)//2\tdistance\tvalue

One output TSV per sample (all chromosomes and IS distances merged).

Usage:
    python scripts/convert_cworld_IS_to_meltron_IS.py \
        --sample HC04_rep1 \
        --chroms "HC04_A01 HC04_A02" \
        --outdir results \
        --is-distances 100000,200000,500000,1000000 \
        --log /tmp/convert_IS.log
"""

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import pandas as pd


def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def natural_sort_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]


def extract_is_from_filename(filename):
    """Extract --is value from cworld bedGraph filename like:
    {prefix}--is100000--nt0--ids50000--ss0--immean.insulation.bedGraph
    """
    m = re.search(r'--is(\d+)', filename)
    if m:
        return int(m.group(1))
    return None


def convert_bedgraph_to_meltron(bedgraph_path, is_distance, log_file=None):
    """
    Convert a single cworld insulation bedGraph to Meltron long format.

    Parameters
    ----------
    bedgraph_path : str
        Path to cworld .insulation.bedGraph file
    is_distance : int
        Actual IS window size used by cworld (bp)
    log_file : str, optional
        Log file path

    Returns
    -------
    pd.DataFrame
        Meltron format: chrom, viewpoint, IS_distance, IS
    """
    rows = []
    with open(bedgraph_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("track type="):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, is_val = parts[0], parts[1], parts[2], parts[3]
            try:
                viewpoint = (int(start) + int(end)) // 2
            except ValueError:
                continue
            if is_val == "NA":
                is_numeric = pd.NA
            else:
                try:
                    is_numeric = float(is_val)
                except ValueError:
                    is_numeric = pd.NA
            rows.append({
                "chrom": chrom,
                "viewpoint": viewpoint,
                "IS_distance": is_distance,
                "IS": is_numeric
            })

    return pd.DataFrame(rows)


def run():
    parser = argparse.ArgumentParser(
        description="Convert cworld insulation bedGraph files to Meltron IS format."
    )
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument(
        "--chroms", required=True,
        help="Space-separated chromosome list"
    )
    parser.add_argument("--outdir", required=True, help="Top-level results directory (OUTDIR)")
    parser.add_argument(
        "--is-distances", required=True,
        help="Comma-separated IS distances in bp"
    )
    parser.add_argument("--log", help="Log file path")
    args = parser.parse_args()

    if args.log:
        Path(args.log).parent.mkdir(parents=True, exist_ok=True)
        Path(args.log).touch()

    distances = [int(d.strip()) for d in args.is_distances.split(",")]
    chrom_list = args.chroms.split()

    log(f"[convert_cworld_IS_to_meltron_IS] sample={args.sample}", args.log)
    log(f"[convert_cworld_IS_to_meltron_IS] distances={args.is_distances}", args.log)

    cworld_dir = os.path.join(args.outdir, "cworld_insulation", args.sample)
    is_out = os.path.join(args.outdir, "IS", args.sample, f"{args.sample}_IS.tsv")

    all_rows = []
    missing = []
    for chrom in chrom_list:
        for dist in distances:
            manifest_path = Path(cworld_dir, f"{args.sample}__{chrom}.is_manifest.tsv")
            if manifest_path.exists():
                manifest_matches = []
                for line in manifest_path.read_text().splitlines():
                    requested_dist, separator, filename = line.partition("\t")
                    if separator and requested_dist == str(dist):
                        manifest_matches.append(filename)
                if len(manifest_matches) != 1:
                    missing.append(
                        f"{chrom} dist={dist}: expected one entry in {manifest_path}, "
                        f"found {len(manifest_matches)}"
                    )
                    continue
                bedgraph_path = Path(cworld_dir, manifest_matches[0])
                if not bedgraph_path.is_file():
                    missing.append(f"{chrom} dist={dist}: manifest target missing: {bedgraph_path}")
                    continue
            else:
                legacy_pattern = os.path.join(
                    cworld_dir,
                    f"{args.sample}__{chrom}--is{dist}--*.insulation.bedGraph",
                )
                matched = glob.glob(legacy_pattern)
                if len(matched) != 1:
                    all_pattern = os.path.join(
                        cworld_dir,
                        f"{args.sample}__{chrom}--is*--*.insulation.bedGraph",
                    )
                    matched = glob.glob(all_pattern)
                if len(matched) != 1:
                    missing.append(
                        f"{chrom} dist={dist}: expected one unmapped CWorld bedGraph, "
                        f"found {len(matched)}"
                    )
                    continue
                bedgraph_path = Path(matched[0])

            log(f"[INFO] Converting {bedgraph_path} (requested dist={dist})", args.log)
            df = convert_bedgraph_to_meltron(bedgraph_path, dist, args.log)
            all_rows.append(df)

    if missing:
        for item in missing:
            log(f"[ERROR] Missing CWorld bedGraph: {item}", args.log)
        raise SystemExit(1)

    if not all_rows:
        log("[ERROR] No insulation bedGraph files found", args.log)
        sys.exit(1)

    merged = pd.concat(all_rows, ignore_index=True)
    merged = merged.sort_values(
        by="chrom",
        key=lambda col: col.map(natural_sort_key)
    ).reset_index(drop=True)

    merged = merged[["chrom", "viewpoint", "IS_distance", "IS"]]

    os.makedirs(os.path.dirname(is_out), exist_ok=True)
    temp_out = f"{is_out}.tmp"
    merged.to_csv(temp_out, sep="\t", index=False)
    os.replace(temp_out, is_out)

    log(f"[INFO] Meltron IS file written: {is_out} ({len(merged)} rows)", args.log)


if __name__ == "__main__":
    run()
