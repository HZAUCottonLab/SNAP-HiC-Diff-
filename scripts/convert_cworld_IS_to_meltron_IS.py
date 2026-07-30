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
        --chroms HC04_A01 HC04_A02 \
        --manifests results/cworld_insulation/HC04_rep1/HC04_rep1__HC04_A01.is_manifest.tsv results/cworld_insulation/HC04_rep1/HC04_rep1__HC04_A02.is_manifest.tsv \
        --outdir results \
        --is-distances 100000,200000,500000,1000000 \
        --log /tmp/convert_IS.log
"""

import argparse
import os
import re
from pathlib import Path

import pandas as pd


def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def natural_sort_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]


def manifest_bedgraph(manifest_path, distance):
    if not manifest_path.is_file():
        raise ValueError(f"{manifest_path}: manifest missing")

    matches = []
    for line in manifest_path.read_text().splitlines():
        requested, separator, filename = line.partition("\t")
        if separator and requested == str(distance):
            matches.append(filename)

    if len(matches) != 1:
        raise ValueError(
            f"{manifest_path}: distance {distance} requires exactly one entry; "
            f"found {len(matches)}"
        )

    root = manifest_path.parent.resolve()
    bedgraph = (root / matches[0]).resolve()
    if root not in bedgraph.parents or not bedgraph.is_file():
        raise ValueError(f"{manifest_path}: manifest target missing: {bedgraph}")
    return bedgraph


def convert_manifest(manifest_path, distances, log_file):
    frames = []
    for distance in distances:
        bedgraph = manifest_bedgraph(manifest_path, distance)
        log(f"[INFO] Converting {bedgraph} (requested dist={distance})", log_file)
        frames.append(convert_bedgraph_to_meltron(bedgraph, distance, log_file))
    return frames


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
    parser.add_argument("--chroms", nargs="+", required=True, help="Chromosome list")
    parser.add_argument("--manifests", nargs="+", required=True, help="Manifest paths")
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
    chrom_list = args.chroms
    manifest_paths = [Path(path) for path in args.manifests]
    resolved_manifest_paths = [path.resolve() for path in manifest_paths]
    if len(set(resolved_manifest_paths)) != len(resolved_manifest_paths):
        parser.error("duplicate manifest path")

    if len(chrom_list) != len(manifest_paths):
        parser.error("--chroms and --manifests must have equal lengths")

    log(f"[convert_cworld_IS_to_meltron_IS] sample={args.sample}", args.log)
    log(f"[convert_cworld_IS_to_meltron_IS] distances={args.is_distances}", args.log)

    frames = []
    try:
        for chrom, manifest_path in zip(chrom_list, manifest_paths, strict=True):
            log(f"[INFO] Converting chromosome {chrom} from {manifest_path}", args.log)
            manifest_frames = convert_manifest(manifest_path, distances, args.log)
            if any(frame.empty for frame in manifest_frames):
                raise ValueError(f"{manifest_path}: empty frame")
            frames.extend(manifest_frames)
    except ValueError as error:
        parser.error(str(error))

    if not frames:
        parser.error("No insulation bedGraph files selected from manifests")

    merged = pd.concat(frames, ignore_index=True)
    merged = merged.sort_values(
        by="chrom",
        key=lambda col: col.map(natural_sort_key)
    ).reset_index(drop=True)

    merged = merged[["chrom", "viewpoint", "IS_distance", "IS"]]

    is_out = Path(args.outdir, "IS", args.sample, f"{args.sample}_IS.tsv")
    is_out.parent.mkdir(parents=True, exist_ok=True)
    temp_out = is_out.with_name(f".{is_out.name}.tmp")
    merged.to_csv(temp_out, sep="\t", index=False)
    os.replace(temp_out, is_out)

    log(f"[INFO] Meltron IS file written: {is_out} ({len(merged)} rows)", args.log)


if __name__ == "__main__":
    run()
