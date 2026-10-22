#!/usr/bin/env python3
"""
Convert cworld-dekker insulation bedGraph output to Meltron IS TSV format.

cworld bedGraph format:
    Line 1: track type=bedGraph name='...' description='...' maxHeightPixels=...
    Data lines: chrom\tstart\tend\tnormalized_insulation

Meltron IS format:
    chrom\tviewpoint\tIS_distance\tIS
    chrom\tstart\tdistance\tvalue

One output TSV per sample (all chromosomes and IS distances merged).

Usage:
    python scripts/convert_cworld_IS_to_meltron_IS.py \
        --sample HC04_rep1 \
        --chroms HC04_A01 HC04_A02 \
        --manifests results/cworld_insulation/HC04_rep1/HC04_rep1__HC04_A01.is_manifest.tsv results/cworld_insulation/HC04_rep1/HC04_rep1__HC04_A02.is_manifest.tsv \
        --outdir results \
        --is-distances 100000,200000,500000,1000000 \
        --expected-bin-size 4000 \
        --log /tmp/convert_IS.log
"""

import argparse
from collections import Counter
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


def positive_int(value):
    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be an integer") from error
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


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


def convert_manifest(manifest_path, distances, expected_bin_size, log_file):
    frames = []
    for distance in distances:
        bedgraph = manifest_bedgraph(manifest_path, distance)
        log(f"[INFO] Converting {bedgraph} (requested dist={distance})", log_file)
        frames.append(convert_bedgraph_to_meltron(bedgraph, distance, expected_bin_size, log_file))
    return frames


def convert_bedgraph_to_meltron(bedgraph_path, is_distance, expected_bin_size, log_file=None):
    """
    Convert a single cworld insulation bedGraph to Meltron long format.

    Parameters
    ----------
    bedgraph_path : str
        Path to cworld .insulation.bedGraph file
    is_distance : int
        Actual IS window size used by cworld (bp)
    expected_bin_size : int
        Configured source resolution in bp
    log_file : str, optional
        Log file path

    Returns
    -------
    pd.DataFrame
        Meltron format: chrom, bin-start viewpoint, IS_distance, IS
    """
    raw_rows = []
    with open(bedgraph_path, "r") as handle:
        for line in handle:
            if line.startswith("track type="):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start_text, end_text, value = parts[:4]
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as error:
                raise ValueError(f"{bedgraph_path}: non-integer interval coordinates") from error
            if end <= start:
                raise ValueError(f"{bedgraph_path}: non-positive bin width at {chrom}:{start}-{end}")
            raw_rows.append((chrom, start, end, value))

    if not raw_rows:
        return pd.DataFrame(columns=["chrom", "viewpoint", "IS_distance", "IS"])

    widths = Counter(end - start for _, start, end, _ in raw_rows)
    modal_count = max(widths.values())
    modal_widths = [width for width, count in widths.items() if count == modal_count]
    if len(modal_widths) != 1:
        raise ValueError(f"{bedgraph_path}: ambiguous modal bin widths {sorted(modal_widths)}")
    observed_bin_size = modal_widths[0]
    if observed_bin_size != expected_bin_size:
        raise ValueError(
            f"{bedgraph_path}: observed complete bin size {observed_bin_size} "
            f"does not match expected bin size {expected_bin_size}"
        )

    terminal_end = {}
    for chrom, _, end, _ in raw_rows:
        terminal_end[chrom] = max(terminal_end.get(chrom, end), end)

    output_rows = []
    for chrom, start, end, value in raw_rows:
        width = end - start
        if width != observed_bin_size:
            if 0 < width < observed_bin_size and end == terminal_end[chrom]:
                log(f"[INFO] Skipped terminal incomplete bin {chrom}:{start}-{end}", log_file)
                continue
            raise ValueError(
                f"{bedgraph_path}: non-terminal or oversized bin {chrom}:{start}-{end} "
                f"has width {width}; expected {observed_bin_size}"
            )
        if start % observed_bin_size != 0:
            raise ValueError(
                f"{bedgraph_path}: {chrom}:{start}-{end} is not aligned to bin size {observed_bin_size}"
            )
        try:
            is_numeric = pd.NA if value == "NA" else float(value)
        except ValueError:
            is_numeric = pd.NA
        output_rows.append({
            "chrom": chrom,
            "viewpoint": start,
            "IS_distance": is_distance,
            "IS": is_numeric,
        })

    return pd.DataFrame(output_rows, columns=["chrom", "viewpoint", "IS_distance", "IS"])


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
    parser.add_argument(
        "--expected-bin-size",
        required=True,
        type=positive_int,
        help="Configured source resolution in bp; must match complete bedGraph bins",
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
            manifest_frames = convert_manifest(
                manifest_path, distances, args.expected_bin_size, args.log
            )
            if any(frame.empty for frame in manifest_frames):
                raise ValueError(f"{manifest_path}: empty frame")
            frames.extend(manifest_frames)
    except ValueError as error:
        parser.error(str(error))

    if not frames:
        parser.error("No insulation bedGraph files selected from manifests")

    merged = pd.concat(frames, ignore_index=True)
    key_columns = ["chrom", "viewpoint", "IS_distance"]
    duplicate_rows = merged.duplicated(key_columns, keep=False)
    if duplicate_rows.any():
        duplicate_keys = merged.loc[duplicate_rows, key_columns].head(5).to_dict("records")
        parser.error(f"duplicate Meltron coordinates: {duplicate_keys}")

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
