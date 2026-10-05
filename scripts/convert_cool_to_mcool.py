#!/usr/bin/env python3
"""Convert .cool files to .mcool via cooler zoomify.

Reads detect_manifest.tsv, performs action per sample:
  - SKIP   -> touch done marker, do nothing
  - CONVERT -> cooler zoomify -> results/mcool/{sample}.mcool
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Convert cool to mcool")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--detect-manifest", required=True, help="detect_manifest.tsv")
    parser.add_argument("--outdir", required=True, help="Output directory for mcool files")
    parser.add_argument("--resolutions", required=True, help="Comma-separated resolutions")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads")
    parser.add_argument("--log-file", required=True, help="Log file path")
    parser.add_argument("--touch-file", required=True, help="Touch file to create on success")
    return parser.parse_args()


def read_action(detect_manifest: str, sample: str) -> tuple:
    """Read action for sample from detect_manifest.tsv.

    Returns:
        Tuple of (cool_uri, assembly, action)
    """
    with open(detect_manifest) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] == sample:
                return parts[1], parts[2], parts[3]
    raise ValueError(f"Sample {sample} not found in {detect_manifest}")


def run_convert(cool_uri: str, mcool_out: str, resolutions: list, threads: int, log_file):
    """Run cooler zoomify."""
    cmd = [
        "cooler", "zoomify",
        cool_uri,
        "-o", mcool_out,
        "-r", ",".join(map(str, resolutions)),
        "--balance",
        "--nproc", str(threads),
    ]
    print(f"[INFO] Converting: {cool_uri} -> {mcool_out}", file=log_file)
    log_file.flush()
    subprocess.run(cmd, check=True)
    print(f"[INFO] Done: {mcool_out}", file=log_file)


def main():
    args = parse_args()
    sample = args.sample
    resolutions = [int(r) for r in args.resolutions.split(",")]

    with open(args.log_file, "w") as log_file:
        print(f"[INFO] Sample: {sample}", file=log_file)

        cool_uri, assembly, action = read_action(args.detect_manifest, sample)
        print(f"[INFO] Action: {action}", file=log_file)

        if action == "SKIP":
            print("[INFO] Already .mcool, skipping", file=log_file)
            Path(args.touch_file).touch()
            return

        if action == "UNKNOWN":
            print(f"[ERROR] Unknown format for sample: {sample}", file=sys.stderr)
            sys.exit(1)

        # ACTION == CONVERT
        os.makedirs(args.outdir, exist_ok=True)
        mcool_out = os.path.join(args.outdir, f"{sample}.mcool")

        run_convert(cool_uri, mcool_out, resolutions, args.threads, log_file)
        Path(args.touch_file).touch()


if __name__ == "__main__":
    main()
