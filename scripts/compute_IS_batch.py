#!/usr/bin/env python3
"""
Compute insulation scores for a sample using cooltools insulation.

Runs cooltools insulation per chromosome (all distances in one call per chrom),
then saves the cooltools-native TSV for each chromosome to insulation/{sample}/.

cooltools insulation command (positional WINDOW args, multiple sizes in one call):
    cooltools insulation sample.cool 100000 200000 500000 1000000 -o out.tsv -p 8

Output per chromosome:
    results/insulation/{sample}/{sample}__{chrom}_insulation.tsv
    (cooltools-native format: chrom, start, end, log2_IS_{WINDOW}, ...)

Usage:
    python scripts/compute_IS_batch.py \
        --sample HC04_rep1 \
        --chroms "Ghir_A01 Ghir_A02 ..." \
        --outdir results \
        --is-distances 100000,200000,500000,1000000 \
        --tmpdir /tmp \
        --n-jobs 8 \
        --log-file results/logs/compute_IS/HC04_rep1.log
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def log(msg, log_file):
    print(msg)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{msg}\n")


def natural_sort_key(s):
    import re
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]


def run():
    parser = argparse.ArgumentParser(
        description="Run cooltools insulation per chromosome for a sample."
    )
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument(
        "--chroms", required=True,
        help="Space-separated chromosome list, e.g. 'Ghir_A01 Ghir_A02 ...'"
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Top-level results directory (OUTDIR)"
    )
    parser.add_argument(
        "--is-distances", required=True,
        help="Comma-separated IS distances in bp, e.g. '100000,200000,500000,1000000'"
    )
    parser.add_argument(
        "--tmpdir", default="/tmp",
        help="Temp directory for cooltools (default: /tmp)"
    )
    parser.add_argument(
        "--n-jobs", type=int, default=1,
        help="Number of parallel processes for cooltools (default: 1)"
    )
    parser.add_argument(
        "--log-file",
        help="Log file path"
    )
    args = parser.parse_args()

    if args.log_file:
        os.makedirs(os.path.dirname(args.log_file), exist_ok=True)
        Path(args.log_file).touch()
        log(f"[compute_IS_batch] sample={args.sample}", args.log_file)
        log(f"[compute_IS_batch] is_distances={args.is_distances}", args.log_file)

    chrom_list = args.chroms.split()
    distances_args = args.is_distances.split(",")

    insulation_dir = os.path.join(args.outdir, "insulation", args.sample)
    os.makedirs(insulation_dir, exist_ok=True)

    for chrom in chrom_list:
        base = f"{args.sample}__{chrom}"
        cool = os.path.join(args.outdir, "balanced", args.sample, f"{base}_iced.cool")
        tmp_out = os.path.join(args.tmpdir, f"{base}_insulation.tsv")
        out_tsv = os.path.join(insulation_dir, f"{base}_insulation.tsv")

        if not os.path.exists(cool):
            log(f"[WARN] Cool file not found: {cool} — skipping", args.log_file)
            continue

        log(f"[INFO] Running cooltools insulation for {base}", args.log_file)
        cmd = [
            "cooltools", "insulation",
            cool,
        ] + distances_args + [
            "-o", tmp_out,
            "-p", str(args.n_jobs),
        ]
        log(f"[INFO] CMD: {' '.join(cmd)}", args.log_file)

        try:
            subprocess.run(cmd, check=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            log(f"[ERROR] cooltools failed for {base}: {e}", args.log_file)
            raise

        os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
        if os.path.exists(tmp_out):
            os.rename(tmp_out, out_tsv)
            log(f"[INFO] Saved: {out_tsv}", args.log_file)
        else:
            log(f"[ERROR] Expected output not found: {tmp_out}", args.log_file)

    log(f"[compute_IS_batch] Done. {len(chrom_list)} chromosomes processed.", args.log_file)


if __name__ == "__main__":
    run()
