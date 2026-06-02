#!/usr/bin/env python3
"""Batch wrapper: convert all cooltools tsv -> npz for one sample.

Receives data via command-line arguments (argparse) to avoid pickle round-trip
issues with Snakemake's dynamically created NamedList subclasses (InputFiles /
OutputFiles), which cannot survive pickle/unpickle across different Python
environments on cluster nodes.

Usage:
    python extract_expected_batch.py \\
        --chroms HC04_A01 HC04_A02 ... \\
        --outdir results \\
        --sample L128_rep2 \\
        --max-dist 10000 \\
        --inputs results/expected/L128_rep2__HC04_A01_cooltools.tsv ... \\
        --outputs results/expected/L128_rep2__HC04_A01_expected.npz ... \\
        --log results/logs/expected/L128_rep2_to_npz.log
"""

import argparse
import os
import sys
import subprocess


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
EXTRACT_SCRIPT = os.path.join(SCRIPT_DIR, "extract_expected.py")


def run_extract_single(tsv_path, npz_path, log_file, max_dist):
    cmd = [
        sys.executable, EXTRACT_SCRIPT,
        "--input", tsv_path,
        "--output", npz_path,
        "--log", log_file,
        "--max-dist", str(max_dist),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0, result.stdout, result.stderr


def main(chroms, outdir, sample, max_dist, inputs, outputs, log):
    log_base = log

    print(f"[INFO] Sample: {sample}")
    print(f"[INFO] Converting TSV -> NPZ for {len(chroms)} chromosomes")

    os.makedirs(os.path.dirname(log_base), exist_ok=True)
    with open(log_base, "w") as lf:
        lf.write(f"[INFO] Sample: {sample}\n")
        lf.write(f"[INFO] Chromosomes: {chroms}\n")
        lf.write(f"[INFO] Max distance: {max_dist} kb\n\n")

    failures = []
    for i, chrom in enumerate(chroms):
        tsv      = inputs[i]
        npz      = outputs[i]
        chrom_log = os.path.join(outdir, "logs", "expected",
                                 f"{sample}__{chrom}_to_npz.log")

        if not os.path.isfile(tsv):
            msg = f"[ERROR] Input not found: {tsv}"
            print(msg)
            with open(log_base, "a") as lf:
                lf.write(msg + "\n")
            failures.append(chrom)
            continue

        ok, stdout, stderr = run_extract_single(tsv, npz, chrom_log, max_dist)
        with open(log_base, "a") as lf:
            lf.write(f"[INFO] Processing chromosome: {chrom}\n")
            lf.write(stdout)
            if stderr:
                lf.write(f"[STDERR] {stderr}\n")
            lf.write(f"[INFO] Finished chromosome: {chrom}\n\n")

        if not ok:
            failures.append(chrom)
            print(f"[ERROR] Failed chromosome: {chrom}")
        else:
            print(f"[INFO] Finished chromosome: {chrom}")

    print()
    if failures:
        print(f"[ERROR] Failed chromosomes: {failures}")
        sys.exit(1)
    else:
        print(f"[INFO] Finished sample: {sample} — all chromosomes done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Batch convert cooltools tsv -> npz for one sample."
    )
    parser.add_argument("--chroms", nargs="+", required=True,
                        help="List of chromosome names")
    parser.add_argument("--outdir", required=True,
                        help="Base output directory")
    parser.add_argument("--sample", required=True,
                        help="Sample name (wildcard value)")
    parser.add_argument("--max-dist", type=int, required=True,
                        help="Maximum distance in bins")
    parser.add_argument("--inputs", nargs="+", required=True,
                        help="List of input TSV paths (same order as --chroms)")
    parser.add_argument("--outputs", nargs="+", required=True,
                        help="List of output NPZ paths (same order as --chroms)")
    parser.add_argument("--log", required=True,
                        help="Path to consolidated log file")
    args = parser.parse_args()
    sys.exit(main(**vars(args)))
