#!/usr/bin/env python3
"""Batch wrapper: compute O/E matrices for all chromosomes of one sample.

Usage:
    python compute_oe_matrix_batch.py \
        --chroms HC04_A01 HC04_A02 ... \
        --sample L128_rep1 \
        --outdir results \
        --inputs results/balanced/L128_rep1__HC04_A01_iced.cool ... \
        --expected results/expected/L128_rep1__HC04_A01_expected.npz ... \
        --log results/logs/oe/L128_rep1_oe.log
"""

import argparse
import os
import subprocess
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
COMPUTE_SCRIPT = os.path.join(SCRIPT_DIR, "compute_oe_matrix.py")


def main(chroms, sample, outdir, inputs, expected, log):
    log_base = log

    print(f"[INFO] Sample: {sample}")
    print(f"[INFO] Computing O/E matrices for {len(chroms)} chromosomes")

    os.makedirs(os.path.dirname(log_base), exist_ok=True)
    with open(log_base, "w") as lf:
        lf.write(f"[INFO] Sample: {sample}\n")
        lf.write(f"[INFO] Chromosomes: {chroms}\n\n")

    failures = []
    for i, chrom in enumerate(chroms):
        cool = inputs[i]
        npz  = expected[i]
        chrom_log = os.path.join(outdir, "logs", "oe",
                                f"{sample}__{chrom}_oe.log")
        output_oe = os.path.join(outdir, "oe", sample,
                                  f"{sample}__{chrom}_oe.npz")

        if not os.path.isfile(cool):
            msg = f"[ERROR] Input not found: {cool}"
            print(msg)
            with open(log_base, "a") as lf:
                lf.write(msg + "\n")
            failures.append(chrom)
            continue

        cmd = [
            sys.executable, COMPUTE_SCRIPT,
            "--input-balanced", cool,
            "--input-expected", npz,
            "--output", output_oe,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        with open(log_base, "a") as lf:
            lf.write(f"[INFO] Processing chromosome: {chrom}\n")
            lf.write(result.stdout)
            if result.stderr:
                lf.write(f"[STDERR] {result.stderr}\n")
            lf.write(f"[INFO] Finished chromosome: {chrom}\n\n")

        if result.returncode != 0:
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
        description="Batch compute O/E matrices for one sample."
    )
    parser.add_argument("--chroms", nargs="+", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--expected", nargs="+", required=True)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()
    sys.exit(main(**vars(args)))
