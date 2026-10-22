#!/usr/bin/env python3
"""
Python wrapper for Meltron R CLI.

Meltron.R calculates melting scores for supplied genomic intervals by comparing
insulation scores between two samples.

Usage:
    python scripts/run_meltron.py \
        --comparison Mutant1_vs_Mutant2 \
        --is-a results/IS/Mutant1/Mutant1_IS.tsv \
        --is-b results/IS/Mutant2/Mutant2_IS.tsv \
        --genomic-intervals results/meltron/Mutant1_vs_Mutant2/Mutant1_vs_Mutant2_intervals.tsv \
        --outfile results/meltron/Mutant1_vs_Mutant2/Mutant1_vs_Mutant2_decondensation.tsv \
        --comparison-direction decondensation \
        --cutoff-na 0.5 \
        --log-file results/logs/meltron/Mutant1_vs_Mutant2.log
"""

import argparse
import subprocess
import sys


def run():
    parser = argparse.ArgumentParser(
        description="Run Meltron R CLI for insulation score comparison."
    )
    parser.add_argument(
        "--is-a", required=True,
        help="Path to reference sample IS table (4-col TSV: chrom, viewpoint, IS_distance, IS)"
    )
    parser.add_argument(
        "--is-b", required=True,
        help="Path to input sample IS table (4-col TSV: chrom, viewpoint, IS_distance, IS)"
    )
    parser.add_argument(
        "--genomic-intervals", required=True,
        help="Path to genomic intervals (4-col TSV: chrom, start, end, ID)"
    )
    parser.add_argument(
        "--outfile", required=True,
        help="Path to output TSV"
    )
    parser.add_argument(
        "--comparison-direction", required=True,
        choices=["decondensation", "condensation"],
        help="'decondensation' for melting events, 'condensation' for chromatin compaction"
    )
    parser.add_argument(
        "--cutoff-na", type=float, default=0.5,
        help="Max fraction of missing IS values per interval (0-0.99, default: 0.5)"
    )
    parser.add_argument(
        "--cores", type=int, default=None,
        help="Number of CPU cores for data.table (default: all available)"
    )
    parser.add_argument(
        "--log-file",
        help="Log file path"
    )
    args = parser.parse_args()

    meltron_r = "env/Meltron/command_line_apps/Meltron.R"

    cmd = [
        "Rscript",
        meltron_r,
    ]
    if args.cores is not None:
        cmd.extend(["--cores", str(args.cores)])
    cmd.extend([
        "--referenceIS", args.is_a,
        "--inputIS", args.is_b,
        "--genomicIntervals", args.genomic_intervals,
        "--outfile", args.outfile,
        "--comparison", args.comparison_direction,
        "--cutoff_NA", str(args.cutoff_na),
    ])

    print(f"[run_meltron] CMD: {' '.join(cmd)}")
    if args.log_file:
        with open(args.log_file, "a") as f:
            f.write(f"[run_meltron] Running {args.comparison_direction}\n")
            f.write(f"[run_meltron] CMD: {' '.join(cmd)}\n")

    result = subprocess.run(cmd, check=True, capture_output=False, text=True)
    print(f"[run_meltron] Done: {args.outfile}")


if __name__ == "__main__":
    run()
