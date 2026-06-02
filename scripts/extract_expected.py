#!/usr/bin/env python3
"""Convert cooltools expected-cis tsv output to 1D expected npz array."""

import os
import sys
import traceback
import argparse

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", "-i", required=True,
                        help="Input cooltools expected-cis TSV file")
    parser.add_argument("--output", "-o", required=True,
                        help="Output NPZ file")
    parser.add_argument("--log", "-l", default="",
                        help="Log file path")
    parser.add_argument("--max-dist", "-m", type=int, default=10000,
                        help="Max distance in kb (default: 10000)")
    return parser.parse_args()


try:
    INPUT_TSV = snakemake.input.tsv
    OUTPUT = snakemake.output.expected
    LOG_FILE = snakemake.log[0]
    MAX_DIST = snakemake.params.get("max_dist", 10000)
except NameError:
    args = parse_args()
    INPUT_TSV = args.input
    OUTPUT = args.output
    LOG_FILE = args.log
    MAX_DIST = args.max_dist


def log_msg(log_file, msg):
    if log_file:
        with open(log_file, 'a') as f:
            f.write(msg + '\n')
    print(msg)


def main():
    if not INPUT_TSV or not OUTPUT:
        log_msg(LOG_FILE, "ERROR: Missing required inputs")
        return False

    log_msg(LOG_FILE, f"Converting cooltools tsv -> npz: {INPUT_TSV}")

    try:
        df = pd.read_csv(INPUT_TSV, sep='\t', usecols=['dist', 'balanced.avg', 'n_valid'])
        df = df[(df['dist'] <= MAX_DIST) & df['balanced.avg'].notna()].copy()
        log_msg(LOG_FILE, f"Loaded {len(df)} valid distance bins from cooltools output")

        n = min(int(df['dist'].max()) + 1, MAX_DIST)
        exp_arr = np.zeros(n, dtype=np.float64)
        cnt_arr = np.zeros(n, dtype=np.float64)

        for _, row in df.iterrows():
            d = int(row['dist'])
            if d < n:
                exp_arr[d] = row['balanced.avg']
                cnt_arr[d] = row['n_valid']

        n_nz = int(np.sum(cnt_arr > 0))
        sample_d = min(100, n - 1)
        log_msg(LOG_FILE, f"Expected array: n_dist={n}, non_zero={n_nz}, "
                          f"sample[d={sample_d}]={exp_arr[sample_d]:.4f}")

        os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
        np.savez_compressed(OUTPUT, expected=exp_arr, counts=cnt_arr)
        log_msg(LOG_FILE, f"Saved to {OUTPUT}")
        log_msg(LOG_FILE, "SUCCESS")
        return True

    except Exception as e:
        log_msg(LOG_FILE, f"ERROR: {e}\n{traceback.format_exc()}")
        return False


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
