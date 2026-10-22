#!/usr/bin/env python3
"""
Detect cool sample formats and determine conversion actions.
Reads input manifest, determines if each sample needs conversion (.cool -> .mcool).
Outputs a TSV with action column: SKIP or CONVERT.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Detect cool sample formats.")
    p.add_argument("--input", required=True, help="Input manifest TSV (sample_id, cool_uri, assembly)")
    p.add_argument("--output", required=True, help="Output detect manifest TSV")
    p.add_argument("--log", required=True, help="Log file path")
    return p.parse_args()


def detect_action(cool_uri):
    """Determine conversion action based on file extension."""
    if cool_uri.endswith(".mcool"):
        return "SKIP"
    elif cool_uri.endswith(".cool"):
        return "CONVERT"
    else:
        return "UNKNOWN"


def main():
    args = parse_args()

    with open(args.log, "w") as log_f:
        log_f.write("")
        print(f"[detect_cool_samples] Input: {args.input}")

        df = pd.read_csv(args.input, sep="\t")
        df["action"] = df["cool_uri"].apply(detect_action)

        unknown = df[df["action"] == "UNKNOWN"]
        if not unknown.empty:
            for _, row in unknown.iterrows():
                msg = f"[WARN] Unknown extension: {row['sample_id']} ({row['cool_uri']})"
                print(msg, file=sys.stderr)
                log_f.write(msg + "\n")

        df_out = df[["sample_id", "cool_uri", "assembly", "action"]]
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df_out.to_csv(args.output, sep="\t", index=False)

        n_samples = len(df)
        n_skip = (df["action"] == "SKIP").sum()
        n_convert = (df["action"] == "CONVERT").sum()
        print(f"[detect_cool_samples] Done: {n_samples} samples ({n_skip} skip, {n_convert} convert)")


if __name__ == "__main__":
    main()
