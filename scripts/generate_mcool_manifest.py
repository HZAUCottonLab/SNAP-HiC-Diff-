#!/usr/bin/env python3
"""Generate manifest_mcool.tsv from detect_manifest.tsv.

Reads detect_manifest.tsv, converts paths based on action:
  - SKIP   -> use original cool_uri
  - CONVERT -> use outdir/{sample}.mcool
"""

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Generate mcool manifest")
    parser.add_argument("--detect-manifest", required=True, help="detect_manifest.tsv")
    parser.add_argument("--outdir", required=True, help="mcool output directory")
    parser.add_argument("--output", required=True, help="Output manifest file")
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.output, "w") as out:
        out.write("sample_id\tcool_uri\tassembly\n")

        with open(args.detect_manifest) as f:
            header = f.readline()  # skip header
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue

                sample_id, cool_uri, assembly, action = parts[0], parts[1], parts[2], parts[3]

                if action == "SKIP":
                    out.write(f"{sample_id}\t{cool_uri}\t{assembly}\n")
                elif action == "CONVERT":
                    mcool_path = f"{args.outdir}/{sample_id}.mcool"
                    out.write(f"{sample_id}\t{mcool_path}\t{assembly}\n")
                else:
                    print(f"[WARN] Skipping sample: {sample_id} (action={action})", file=sys.stderr)

    with open(args.output) as f:
        lines = sum(1 for _ in f)

    print(f"[INFO] manifest_mcool.tsv: {lines} lines written")


if __name__ == "__main__":
    main()
