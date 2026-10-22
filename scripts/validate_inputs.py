#!/usr/bin/env python3
"""
Validate mcool files and extract chromosome list.
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Validate mcool files.")
    p.add_argument("--input", required=True, help="Input manifest TSV")
    p.add_argument("--output", required=True, help="Output validated TSV")
    p.add_argument("--log", required=True, help="Log file path")
    return p.parse_args()


def log(msg, fh):
    fh.write(f"{msg}\n")
    fh.flush()
    print(msg)


def run_hictk_validate(mcool_path):
    cmd = ["hictk", "validate", mcool_path]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        return result.returncode == 0, result.stdout + result.stderr
    except Exception as e:
        return False, str(e)


def get_uri(cool_path, resolution):
    """Build hictk URI. For .cool files use direct path (no resolution suffix).
    For .mcool files append ::/resolutions/{resolution}."""
    if cool_path.endswith(".cool"):
        return cool_path
    return f"{cool_path}::/resolutions/{resolution}"


def run_hictk_dump_chroms(mcool_path, resolution):
    uri = get_uri(mcool_path, resolution)
    cmd = ["hictk", "dump", uri, "--table", "chroms"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode == 0:
            lines = result.stdout.strip().split("\n")
            chroms = []
            for line in lines:
                parts = line.split()
                if len(parts) >= 2:
                    chroms.append(parts[0])
            return True, chroms
        # Failed: return actual output for diagnostics
        return False, result.stdout.strip().split("\n") if result.stdout else []
    except Exception as e:
        return False, []


def run_hictk_dump_resolutions(mcool_path):
    cmd = ["hictk", "dump", mcool_path, "--table", "resolutions"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode == 0:
            lines = result.stdout.strip().split("\n")
            resolutions = []
            for line in lines:
                parts = line.split()
                if len(parts) >= 1:
                    resolutions.append(int(parts[0]))
            return True, resolutions
        return False, []
    except Exception:
        return False, []


def main():
    args = parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.log).parent.mkdir(parents=True, exist_ok=True)

    log(f"Starting validation: {args.input}", sys.stdout)

    manifest = pd.read_csv(args.input, sep="\t")
    log(f"Found {len(manifest)} samples in manifest", sys.stdout)

    import yaml
    config_path = os.path.join(os.path.dirname(__file__), "..", "config", "config.yaml")
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    target_res = cfg.get("resolution", 32000)
    log(f"Target resolution: {target_res}", sys.stdout)

    validated_rows = []

    for idx, row in manifest.iterrows():
        sample_id = row["sample_id"]
        cool_uri = row["cool_uri"]

        log(f"\nValidating sample: {sample_id}", sys.stdout)
        log(f"  URI: {cool_uri}", sys.stdout)

        if not os.path.exists(cool_uri):
            log(f"  ERROR: File not found", sys.stdout)
            continue

        valid, msg = run_hictk_validate(cool_uri)
        if not valid:
            log(f"  ERROR: hictk validate failed - {msg}", sys.stdout)
            continue
        log(f"  OK: hictk validate passed", sys.stdout)

        has_resolution, resolutions = run_hictk_dump_resolutions(cool_uri)
        if not has_resolution:
            log(f"  ERROR: Could not extract resolutions", sys.stdout)
            continue
        log(f"  Resolutions: {resolutions}", sys.stdout)

        has_chroms, chroms = run_hictk_dump_chroms(cool_uri, target_res)
        if not has_chroms:
            log(f"  ERROR: Could not extract chromosomes — actual: {chroms[:20]}", sys.stdout)
            continue

        log(f"  Chromosomes ({len(chroms)}): {chroms[:5]}...", sys.stdout)

        validated_rows.append({
            "sample_id": sample_id,
            "cool_uri": cool_uri,
            "assembly": row.get("assembly", "unknown"),
            "chrom_list": ",".join(chroms),
            "resolution_list": ",".join(map(str, resolutions)),
            "status": "PASS"
        })

    validated_df = pd.DataFrame(validated_rows)
    validated_df.to_csv(args.output, sep="\t", index=False)

    log(f"\nValidation complete: {len(validated_rows)}/{len(manifest)} samples passed", sys.stdout)
    log(f"Output: {args.output}", sys.stdout)

    return len(validated_rows) > 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
