#!/usr/bin/env python3
"""
Validate cworld-dekker IS computation pipeline.

Single-sample × single-chromosome × single-window validation:
1. cooler → cworld matrix format correct
2. cworld insulation runs successfully
3. bedGraph output format correct
4. Meltron TSV conversion correct
5. Compare with existing cooltools IS (if available)

Usage:
    python scripts/validate_cworld_is.py \
        --cool results/balanced/H1394/H1394__HC04_A01_iced.cool \
        --sample H1394 \
        --chrom HC04_A01 \
        --is-distance 100000 \
        --cworld-dir env/cworld-dekker \
        --outdir /tmp/cworld_validate
"""

import argparse
import glob
import gzip
import os
import re
import subprocess
import sys
from pathlib import Path

import cooler
import numpy as np
import pandas as pd


def log(msg):
    print(msg)
    sys.stdout.flush()


def validate_matrix_format(matrix_path, chrom):
    """Validate cworld matrix format matches expected structure."""
    log(f"[CHECK] Matrix format: {matrix_path}")
    with gzip.open(matrix_path, "rt") as f:
        header = f.readline()
        parts = header.rstrip("\n").split("\t")
        if parts[0] != "":
            log(f"[FAIL] Header first cell should be empty, got: {parts[0]}")
            return False
        if not parts[1].startswith("HiC_1|aaa|"):
            log(f"[FAIL] Header first bin label wrong: {parts[1]}")
            return False
        if chrom not in parts[1]:
            log(f"[FAIL] Chromosome {chrom} not in header: {parts[1]}")
            return False
        log(f"[PASS] Header format correct. {len(parts)-1} bins.")
    return True


def validate_bedgraph(bedgraph_path, n_bins_expected):
    """Validate cworld bedGraph format."""
    log(f"[CHECK] bedGraph format: {bedgraph_path}")
    data_lines = 0
    has_track_header = False
    with open(bedgraph_path, "r") as f:
        for line in f:
            if line.startswith("track type="):
                has_track_header = True
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                data_lines += 1
    if not has_track_header:
        log(f"[FAIL] Missing track header line")
        return False
    if data_lines == 0:
        log(f"[FAIL] No data lines found")
        return False
    log(f"[PASS] bedGraph: {data_lines} data lines, track header present.")
    return True


def validate_meltron_tsv(tsv_path):
    """Validate Meltron IS TSV format."""
    log(f"[CHECK] Meltron TSV: {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t")
    required_cols = {"chrom", "viewpoint", "IS_distance", "IS"}
    if not required_cols.issubset(set(df.columns)):
        missing = required_cols - set(df.columns)
        log(f"[FAIL] Missing columns: {missing}")
        return False
    if len(df) == 0:
        log(f"[FAIL] Empty TSV")
        return False
    na_count = df["IS"].isna().sum()
    valid_count = (~df["IS"].isna()).sum()
    log(f"[PASS] TSV: {len(df)} rows, {valid_count} valid IS, {na_count} NA")
    return True


def compare_with_cooltools(cooltools_tsv, meltron_tsv):
    """Compare cworld-derived IS with existing cooltools IS."""
    log(f"[CHECK] Comparing with cooltools IS")
    ct = pd.read_csv(cooltools_tsv, sep="\t")
    mt = pd.read_csv(meltron_tsv, sep="\t")

    ct_view = ct[["chrom", "start", "end", "log2_insulation_score_100000"]].copy()
    ct_view["viewpoint"] = (ct_view["start"] + ct_view["end"]) // 2
    ct_view = ct_view.rename(columns={"log2_insulation_score_100000": "IS_ct"})
    ct_view = ct_view[ct_view["chrom"].str.startswith("HC04_A01")]

    mt_100k = mt[mt["IS_distance"] == 100000].copy()

    merged = pd.merge(ct_view, mt_100k, on=["chrom", "viewpoint"], how="inner")
    valid = merged.dropna(subset=["IS_ct", "IS"])

    if len(valid) < 10:
        log(f"[WARN] Too few overlapping points to compare: {len(valid)}")
        return None

    corr_pearson = valid["IS_ct"].corr(valid["IS"])
    corr_spearman = valid["IS_ct"].corr(valid["IS"], method="spearman")

    log(f"[COMPARE] Pearson r={corr_pearson:.4f} (n={len(valid)})")
    log(f"[COMPARE] Spearman rho={corr_spearman:.4f}")
    log(f"[COMPARE] mean_ct={valid['IS_ct'].mean():.4f}, mean_cworld={valid['IS'].mean():.4f}")

    return {"pearson": corr_pearson, "spearman": corr_spearman, "n": len(valid)}


def run():
    parser = argparse.ArgumentParser(description="Validate cworld IS pipeline")
    parser.add_argument("--cool", required=True, help="Input .cool file")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--chrom", required=True, help="Chromosome name")
    parser.add_argument("--is-distance", type=int, required=True, help="IS window size (bp)")
    parser.add_argument("--cworld-dir", required=True, help="Path to cworld-dekker directory")
    parser.add_argument("--outdir", required=True, help="Output directory for validation")
    parser.add_argument(
        "--cooltools-is",
        help="Existing cooltools IS TSV for comparison (optional)"
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    cworld_sample_dir = os.path.join(args.outdir, "cworld_insulation", args.sample)
    os.makedirs(cworld_sample_dir, exist_ok=True)

    matrix_out = os.path.join(cworld_sample_dir, f"{args.sample}__{args.chrom}.matrix.gz")
    prefix_out = os.path.join(cworld_sample_dir, f"{args.sample}__{args.chrom}")
    bedgraph_out = os.path.join(args.outdir, f"{args.sample}__{args.chrom}.insulation.bedGraph")

    log(f"=== cworld IS Validation ===")
    log(f"Sample: {args.sample}, Chrom: {args.chrom}, IS: {args.is_distance}")
    log(f"Cool: {args.cool}")
    log(f"Cworld dir: {args.cworld_dir}")
    log(f"Output dir: {args.outdir}")

    log(f"\n--- Step 1: Convert cool to cworld matrix ---")
    cmd_convert = [
        sys.executable,
        "scripts/convert_cool_to_cworld.py",
        "--cool", args.cool,
        "--chrom", args.chrom,
        "--out", matrix_out,
        "--log", os.path.join(args.outdir, "convert_cool.log"),
    ]
    log(f"CMD: {' '.join(cmd_convert)}")
    result = subprocess.run(cmd_convert, check=True, capture_output=True, text=True)
    log(result.stdout)
    if result.returncode != 0:
        log(f"[FAIL] convert_cool_to_cworld failed: {result.stderr}")
        sys.exit(1)
    log(f"[PASS] Matrix conversion done: {matrix_out}")

    if not validate_matrix_format(matrix_out, args.chrom):
        log(f"[FAIL] Matrix format validation failed")
        sys.exit(1)

    log(f"\n--- Step 2: Run cworld insulation ---")
    perl_script = os.path.join(args.cworld_dir, "scripts/perl/matrix2insulation.pl")
    if not os.path.exists(perl_script):
        log(f"[FAIL] Perl script not found: {perl_script}")
        sys.exit(1)

    cmd_cworld = [
        "perl", "-I", os.path.join(args.cworld_dir, "lib"),
        perl_script,
        "-v",
        "--is", str(args.is_distance),
        "--ids", str(args.is_distance // 2),
        "-i", matrix_out,
        "-o", prefix_out,
    ]
    log(f"CMD: {' '.join(cmd_cworld)}")
    env = os.environ.copy()
    env["PERL5LIB"] = f"{args.cworld_dir}/lib:{env.get('PERL5LIB', '')}"
    result = subprocess.run(cmd_cworld, env=env, check=False, capture_output=True, text=True)
    Path(args.outdir, "cworld_stdout.log").write_text(result.stdout)
    Path(args.outdir, "cworld_stderr.log").write_text(result.stderr)
    log(result.stdout)
    if result.stderr:
        log(f"[CWorld STDERR]\n{result.stderr}")
    if result.returncode != 0:
        log(f"[FAIL] cworld insulation failed: {result.stderr}")
        sys.exit(1)

    prefix_artifacts = sorted(glob.glob(prefix_out + "*"))
    Path(args.outdir, "cworld_artifacts.txt").write_text(
        "\n".join(prefix_artifacts) + ("\n" if prefix_artifacts else "")
    )
    log(f"[CHECK] CWorld prefix artifacts: {len(prefix_artifacts)}")
    for artifact in prefix_artifacts:
        log(f"[ARTIFACT] {artifact}")

    matched_bedgraph = glob.glob(prefix_out + "--is" + str(args.is_distance) + "--*--immean.insulation.bedGraph")
    if not matched_bedgraph:
        log(f"[FAIL] No bedGraph output found for prefix: {prefix_out}")
        sys.exit(1)
    bedgraph_path = matched_bedgraph[0]
    log(f"[PASS] cworld insulation done: {bedgraph_path}")

    with cooler.Cooler(args.cool) as clr:
        n_bins = len(clr.bins()[:][clr.bins()[:]["chrom"] == args.chrom])
    if not validate_bedgraph(bedgraph_path, n_bins):
        log(f"[FAIL] bedGraph validation failed")
        sys.exit(1)

    log(f"\n--- Step 3: Convert to Meltron TSV ---")
    cworld_dir_for_convert = args.outdir
    cmd_convert_is = [
        sys.executable,
        "scripts/convert_cworld_IS_to_meltron_IS.py",
        "--sample", args.sample,
        "--chroms", args.chrom,
        "--outdir", cworld_dir_for_convert,
        "--is-distances", str(args.is_distance),
        "--log", os.path.join(args.outdir, "convert_IS.log"),
    ]
    log(f"CMD: {' '.join(cmd_convert_is)}")
    result = subprocess.run(cmd_convert_is, check=True, capture_output=True, text=True)
    log(result.stdout)
    if result.returncode != 0:
        log(f"[FAIL] convert_cworld_IS_to_meltron_IS failed: {result.stderr}")
        sys.exit(1)

    meltron_tsv = os.path.join(cworld_dir_for_convert, "IS", args.sample, f"{args.sample}_IS.tsv")
    if not os.path.exists(meltron_tsv):
        log(f"[FAIL] Meltron TSV not found: {meltron_tsv}")
        sys.exit(1)
    if not validate_meltron_tsv(meltron_tsv):
        log(f"[FAIL] Meltron TSV validation failed")
        sys.exit(1)

    log(f"\n--- Step 4: Compare with cooltools IS (if available) ---")
    if args.cooltools_is and os.path.exists(args.cooltools_is):
        result = compare_with_cooltools(args.cooltools_is, meltron_tsv)
        if result and result["pearson"] < 0.9:
            log(f"[WARN] Pearson correlation {result['pearson']:.4f} < 0.9 — investigate difference")
        elif result:
            log(f"[PASS] Good correlation with cooltools IS")
    else:
        log(f"[SKIP] No cooltools IS provided for comparison")

    log(f"\n=== Validation Complete ===")


if __name__ == "__main__":
    run()
