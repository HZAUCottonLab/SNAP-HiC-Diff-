#!/usr/bin/env python3
"""
Split QC-pass windows into CPU-balanced chunks per sample pair.
"""

import os
import sys
import pandas as pd
from pathlib import Path

try:
    INPUT_QC = snakemake.input.qc
    OUTPUT_DIR = str(Path(snakemake.output[0]).parent)
    N_CHUNKS = snakemake.params.n_chunks
    PAIR_ID = snakemake.wildcards.pair_id
    LOG_FILE = snakemake.log[0]
except NameError:
    INPUT_QC = os.environ.get("INPUT_QC", "")
    OUTPUT_DIR = os.environ.get("OUTPUT_DIR", "")
    N_CHUNKS = int(os.environ.get("N_CHUNKS", "10"))
    PAIR_ID = os.environ.get("PAIR_ID", "")
    LOG_FILE = os.environ.get("LOG_FILE", "")


def log(msg):
    print(msg)


def split_into_chunks(qc_df, n_chunks):
    """Split QC-pass windows evenly across chunks."""
    # Filter PASS windows only
    qc_pass = qc_df[qc_df["status"] == "PASS"].copy()
    n_windows = len(qc_pass)

    log(f"Splitting {n_windows} PASS windows into {n_chunks} chunks")

    if n_windows == 0:
        log("WARNING: No PASS windows found!")
        return []

    # Adjust n_chunks if needed
    n_chunks = min(n_chunks, n_windows)

    # Calculate chunk sizes
    base_size = n_windows // n_chunks
    remainder = n_windows % n_chunks

    chunks = []
    idx = 0
    for i in range(n_chunks):
        chunk_size = base_size + (1 if i < remainder else 0)
        if chunk_size > 0:
            chunk_windows = qc_pass.iloc[idx:idx + chunk_size]
            chunks.append(chunk_windows)
            idx += chunk_size

    return chunks


def main():
    input_qc = os.environ.get("INPUT_QC", "")
    output_dir = os.environ.get("OUTPUT_DIR", "")
    n_chunks = int(os.environ.get("N_CHUNKS", "10"))
    pair_id = os.environ.get("PAIR_ID", "")
    log_file = os.environ.get("LOG_FILE", "")

    log(f"Splitting windows into chunks for pair {pair_id}")
    log(f"  Input: {input_qc}")
    log(f"  Output dir: {output_dir}")
    log(f"  N chunks: {n_chunks}")

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Load QC results
    if not input_qc or not os.path.exists(input_qc):
        log(f"ERROR: QC file not found: {input_qc}")
        return False

    qc_df = pd.read_parquet(input_qc)
    log(f"Loaded {len(qc_df)} windows from QC")

    # Split into chunks
    chunks = split_into_chunks(qc_df, n_chunks)

    log(f"Created {len(chunks)} chunks")

    # Save each chunk
    for i, chunk_df in enumerate(chunks):
        chunk_id = f"{i:03d}"
        output_path = Path(output_dir) / f"chunk_{chunk_id}_windows.parquet"
        chunk_df.to_parquet(output_path, index=False)
        log(f"  Chunk {chunk_id}: {len(chunk_df)} windows -> {output_path}")

    # Save chunk manifest
    manifest_path = Path(output_dir) / "chunk_manifest.tsv"
    manifest_data = []
    for i, chunk_df in enumerate(chunks):
        chunk_id = f"{i:03d}"
        manifest_data.append({
            "chunk_id": chunk_id,
            "n_windows": len(chunk_df),
            "path": f"chunk_{chunk_id}_windows.parquet"
        })
    manifest_df = pd.DataFrame(manifest_data)
    manifest_df.to_csv(manifest_path, sep="\t", index=False)
    log(f"Saved chunk manifest: {manifest_path}")

    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
