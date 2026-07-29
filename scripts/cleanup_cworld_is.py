#!/usr/bin/env python3
"""Safely remove per-sample cworld insulation-score artifacts."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import sys

SAFE_SAMPLE_CHARS = frozenset(
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remove selected sample cworld insulation-score outputs."
    )
    parser.add_argument("--outdir", required=True, type=Path, help="Results directory")
    parser.add_argument(
        "--samples",
        required=True,
        nargs="+",
        help="Sample names to remove",
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Perform deletion; without this flag the command is a dry run",
    )
    return parser.parse_args()


def validate_sample(sample: str) -> None:
    if not sample or not set(sample) <= SAFE_SAMPLE_CHARS:
        raise ValueError(f"unsafe sample name: {sample!r}")


def collect_targets(root: Path, sample: str) -> list[Path]:
    return [
        root / "cworld_insulation" / sample,
        root / "IS" / sample,
        *sorted((root / "logs" / "compute_IS").glob(f"{sample}__*")),
        *sorted((root / "benchmark" / "compute_IS").glob(f"{sample}__*")),
        root / "logs" / "compute_IS" / f"{sample}_merge_IS.log",
        root / "benchmark" / "compute_IS" / f"{sample}_merge_IS.tsv",
    ]


def validate_target(root: Path, target: Path) -> None:
    resolved_target = target.resolve()
    if resolved_target == root or not resolved_target.is_relative_to(root):
        raise ValueError(f"unsafe cleanup target: {target}")


def remove_target(target: Path) -> None:
    if target.is_symlink() or target.is_file():
        target.unlink()
    elif target.is_dir():
        shutil.rmtree(target)


def main() -> int:
    args = parse_args()
    root = args.outdir.resolve()

    try:
        for sample in args.samples:
            validate_sample(sample)
        targets = [
            target
            for sample in args.samples
            for target in collect_targets(root, sample)
        ]
        for target in targets:
            validate_target(root, target)
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    for target in targets:
        if not target.exists() and not target.is_symlink():
            print(f"[SKIP] missing: {target}")
        elif args.execute:
            remove_target(target)
            print(f"[DELETE] removed: {target}")
        else:
            print(f"[DRY-RUN] would remove: {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
