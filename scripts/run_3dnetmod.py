#!/usr/bin/env python3
from argparse import ArgumentParser
import os
from pathlib import Path
import subprocess
import sys

import cooler

from tad_errors import TadError


STAGES = [
    "3DNetMod_preprocess.py",
    "3DNetMod_GPS_MMCP.py",
    "3DNetMod_HSVM.py",
]


def build_env(code_dir, run_dir):
    code_dir = Path(code_dir).resolve()
    run_dir = Path(run_dir).resolve()
    env = os.environ.copy()
    env["PY3_OUTDIR"] = str(run_dir)
    env["PY3_RUNDIR"] = str(run_dir)
    env["PY3_CODEDIR"] = str(code_dir)
    env["PYTHONPATH"] = ":".join(
        [str(code_dir), env.get("PYTHONPATH", "")]
    ).rstrip(":")
    return env


def run_stages(settings, code_dir, run_dir, *, sample, signal, resolution):
    settings = Path(settings).resolve()
    code_dir = Path(code_dir).resolve()
    run_dir = Path(run_dir).resolve()
    env = build_env(code_dir, run_dir)
    for stage in STAGES:
        try:
            subprocess.run(
                [
                    str(Path(sys.executable).resolve()),
                    str((code_dir / stage).resolve()),
                    str(settings),
                ],
                cwd=run_dir,
                env=env,
                check=True,
            )
        except subprocess.CalledProcessError as error:
            raise TadError(
                "TAD_STAGE_FAILED",
                sample=sample,
                signal=signal,
                resolution=resolution,
                stage=stage,
                returncode=error.returncode,
                path=run_dir,
            ) from error


def read_chrom_sizes(mcool, resolution):
    clr = cooler.Cooler(f"{Path(mcool).resolve()}::/resolutions/{resolution}")
    return {str(chrom): int(length) for chrom, length in clr.chromsizes.items()}


def atomic_write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text)
    temporary.replace(path)


def collect_final_calls(
    run_dir, bed_output, bedpe_output, chrom_sizes,
    *, sample, signal, resolution,
):
    run_dir = Path(run_dir).resolve()
    final_dir = run_dir / "FINAL_DOMAIN_CALLS"
    domain_files = sorted(final_dir.glob("*_master_consistent_domains.txt"))
    chrom_order = {chrom: index for index, chrom in enumerate(chrom_sizes)}
    calls = []
    for domain_file in domain_files:
        for raw_line in domain_file.read_text().splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            columns = line.split()
            try:
                chrom = columns[0]
                start = int(columns[1])
                end = int(columns[2])
            except (IndexError, ValueError) as error:
                raise TadError(
                    "TAD_INVALID_COORDINATE",
                    sample=sample,
                    signal=signal,
                    resolution=resolution,
                    path=domain_file,
                    line=line,
                ) from error
            if (
                chrom not in chrom_sizes
                or start < 0
                or start >= end
                or end > chrom_sizes.get(chrom, -1)
            ):
                raise TadError(
                    "TAD_INVALID_COORDINATE",
                    sample=sample,
                    signal=signal,
                    resolution=resolution,
                    path=domain_file,
                    line=line,
                )
            calls.append((chrom, start, end))
    calls.sort(key=lambda row: (chrom_order[row[0]], row[1], row[2]))
    bed_text = "".join(f"{chrom}\t{start}\t{end}\n" for chrom, start, end in calls)
    bedpe_text = "".join(
        f"{chrom}\t{start}\t{end}\t{chrom}\t{start}\t{end}\n"
        for chrom, start, end in calls
    )
    atomic_write(bed_output, bed_text)
    atomic_write(bedpe_output, bedpe_text)


def build_parser():
    parser = ArgumentParser()
    parser.add_argument("--settings", required=True)
    parser.add_argument("--code-dir", required=True)
    parser.add_argument("--mcool", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--signal", required=True, choices=("raw", "balanced"))
    parser.add_argument("--resolution", required=True, type=int)
    parser.add_argument("--bed-output", required=True)
    parser.add_argument("--bedpe-output", required=True)
    return parser


def run_from_args(args):
    settings = Path(args.settings).resolve()
    run_dir = settings.parent
    run_stages(
        settings,
        args.code_dir,
        run_dir,
        sample=args.sample,
        signal=args.signal,
        resolution=args.resolution,
    )
    collect_final_calls(
        run_dir,
        args.bed_output,
        args.bedpe_output,
        read_chrom_sizes(args.mcool, args.resolution),
        sample=args.sample,
        signal=args.signal,
        resolution=args.resolution,
    )


def cli(argv=None):
    args = build_parser().parse_args(argv)
    try:
        run_from_args(args)
    except TadError as error:
        print(error, file=sys.stderr)
        return 1
    except Exception as error:
        wrapped = TadError(
            "TAD_STAGE_FAILED",
            sample=args.sample,
            signal=args.signal,
            resolution=args.resolution,
            stage="runner",
            returncode="NA",
            path=Path(args.settings).resolve().parent,
            error=type(error).__name__,
        )
        print(wrapped, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
