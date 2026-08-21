from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
import sys

import cooler

from tad_errors import TadError


@dataclass(eq=True)
class Interval:
    chrom: str
    start: int
    end: int
    resolutions: set
    n_source_boundaries: int


def load_chrom_sizes(mcool, resolution):
    clr = cooler.Cooler(f"{Path(mcool).resolve()}::/resolutions/{resolution}")
    return {str(chrom): int(length) for chrom, length in clr.chromsizes.items()}


def read_tads(path, *, sample, signal, resolution):
    rows = []
    for line_number, raw_line in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        columns = line.split()
        try:
            rows.append((columns[0], int(columns[1]), int(columns[2])))
        except (IndexError, ValueError) as error:
            raise TadError(
                "TAD_INVALID_COORDINATE",
                sample=sample,
                signal=signal,
                resolution=resolution,
                path=path,
                line_number=line_number,
                line=line,
            ) from error
    return rows


def boundaries_from_tads(tads, chrom_sizes, resolution, *, sample, signal):
    boundaries = []
    for chrom, start, end in tads:
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
                chrom=chrom,
                start=start,
                end=end,
            )
        chrom_length = chrom_sizes[chrom]
        for position in (start, end):
            boundaries.append(
                Interval(
                    chrom,
                    max(0, position - resolution),
                    min(chrom_length, position + resolution),
                    {resolution},
                    1,
                )
            )
    return boundaries


def merge_intervals(intervals, chrom_order):
    ordered = sorted(
        intervals,
        key=lambda item: (chrom_order[item.chrom], item.start, item.end),
    )
    merged = []
    for item in ordered:
        if (
            merged
            and merged[-1].chrom == item.chrom
            and item.start < merged[-1].end
        ):
            merged[-1].end = max(merged[-1].end, item.end)
            merged[-1].resolutions.update(item.resolutions)
            merged[-1].n_source_boundaries += item.n_source_boundaries
        else:
            merged.append(
                Interval(
                    item.chrom,
                    item.start,
                    item.end,
                    set(item.resolutions),
                    item.n_source_boundaries,
                )
            )
    return merged


def atomic_write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text)
    temporary.replace(path)


def write_resolution_boundaries(path, intervals, resolution):
    lines = []
    for index, item in enumerate(intervals, 1):
        lines.append(
            f"{item.chrom}\t{item.start}\t{item.end}\t"
            f"boundary_{index:06d}\t{resolution}\t"
            f"{item.n_source_boundaries}\n"
        )
    atomic_write(path, "".join(lines))


def read_resolution_boundaries(path, *, sample, signal):
    intervals = []
    for line_number, line in enumerate(Path(path).read_text().splitlines(), 1):
        columns = line.split("\t")
        try:
            chrom, start, end = columns[0], int(columns[1]), int(columns[2])
            resolution = int(columns[4])
            count = int(columns[5])
        except (IndexError, ValueError) as error:
            raise TadError(
                "TAD_BOUNDARY_MERGE_FAILED",
                sample=sample,
                signal=signal,
                path=path,
                line_number=line_number,
                line=line,
            ) from error
        intervals.append(Interval(chrom, start, end, {resolution}, count))
    return intervals


def write_multires_boundaries(path, intervals):
    lines = []
    for index, item in enumerate(intervals, 1):
        resolutions = sorted(item.resolutions)
        lines.append(
            f"{item.chrom}\t{item.start}\t{item.end}\t"
            f"boundary_{index:06d}\t"
            f"{','.join(map(str, resolutions))}\t{len(resolutions)}\t"
            f"{item.n_source_boundaries}\n"
        )
    atomic_write(path, "".join(lines))


def build_parser():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    extract = subparsers.add_parser("extract")
    extract.add_argument("--tad", required=True)
    extract.add_argument("--mcool", required=True)
    extract.add_argument("--resolution", required=True, type=int)
    extract.add_argument("--sample", required=True)
    extract.add_argument("--signal", required=True, choices=("raw", "balanced"))
    extract.add_argument("--output", required=True)

    merge = subparsers.add_parser("merge")
    merge.add_argument("--inputs", required=True, nargs="+")
    merge.add_argument("--mcool", required=True)
    merge.add_argument("--chrom-resolution", required=True, type=int)
    merge.add_argument("--sample", required=True)
    merge.add_argument("--signal", required=True, choices=("raw", "balanced"))
    merge.add_argument("--output", required=True)
    return parser


def run_from_args(args):
    chrom_resolution = (
        args.resolution if args.command == "extract" else args.chrom_resolution
    )
    chrom_sizes = load_chrom_sizes(args.mcool, chrom_resolution)
    chrom_order = {chrom: index for index, chrom in enumerate(chrom_sizes)}
    if args.command == "extract":
        intervals = boundaries_from_tads(
            read_tads(
                args.tad,
                sample=args.sample,
                signal=args.signal,
                resolution=args.resolution,
            ),
            chrom_sizes,
            args.resolution,
            sample=args.sample,
            signal=args.signal,
        )
        write_resolution_boundaries(
            args.output,
            merge_intervals(intervals, chrom_order),
            args.resolution,
        )
        return

    intervals = []
    for input_path in args.inputs:
        intervals.extend(
            read_resolution_boundaries(
                input_path,
                sample=args.sample,
                signal=args.signal,
            )
        )
    for item in intervals:
        if (
            item.chrom not in chrom_sizes
            or item.start < 0
            or item.start >= item.end
            or item.end > chrom_sizes.get(item.chrom, -1)
        ):
            raise TadError(
                "TAD_BOUNDARY_MERGE_FAILED",
                sample=args.sample,
                signal=args.signal,
                path=args.output,
                chrom=item.chrom,
                start=item.start,
                end=item.end,
            )
    write_multires_boundaries(
        args.output,
        merge_intervals(intervals, chrom_order),
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
            "TAD_BOUNDARY_MERGE_FAILED",
            sample=args.sample,
            signal=args.signal,
            path=args.output,
            error=type(error).__name__,
        )
        print(wrapped, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
