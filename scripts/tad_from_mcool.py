from argparse import ArgumentParser
import json
from pathlib import Path
import sys

import cooler
import h5py
import numpy as np

from tad_errors import TadError


VALID_SIGNALS = {"raw", "balanced"}


def validate_mcool_resolution(mcool_path, resolution, sample, signal):
    bins_path = f"resolutions/{resolution}/bins"
    weight_path = f"{bins_path}/weight"
    with h5py.File(mcool_path, "r") as handle:
        if bins_path not in handle:
            raise TadError(
                "TAD_MISSING_RESOLUTION",
                sample=sample,
                signal=signal,
                resolution=resolution,
                path=mcool_path,
            )
        if signal == "balanced" and weight_path not in handle:
            raise TadError(
                "TAD_MISSING_WEIGHT",
                sample=sample,
                signal=signal,
                resolution=resolution,
                path=mcool_path,
            )


def iter_pixels_as_triples(clr, signal):
    if signal not in VALID_SIGNALS:
        raise TadError("TAD_CONFIG_INVALID", signal=signal)
    pixels = clr.pixels()[:]
    weights = None
    if signal == "balanced":
        weights = clr.bins()[:]["weight"].to_numpy(dtype=float)
    for bin1, bin2, count in zip(
        pixels["bin1_id"].to_numpy(),
        pixels["bin2_id"].to_numpy(),
        pixels["count"].to_numpy(dtype=float),
    ):
        if not np.isfinite(count) or count <= 0:
            continue
        value = count
        if weights is not None:
            value = count * weights[int(bin1)] * weights[int(bin2)]
        if np.isfinite(value) and value > 0:
            yield int(bin1) + 1, int(bin2) + 1, float(value)


def write_bed_from_bins(path: Path, clr: cooler.Cooler) -> int:
    bins = clr.bins()[:]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for bin_id, (chrom, start, end) in enumerate(
            zip(bins["chrom"], bins["start"], bins["end"]), start=1
        ):
            handle.write(f"{chrom}\t{start}\t{end}\t{bin_id}\n")
    return len(bins)


def write_interactions(path, triples, *, sample, signal, resolution):
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w") as handle:
        for bin1, bin2, value in triples:
            handle.write(f"{bin1}\t{bin2}\t{value}\n")
            count += 1
    if count == 0:
        raise TadError(
            "TAD_EMPTY_INTERACTIONS",
            sample=sample,
            signal=signal,
            resolution=resolution,
            path=path,
        )
    return count


def effective_d_config(d_config, resolution: int):
    effective = dict(d_config)
    resolution_overrides = effective.pop("resolution_overrides", {})
    effective.update(
        resolution_overrides.get(
            str(resolution), resolution_overrides.get(resolution, {})
        )
    )
    return effective


def render_settings(
    path: Path,
    *,
    sample: str,
    resolution: int,
    d_config,
    bad_region_bed,
    badregionfilter: bool,
    excluded_chromosomes,
    processors: int,
    bed_file=None,
    counts_file=None,
) -> None:
    effective = effective_d_config(d_config, resolution)
    sizes = effective["size_bins"]
    variances = effective["variance_thresholds"]
    bed_file = bed_file or f"{sample}.bed"
    counts_file = counts_file or f"{sample}.interactions"

    settings = [
        ("bed_file", bed_file),
        ("sample_1", sample),
        ("counts_file_1", counts_file),
        ("region_size", effective["region_size_bins"]),
        ("overlap", effective["overlap_bins"]),
        ("logged", effective["logged"]),
        ("scale", effective["scale"]),
        ("plateau", effective["plateau"]),
        ("chaosfilter", effective["chaosfilter"]),
        ("chaos_pct", effective["chaos_pct"]),
        ("max_gamma", effective["max_gamma"]),
        ("num_part", effective["num_part"]),
        ("plots", effective["plots"]),
        ("pctile_threshold", effective["pctile_threshold"]),
        ("pct_value", effective["pct_value"]),
        ("final_consistent_domains", effective["final_consistent_domains"]),
        ("size_threshold", effective["size_threshold_bins"]),
    ]
    settings.extend(
        (f"size_s{index}", size * resolution)
        for index, size in enumerate(sizes, start=1)
    )
    settings.extend(
        (f"var_thresh{index}", variance)
        for index, variance in enumerate(variances, start=1)
    )
    settings.extend(
        [
            ("boundary_buffer", effective["boundary_buffer_bins"] * resolution),
            ("badregionfile", bad_region_bed),
            ("badregionfilter", badregionfilter),
            ("excluded_chromosomes", ",".join(excluded_chromosomes)),
            ("processors", processors),
        ]
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for key, value in settings:
            handle.write(f"{key}\t{value}\n")


def parse_boolean(value: str) -> bool:
    normalized = value.lower()
    if normalized in {"true", "1"}:
        return True
    if normalized in {"false", "0"}:
        return False
    raise ValueError(f"expected a boolean value, got: {value}")


def build_parser() -> ArgumentParser:
    parser = ArgumentParser()
    parser.add_argument("--mcool", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--signal", required=True, choices=sorted(VALID_SIGNALS))
    parser.add_argument("--resolution", required=True, type=int)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--ready-flag", required=True)
    parser.add_argument("--bad-region-bed", required=True)
    parser.add_argument("--badregionfilter", required=True)
    parser.add_argument("--excluded-chromosomes", default="")
    parser.add_argument("--d-config-json", required=True)
    parser.add_argument("--processors", required=True, type=int)
    return parser


def run_from_args(args) -> None:
    mcool_path = Path(args.mcool).resolve()
    validate_mcool_resolution(mcool_path, args.resolution, args.sample, args.signal)

    run_dir = Path(args.run_dir).resolve()
    ready_flag = Path(args.ready_flag).resolve()
    input_dir = run_dir / "input"
    bed_path = input_dir / f"{args.sample}.bed"
    interactions_path = input_dir / f"{args.sample}.{args.signal}.interactions"
    settings_path = run_dir / f"settings_{args.sample}.{args.signal}.txt"
    clr = cooler.Cooler(f"{mcool_path}::/resolutions/{args.resolution}")

    n_bins = write_bed_from_bins(bed_path, clr)
    n_interactions = write_interactions(
        interactions_path,
        iter_pixels_as_triples(clr, args.signal),
        sample=args.sample,
        signal=args.signal,
        resolution=args.resolution,
    )
    render_settings(
        settings_path,
        sample=args.sample,
        resolution=args.resolution,
        d_config=json.loads(args.d_config_json),
        bad_region_bed=args.bad_region_bed,
        badregionfilter=parse_boolean(args.badregionfilter),
        excluded_chromosomes=(
            args.excluded_chromosomes.split(",") if args.excluded_chromosomes else []
        ),
        processors=args.processors,
        bed_file=bed_path.name,
        counts_file=interactions_path.name,
    )

    ready_flag.parent.mkdir(parents=True, exist_ok=True)
    ready_text = "\n".join(
        [
            f"sample={args.sample}",
            f"signal={args.signal}",
            f"mcool={mcool_path}",
            f"resolution={args.resolution}",
            f"run_dir={run_dir}",
            f"settings={settings_path}",
            f"n_bins={n_bins}",
            f"n_interactions={n_interactions}",
        ]
    ) + "\n"
    ready_tmp = ready_flag.with_name(ready_flag.name + ".tmp")
    ready_tmp.write_text(ready_text)
    ready_tmp.replace(ready_flag)


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
            stage="prepare",
            returncode="NA",
            path=Path(args.run_dir).resolve(),
            error=type(error).__name__,
        )
        print(wrapped, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
