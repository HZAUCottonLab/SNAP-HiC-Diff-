from argparse import ArgumentError, ArgumentParser
import json
from pathlib import Path
import sys
from uuid import uuid4

import cooler
import h5py
import numpy as np

from tad_errors import TadError


VALID_SIGNALS = {"raw", "balanced"}


class TadArgumentError(ValueError):
    pass


class TadArgumentParser(ArgumentParser):
    def error(self, message):
        raise TadArgumentError(message)


def configuration_context(sample=None, signal=None, resolution=None, run_dir=None):
    return {
        "sample": sample,
        "signal": signal,
        "resolution": resolution,
        "path": Path(run_dir).resolve() if run_dir else None,
    }


def configuration_context_from_argv(argv):
    options = {
        "--sample": "sample",
        "--signal": "signal",
        "--resolution": "resolution",
        "--run-dir": "run_dir",
    }
    values = {}
    for index, argument in enumerate(argv):
        for option, key in options.items():
            if argument == option and index + 1 < len(argv):
                values[key] = argv[index + 1]
            elif argument.startswith(f"{option}="):
                values[key] = argument.split("=", 1)[1]
    return configuration_context(**values)


def configuration_error(context, error):
    return TadError("TAD_CONFIG_INVALID", **context, error=type(error).__name__)



def archive_previous_workspace(run_dir: Path):
    """Move old stage state aside so a prepare cannot reuse derived artifacts."""
    run_dir.mkdir(parents=True, exist_ok=True)
    archive_root = run_dir / ".previous_attempts"
    stale_artifacts = [
        artifact for artifact in run_dir.iterdir()
        if artifact.name != archive_root.name
    ]
    if not stale_artifacts:
        return None
    archive_root.mkdir(exist_ok=True)
    archive_dir = archive_root / f"prepare-{uuid4().hex}"
    archive_dir.mkdir()
    for artifact in stale_artifacts:
        artifact.replace(archive_dir / artifact.name)
    return archive_dir


def parse_preparation_config(args):
    try:
        d_config = json.loads(args.d_config_json)
        if not isinstance(d_config, dict):
            raise ValueError("D configuration must be a JSON object")
        badregionfilter = parse_boolean(args.badregionfilter)
        effective = effective_d_config(d_config, args.resolution)
        for key in (
            "size_bins",
            "variance_thresholds",
            "region_size_bins",
            "overlap_bins",
            "logged",
            "scale",
            "plateau",
            "chaosfilter",
            "chaos_pct",
            "max_gamma",
            "num_part",
            "plots",
            "pctile_threshold",
            "pct_value",
            "final_consistent_domains",
            "size_threshold_bins",
            "boundary_buffer_bins",
        ):
            effective[key]
    except (AttributeError, KeyError, TypeError, ValueError) as error:
        raise configuration_error(
            configuration_context(
                sample=args.sample,
                signal=args.signal,
                resolution=args.resolution,
                run_dir=args.run_dir,
            ),
            error,
        ) from error
    return d_config, badregionfilter


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
            weight1 = weights[int(bin1)]
            weight2 = weights[int(bin2)]
            if (
                not np.isfinite(weight1)
                or not np.isfinite(weight2)
                or weight1 <= 0
                or weight2 <= 0
            ):
                continue
            value = count * weight1 * weight2
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
    parser = TadArgumentParser(exit_on_error=False)
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

    run_dir = Path(args.run_dir).resolve()
    ready_flag = Path(args.ready_flag).resolve()
    ready_flag.unlink(missing_ok=True)
    d_config, badregionfilter = parse_preparation_config(args)
    validate_mcool_resolution(mcool_path, args.resolution, args.sample, args.signal)
    archived_workspace = archive_previous_workspace(run_dir)

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
        d_config=d_config,
        bad_region_bed=args.bad_region_bed,
        badregionfilter=badregionfilter,
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
            f"archived_workspace={archived_workspace or ''}",
        ]
    ) + "\n"
    ready_tmp = ready_flag.with_name(ready_flag.name + ".tmp")
    ready_tmp.write_text(ready_text)
    ready_tmp.replace(ready_flag)


def cli(argv=None):
    arguments = list(sys.argv[1:] if argv is None else argv)
    try:
        args = build_parser().parse_args(arguments)
    except (ArgumentError, TadArgumentError) as error:
        print(
            configuration_error(configuration_context_from_argv(arguments), error),
            file=sys.stderr,
        )
        return 1
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
