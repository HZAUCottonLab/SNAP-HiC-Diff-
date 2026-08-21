from dataclasses import dataclass


class TadConfigError(ValueError):
    pass


REQUIRED_D_KEYS = {
    "region_size_bins", "overlap_bins", "logged", "scale", "plateau",
    "chaosfilter", "chaos_pct", "max_gamma", "num_part", "plots",
    "pctile_threshold", "pct_value", "final_consistent_domains",
    "size_threshold_bins", "boundary_buffer_bins", "size_bins",
    "variance_thresholds", "resolution_overrides",
}


@dataclass(frozen=True)
class TadSettings:
    resolutions: tuple
    signals: tuple
    threads: int
    parameter_set: str
    run_root: str
    code_dir: str
    badregionfilter: bool
    bad_region_bed: object
    excluded_chromosomes: tuple
    d: dict


def parse_tad_settings(config):
    try:
        return _parse_tad_settings(config)
    except TadConfigError:
        raise
    except Exception as error:
        raise TadConfigError(
            f"ERROR[TAD_CONFIG_INVALID] invalid={type(error).__name__}"
        ) from error


def _parse_tad_settings(config):
    try:
        section = config["tad"]
        resolutions = tuple(section["resolutions"])
        signals = tuple(section["signals"])
        threads = section["threads"]
        parameter_set = section["parameter_set"]
        d = section["d"]
    except (KeyError, TypeError) as error:
        raise TadConfigError(f"ERROR[TAD_CONFIG_INVALID] missing={error}") from error

    if set(resolutions) != {4000, 8000, 20000} or len(resolutions) != 3:
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] tad.resolutions must be 4000,8000,20000")
    if any(resolution not in config["convert"]["resolutions"] for resolution in resolutions):
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] tad.resolutions require convert.resolutions")
    if set(signals) != {"raw", "balanced"} or len(signals) != 2:
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] tad.signals must be raw,balanced")
    if isinstance(threads, bool) or not isinstance(threads, int) or threads <= 0:
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] tad.threads must be positive")
    if parameter_set != "D":
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] tad.parameter_set must be D")
    missing = REQUIRED_D_KEYS.difference(d)
    if missing:
        raise TadConfigError(
            "ERROR[TAD_CONFIG_INVALID] tad.d missing=" + ",".join(sorted(missing))
        )
    if len(d["size_bins"]) != 5 or len(d["variance_thresholds"]) != 5:
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] D arrays must contain five values")
    excluded = section.get("excluded_chromosomes", [])
    if not isinstance(excluded, list):
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] excluded_chromosomes must be a list")
    badregionfilter = section.get("badregionfilter", False)
    bad_region_bed = section.get("bad_region_bed")
    if badregionfilter and not bad_region_bed:
        raise TadConfigError("ERROR[TAD_CONFIG_INVALID] bad_region_bed is required")
    return TadSettings(
        resolutions=tuple(sorted(resolutions)),
        signals=tuple(signal for signal in ("raw", "balanced") if signal in signals),
        threads=threads,
        parameter_set=parameter_set,
        run_root=section["run_root"],
        code_dir=section["code_dir"],
        badregionfilter=badregionfilter,
        bad_region_bed=bad_region_bed,
        excluded_chromosomes=tuple(excluded),
        d=d,
    )
