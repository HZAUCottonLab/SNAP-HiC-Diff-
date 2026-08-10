from dataclasses import dataclass
from typing import Any, Mapping


EXPECTED_SCALE_COUNT = 10
MELTRON_MIN_RESOLUTION_BP = 1000


class ConfigError(ValueError):
    pass


@dataclass(frozen=True)
class CworldISSettings:
    source_resolution_bp: int
    is_distances_bp: tuple[int, ...]
    ids_divisor: int
    chrom_job_threads: int
    merge_job_threads: int

    @property
    def ids_distances_bp(self) -> tuple[int, ...]:
        return tuple(distance // self.ids_divisor for distance in self.is_distances_bp)


def parse_cworld_is_settings(config: Mapping[str, Any]) -> CworldISSettings:
    try:
        raw = config["meltron"]["cworld_is"]
    except (KeyError, TypeError) as exc:
        raise ConfigError("meltron.cworld_is is required") from exc
    if not isinstance(raw, Mapping):
        raise ConfigError("meltron.cworld_is must be a mapping")

    keys = (
        "source_resolution_bp", "is_distances_bp", "ids_divisor",
        "chrom_job_threads", "merge_job_threads",
    )
    missing = [key for key in keys if key not in raw]
    if missing:
        raise ConfigError(f"meltron.cworld_is missing keys: {', '.join(missing)}")

    resolution = raw["source_resolution_bp"]
    distances = raw["is_distances_bp"]
    divisor = raw["ids_divisor"]
    chrom_threads = raw["chrom_job_threads"]
    merge_threads = raw["merge_job_threads"]
    scalar_values = (resolution, divisor, chrom_threads, merge_threads)
    if any(isinstance(value, bool) or not isinstance(value, int) for value in scalar_values):
        raise ConfigError("resolution, ids_divisor, and thread counts must be integers")
    if resolution <= 0:
        raise ConfigError("source_resolution_bp must be positive")
    if resolution < MELTRON_MIN_RESOLUTION_BP:
        raise ConfigError(
            "source_resolution_bp must be at least 1000 bp for Meltron"
        )
    if divisor != 2:
        raise ConfigError("ids_divisor must be 2 so --ids equals --is / 2")
    if chrom_threads <= 0 or merge_threads <= 0:
        raise ConfigError("thread counts must be positive")
    if not isinstance(distances, list) or any(
        isinstance(value, bool) or not isinstance(value, int) for value in distances
    ):
        raise ConfigError("is_distances_bp must be a list of integers")
    if len(distances) != EXPECTED_SCALE_COUNT:
        raise ConfigError("is_distances_bp must contain exactly 10 values")
    if len(set(distances)) != len(distances):
        raise ConfigError("is_distances_bp values must be unique")
    if any(
        value <= 0 or value % resolution or value % divisor for value in distances
    ):
        raise ConfigError(
            "every IS distance must be positive and divisible by source_resolution_bp and ids_divisor"
        )
    return CworldISSettings(
        resolution, tuple(distances), divisor, chrom_threads, merge_threads
    )
