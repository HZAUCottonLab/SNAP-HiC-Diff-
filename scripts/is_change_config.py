from dataclasses import dataclass
from typing import Any, Mapping

from cworld_is_config import ConfigError, CworldISSettings


_BAND_NAMES = ("small", "medium", "large")
_SCALAR_KEYS = (
    "raw_knn_pvalue_threshold",
    "meltron_score_threshold",
    "merge_gap_bp",
    "context_bp",
    "local_radius_bins",
    "min_local_valid_fraction",
    "flank_fraction_of_scale",
    "min_flank_bins",
    "min_flank_valid_fraction",
    "max_missing_run_bins",
    "high_percentile",
    "broad_percentile",
    "spatial_cluster_tolerance_bp",
    "min_primary_scale_bands",
    "diagnostic_plot_top_n",
)
_FRACTION_KEYS = (
    "raw_knn_pvalue_threshold",
    "min_local_valid_fraction",
    "flank_fraction_of_scale",
    "min_flank_valid_fraction",
    "high_percentile",
    "broad_percentile",
)
_POSITIVE_INTEGER_KEYS = (
    "context_bp",
    "local_radius_bins",
    "min_flank_bins",
    "max_missing_run_bins",
    "spatial_cluster_tolerance_bp",
    "min_primary_scale_bands",
    "diagnostic_plot_top_n",
)


@dataclass(frozen=True)
class ISChangeSettings:
    resolution_bp: int
    scales_bp: tuple
    scale_bands: dict[str, tuple]
    raw_knn_pvalue_threshold: float
    meltron_score_threshold: float
    merge_gap_bp: int
    context_bp: int
    local_radius_bins: int
    min_local_valid_fraction: float
    flank_fraction_of_scale: float
    min_flank_bins: int
    min_flank_valid_fraction: float
    max_missing_run_bins: int
    high_percentile: float
    broad_percentile: float
    spatial_cluster_tolerance_bp: int
    min_primary_scale_bands: int
    diagnostic_plot_top_n: int


def _is_number(value: Any) -> bool:
    return not isinstance(value, bool) and isinstance(value, (int, float))


def validate_is_change_config(
    raw: Mapping[str, Any], cworld_settings: CworldISSettings
) -> None:
    if not isinstance(raw, Mapping):
        raise ConfigError("is_change must be a mapping")
    missing = [key for key in (*_SCALAR_KEYS, "scale_bands") if key not in raw]
    if missing:
        raise ConfigError(f"is_change missing keys: {', '.join(missing)}")

    for key in _FRACTION_KEYS:
        value = raw[key]
        if not _is_number(value) or not 0 < value < 1:
            raise ConfigError(f"{key} must be strictly between zero and one")
    if raw["min_flank_valid_fraction"] < 0.80:
        raise ConfigError("min_flank_valid_fraction must be at least 0.80")
    if not _is_number(raw["meltron_score_threshold"]):
        raise ConfigError("meltron_score_threshold must be numeric")
    if isinstance(raw["merge_gap_bp"], bool) or raw["merge_gap_bp"] != 0:
        raise ConfigError("merge_gap_bp must equal zero")

    for key in _POSITIVE_INTEGER_KEYS:
        value = raw[key]
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ConfigError(f"{key} must be a positive integer")
    if raw["min_primary_scale_bands"] != 2:
        raise ConfigError("min_primary_scale_bands must equal exactly two")
    if raw["context_bp"] < max(cworld_settings.is_distances_bp):
        raise ConfigError("context_bp must cover the largest scale")
    if raw["high_percentile"] <= raw["broad_percentile"]:
        raise ConfigError("high_percentile must exceed broad_percentile")
    if raw["spatial_cluster_tolerance_bp"] % cworld_settings.source_resolution_bp:
        raise ConfigError(
            "spatial_cluster_tolerance_bp must be divisible by resolution_bp"
        )

    bands = raw["scale_bands"]
    if not isinstance(bands, Mapping):
        raise ConfigError("scale_bands must be a mapping")
    if set(bands) != set(_BAND_NAMES):
        raise ConfigError("scale_bands must partition scales into small, medium, and large")
    if any(not isinstance(bands[name], list) for name in _BAND_NAMES):
        raise ConfigError("scale_bands values must be lists")
    flattened = tuple(scale for name in _BAND_NAMES for scale in bands[name])
    if flattened != cworld_settings.is_distances_bp:
        raise ConfigError("scale_bands must partition cworld IS distances exactly")


def parse_is_change_settings(
    config: Mapping[str, Any], cworld_settings: CworldISSettings
) -> ISChangeSettings:
    try:
        raw = config["is_change"]
    except (KeyError, TypeError) as exc:
        raise ConfigError("is_change is required") from exc
    validate_is_change_config(raw, cworld_settings)
    return ISChangeSettings(
        resolution_bp=cworld_settings.source_resolution_bp,
        scales_bp=tuple(cworld_settings.is_distances_bp),
        scale_bands={name: tuple(scales) for name, scales in raw["scale_bands"].items()},
        **{key: raw[key] for key in _SCALAR_KEYS},
    )
