"""Validate bin-derived window scales and exact KNN background presets."""

from math import ceil, isfinite
from typing import Any, Mapping


def _positive_integer(value, label):
    if type(value) is not int or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def build_window_scales(
    config: Mapping[str, Any],
) -> dict[str, dict[str, int | str]]:
    for legacy_key in ("window_size", "step_size"):
        if legacy_key in config:
            raise ValueError(f"remove legacy config key: {legacy_key}")

    resolution = _positive_integer(config.get("resolution"), "resolution")
    raw_windows = config.get("window_analysis")
    if not isinstance(raw_windows, Mapping):
        raise ValueError("window_analysis must be a mapping")

    window_bins_values = raw_windows.get("window_bins")
    if not isinstance(window_bins_values, list) or not window_bins_values:
        raise ValueError("window_bins must be a non-empty list")
    if any(type(value) is not int or value <= 0 for value in window_bins_values):
        raise ValueError("window_bins must contain positive integers")
    if len(window_bins_values) != len(set(window_bins_values)):
        raise ValueError("window_bins must not contain duplicates")

    step_fraction = raw_windows.get("step_fraction")
    if isinstance(step_fraction, bool) or not isinstance(step_fraction, (int, float)):
        raise ValueError("step_fraction must be numeric")
    if not isfinite(step_fraction):
        raise ValueError("step_fraction must be finite")
    if not 0 < step_fraction <= 1:
        raise ValueError("step_fraction must be greater than 0 and at most 1")

    raw_presets = config.get("knn_empirical_pvalue", {}).get("background_presets")
    if not isinstance(raw_presets, list) or not raw_presets:
        raise ValueError("background_presets must be a non-empty list")

    active_presets = {}
    for preset in raw_presets:
        if not isinstance(preset, Mapping):
            raise ValueError("each background preset must be a mapping")
        preset_resolution = _positive_integer(
            preset.get("resolution_bp"), "background resolution_bp"
        )
        if preset_resolution != resolution:
            continue
        preset_bins = _positive_integer(
            preset.get("window_bins"), "background window_bins"
        )
        preset_step = _positive_integer(
            preset.get("step_bins"), "background step_bins"
        )
        path = preset.get("file")
        if not isinstance(path, str) or not path.strip():
            raise ValueError("background preset file must be a non-empty string")
        key = (preset_resolution, preset_bins, preset_step)
        if key in active_presets:
            raise ValueError(f"duplicate background preset: {key}")
        active_presets[key] = path

    expected = [
        (resolution, window_bins, ceil(window_bins * step_fraction))
        for window_bins in window_bins_values
    ]
    unused = sorted(set(active_presets).difference(expected))
    if unused:
        raise ValueError(f"unused background preset for active resolution: {unused[0]}")

    scales = {}
    for key in expected:
        _, window_bins, step_bins = key
        if key not in active_presets:
            raise ValueError(
                "missing background preset for "
                f"resolution_bp={resolution}, window_bins={window_bins}, "
                f"step_bins={step_bins}"
            )
        scales[f"w{window_bins}"] = {
            "resolution_bp": resolution,
            "window_bins": window_bins,
            "window_size_bp": resolution * window_bins,
            "step_bins": step_bins,
            "step_size_bp": resolution * step_bins,
            "background_file": active_presets[key],
        }
    return scales
