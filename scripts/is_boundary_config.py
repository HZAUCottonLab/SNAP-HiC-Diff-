from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from math import isfinite
from typing import Any

from cworld_is_config import CworldISSettings
from tad_config import TadSettings


@dataclass(frozen=True)
class StrengthSettings:
    local_radius_scale_factor: float
    min_reference_boundaries: int
    min_valid_scales_per_band: int
    robust_z_threshold: float
    stable_abs_z_threshold: float
    sensitivity_thresholds: Sequence[float]


@dataclass(frozen=True)
class PermutationSettings:
    enabled: bool
    n_permutations: int
    seed: int


@dataclass(frozen=True)
class ISBoundarySettings:
    enabled: bool
    association_padding_bp: int
    raw_balanced_overlap_only: bool
    base_bin_bp: int
    jitter_quantile: float
    jitter_min_pairs: int
    jitter_fallback_bp: int
    max_shift_bp: int
    shift_sensitivity_bp: Sequence[int]
    match_tie_tolerance_bp: int
    strength: StrengthSettings
    permutation: PermutationSettings


def _invalid(message: str) -> ValueError:
    return ValueError(f"ERROR[IS_BOUNDARY_CONFIG_INVALID] {message}")


def _required(raw: Mapping[str, Any], name: str) -> Any:
    try:
        return raw[name]
    except KeyError as error:
        raise _invalid(f"missing={name}") from error


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise _invalid(f"{name} must be a mapping")
    return value


def _sequence(value: Any, name: str) -> tuple[Any, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise _invalid(f"{name} must be a sequence")
    return tuple(value)


def _positive_int(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise _invalid(f"{name} must be positive")
    return value


def _positive_float(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise _invalid(f"{name} must be finite and positive")
    try:
        value = float(value)
    except OverflowError as error:
        raise _invalid(f"{name} must be finite and positive") from error
    if not isfinite(value) or value <= 0:
        raise _invalid(f"{name} must be finite and positive")
    return value


def _distance(value: Any, name: str, base: int) -> int:
    value = _positive_int(value, name)
    if value % base:
        raise _invalid(f"{name} must align to the 4 kb grid")
    return value


def parse_is_boundary_settings(
    config: Mapping[str, Any],
    cworld: CworldISSettings,
    tad: TadSettings,
) -> ISBoundarySettings:
    config = _mapping(config, "config")
    raw = _mapping(
        _required(config, "is_boundary_annotation"),
        "is_boundary_annotation",
    )
    strength_raw = _mapping(_required(raw, "strength"), "strength")
    permutation_raw = _mapping(_required(raw, "permutation"), "permutation")

    base = _positive_int(_required(raw, "base_bin_bp"), "base_bin_bp")
    if base != cworld.source_resolution_bp or base != min(tad.resolutions):
        raise _invalid(
            "base_bin_bp must equal the 4 kb IS and finest TAD resolution"
        )
    padding = _required(raw, "association_padding_bp")
    if isinstance(padding, bool) or not isinstance(padding, int) or padding != 0:
        raise _invalid("association_padding_bp must equal zero")
    if _required(raw, "raw_balanced_overlap_only") is not True:
        raise _invalid("raw_balanced_overlap_only must be true")
    quantile = _positive_float(_required(raw, "jitter_quantile"), "jitter_quantile")
    if quantile >= 1:
        raise _invalid("jitter_quantile must be between zero and one")

    sensitivity = _sequence(
        _required(raw, "shift_sensitivity_bp"),
        "shift_sensitivity_bp",
    )
    if not sensitivity:
        raise _invalid("shift_sensitivity_bp must not be empty")
    sensitivity = tuple(
        _distance(value, "shift_sensitivity_bp", base) for value in sensitivity
    )
    if sensitivity != (60000, 80000, 100000):
        raise _invalid("shift_sensitivity_bp must use the fixed publication contract [60000, 80000, 100000]")
    max_shift = _positive_int(_required(raw, "max_shift_bp"), "max_shift_bp")
    if max_shift != 100000:
        raise _invalid("max_shift_bp must equal the fixed primary value 100000")
    if max_shift not in sensitivity:
        raise _invalid("shift_sensitivity_bp must contain max_shift_bp")
    jitter_fallback = _distance(
        _required(raw, "jitter_fallback_bp"),
        "jitter_fallback_bp",
        base,
    )
    match_tie_tolerance = _distance(
        _required(raw, "match_tie_tolerance_bp"),
        "match_tie_tolerance_bp",
        base,
    )
    if jitter_fallback != max(tad.resolutions):
        raise _invalid(
            "jitter_fallback_bp must equal the maximum TAD resolution"
        )

    sensitivity_thresholds = _sequence(
        _required(strength_raw, "sensitivity_thresholds"),
        "sensitivity_thresholds",
    )
    if not sensitivity_thresholds:
        raise _invalid("sensitivity_thresholds must not be empty")
    sensitivity_thresholds = tuple(
        _positive_float(value, "sensitivity_thresholds")
        for value in sensitivity_thresholds
    )
    strength = StrengthSettings(
        local_radius_scale_factor=_positive_float(
            _required(strength_raw, "local_radius_scale_factor"),
            "local_radius_scale_factor",
        ),
        min_reference_boundaries=_positive_int(
            _required(strength_raw, "min_reference_boundaries"),
            "min_reference_boundaries",
        ),
        min_valid_scales_per_band=_positive_int(
            _required(strength_raw, "min_valid_scales_per_band"),
            "min_valid_scales_per_band",
        ),
        robust_z_threshold=_positive_float(
            _required(strength_raw, "robust_z_threshold"),
            "robust_z_threshold",
        ),
        stable_abs_z_threshold=_positive_float(
            _required(strength_raw, "stable_abs_z_threshold"),
            "stable_abs_z_threshold",
        ),
        sensitivity_thresholds=sensitivity_thresholds,
    )
    if strength.robust_z_threshold != 2.0:
        raise _invalid("robust_z_threshold must equal the fixed primary value 2.0")
    if strength.stable_abs_z_threshold >= strength.robust_z_threshold:
        raise _invalid(
            "stable_abs_z_threshold must be smaller than robust_z_threshold"
        )
    if sensitivity_thresholds != (1.5, 2.0, 2.5):
        raise _invalid("sensitivity_thresholds must use the fixed publication contract [1.5, 2.0, 2.5]")

    enabled = _required(raw, "enabled")
    permutation_enabled = _required(permutation_raw, "enabled")
    if not isinstance(enabled, bool) or not isinstance(permutation_enabled, bool):
        raise _invalid("enabled fields must be boolean")
    if permutation_enabled and not enabled:
        raise _invalid(
            "permutation.enabled requires is_boundary_annotation.enabled"
        )
    permutation = PermutationSettings(
        enabled=permutation_enabled,
        n_permutations=_positive_int(
            _required(permutation_raw, "n_permutations"),
            "n_permutations",
        ),
        seed=_positive_int(_required(permutation_raw, "seed"), "seed"),
    )

    return ISBoundarySettings(
        enabled=enabled,
        association_padding_bp=padding,
        raw_balanced_overlap_only=True,
        base_bin_bp=base,
        jitter_quantile=quantile,
        jitter_min_pairs=_positive_int(
            _required(raw, "jitter_min_pairs"),
            "jitter_min_pairs",
        ),
        jitter_fallback_bp=jitter_fallback,
        max_shift_bp=max_shift,
        shift_sensitivity_bp=sensitivity,
        match_tie_tolerance_bp=match_tie_tolerance,
        strength=strength,
        permutation=permutation,
    )
