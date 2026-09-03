from dataclasses import dataclass
from typing import Any, Mapping, Sequence

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


def _positive_int(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"ERROR[IS_BOUNDARY_CONFIG_INVALID] {name} must be positive")
    return value


def parse_is_boundary_settings(
    config: Mapping[str, Any],
    cworld: CworldISSettings,
    tad: TadSettings,
) -> ISBoundarySettings:
    try:
        raw = config["is_boundary_annotation"]
        strength_raw = raw["strength"]
        permutation_raw = raw["permutation"]
    except (KeyError, TypeError) as error:
        raise ValueError(
            f"ERROR[IS_BOUNDARY_CONFIG_INVALID] missing={error}"
        ) from error
    base = _positive_int(raw["base_bin_bp"], "base_bin_bp")
    if base != cworld.source_resolution_bp or base != min(tad.resolutions):
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] base_bin_bp must equal the 4 kb IS and finest TAD resolution"
        )
    if raw["association_padding_bp"] != 0:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] association_padding_bp must equal zero"
        )
    if raw["raw_balanced_overlap_only"] is not True:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] raw_balanced_overlap_only must be true"
        )
    quantile = raw["jitter_quantile"]
    if isinstance(quantile, bool) or not isinstance(quantile, (int, float)) or not 0 < quantile < 1:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] jitter_quantile must be between zero and one"
        )
    sensitivity = tuple(raw["shift_sensitivity_bp"])
    if sensitivity != tuple(sorted(set(sensitivity))):
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] shift_sensitivity_bp must be sorted and unique"
        )
    max_shift = _positive_int(raw["max_shift_bp"], "max_shift_bp")
    if max_shift not in sensitivity:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] shift_sensitivity_bp must contain max_shift_bp"
        )
    grid_values = sensitivity + (
        raw["jitter_fallback_bp"], raw["match_tie_tolerance_bp"], max_shift,
    )
    if any(value % base for value in grid_values):
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] distance values must align to the 4 kb grid"
        )
    strength = StrengthSettings(
        local_radius_scale_factor=float(strength_raw["local_radius_scale_factor"]),
        min_reference_boundaries=_positive_int(
            strength_raw["min_reference_boundaries"], "min_reference_boundaries"
        ),
        min_valid_scales_per_band=_positive_int(
            strength_raw["min_valid_scales_per_band"], "min_valid_scales_per_band"
        ),
        robust_z_threshold=float(strength_raw["robust_z_threshold"]),
        stable_abs_z_threshold=float(strength_raw["stable_abs_z_threshold"]),
        sensitivity_thresholds=tuple(strength_raw["sensitivity_thresholds"]),
    )
    if strength.local_radius_scale_factor <= 0:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] local_radius_scale_factor must be positive"
        )
    if strength.stable_abs_z_threshold >= strength.robust_z_threshold:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] stable_abs_z_threshold must be smaller than robust_z_threshold"
        )
    if (
        tuple(strength.sensitivity_thresholds)
        != tuple(sorted(set(strength.sensitivity_thresholds)))
        or strength.robust_z_threshold not in strength.sensitivity_thresholds
    ):
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] sensitivity_thresholds must be sorted, unique, and contain robust_z_threshold"
        )
    permutation = PermutationSettings(
        enabled=permutation_raw["enabled"],
        n_permutations=_positive_int(
            permutation_raw["n_permutations"], "n_permutations"
        ),
        seed=_positive_int(permutation_raw["seed"], "seed"),
    )
    if not isinstance(raw["enabled"], bool) or not isinstance(permutation.enabled, bool):
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] enabled fields must be boolean"
        )
    if permutation.enabled and not raw["enabled"]:
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] permutation.enabled requires is_boundary_annotation.enabled"
        )
    jitter_fallback = _positive_int(raw["jitter_fallback_bp"], "jitter_fallback_bp")
    if jitter_fallback != max(tad.resolutions):
        raise ValueError(
            "ERROR[IS_BOUNDARY_CONFIG_INVALID] jitter_fallback_bp must equal the maximum TAD resolution"
        )
    return ISBoundarySettings(
        enabled=raw["enabled"],
        association_padding_bp=raw["association_padding_bp"],
        raw_balanced_overlap_only=True,
        base_bin_bp=base,
        jitter_quantile=float(quantile),
        jitter_min_pairs=_positive_int(raw["jitter_min_pairs"], "jitter_min_pairs"),
        jitter_fallback_bp=jitter_fallback,
        max_shift_bp=max_shift,
        shift_sensitivity_bp=sensitivity,
        match_tie_tolerance_bp=_positive_int(
            raw["match_tie_tolerance_bp"], "match_tie_tolerance_bp"
        ),
        strength=strength,
        permutation=permutation,
    )
