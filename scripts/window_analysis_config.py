from dataclasses import dataclass
from math import ceil, isfinite
from pathlib import Path


class WindowAnalysisConfigError(ValueError):
    pass


@dataclass(frozen=True)
class WindowScale:
    resolution_bp: int
    window_bins: int
    window_size_bp: int
    step_bins: int
    step_size_bp: int

    @property
    def tag(self):
        return f"r{self.resolution_bp}/w{self.window_bins}"


@dataclass(frozen=True)
class BackgroundPreset:
    resolution_bp: int
    window_bins: int
    step_bins: int
    file: Path


@dataclass(frozen=True)
class WindowAnalysisSettings:
    resolution_bp: int
    step_fraction: float
    scales: tuple[WindowScale, ...]
    backgrounds_by_window_bins: dict[int, BackgroundPreset]

    def scale(self, window_bins):
        value = int(window_bins)
        return next(scale for scale in self.scales if scale.window_bins == value)


def parse_window_analysis_settings(
    config: dict, check_files: bool = True
) -> WindowAnalysisSettings:
    if "window_analysis" not in config:
        raise WindowAnalysisConfigError(
            "window_analysis.window_bins and window_analysis.step_fraction are required; "
            "legacy window_size/step_size are not supported"
        )
    resolution = config.get("resolution")
    bins = config["window_analysis"].get("window_bins")
    fraction = config["window_analysis"].get("step_fraction")
    if isinstance(resolution, bool) or not isinstance(resolution, int) or resolution <= 0:
        raise WindowAnalysisConfigError("resolution must be a positive integer")
    if not isinstance(bins, list) or not bins or any(
        isinstance(value, bool) or not isinstance(value, int) or value <= 0
        for value in bins
    ):
        raise WindowAnalysisConfigError(
            "window_analysis.window_bins must be a non-empty list of positive integers"
        )
    if len(set(bins)) != len(bins):
        raise WindowAnalysisConfigError("window_analysis.window_bins contains duplicates")
    if (
        isinstance(fraction, bool)
        or not isinstance(fraction, (int, float))
        or not isfinite(float(fraction))
        or fraction <= 0
    ):
        raise WindowAnalysisConfigError(
            "window_analysis.step_fraction must be finite and positive"
        )
    scales = tuple(
        WindowScale(
            resolution,
            window_bins,
            resolution * window_bins,
            ceil(window_bins * float(fraction)),
            resolution * ceil(window_bins * float(fraction)),
        )
        for window_bins in bins
    )
    records = config.get("knn_empirical_pvalue", {}).get("background_presets", [])
    active = {}
    expected_keys = {
        (resolution, scale.window_bins, scale.step_bins) for scale in scales
    }
    for record in records:
        key = (
            record.get("resolution_bp"),
            record.get("window_bins"),
            record.get("step_bins"),
        )
        if key[0] != resolution:
            continue
        if key in expected_keys:
            if key[1] in active:
                raise WindowAnalysisConfigError(
                    f"duplicate background preset for r{resolution}/w{key[1]}"
                )
            path = Path(record.get("file", ""))
            if check_files and not path.is_file():
                raise WindowAnalysisConfigError(
                    f"background file missing for r{resolution}/w{key[1]}: {path}"
                )
            active[key[1]] = BackgroundPreset(*key, path)
    for scale in scales:
        if scale.window_bins not in active:
            raise WindowAnalysisConfigError(
                f"missing background preset for r{resolution}/w{scale.window_bins} "
                f"step_bins={scale.step_bins}"
            )
    return WindowAnalysisSettings(resolution, float(fraction), scales, active)
