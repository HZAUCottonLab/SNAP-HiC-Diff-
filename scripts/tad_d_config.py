from math import isfinite

BOOLEAN_FIELDS = ("logged", "chaosfilter", "plots", "final_consistent_domains")


def _invalid(key):
    raise ValueError(f"invalid tad.d.{key}")


def _integer(values, key, minimum):
    value = values[key]
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        _invalid(key)


def _number(values, key, minimum, maximum=None):
    value = values[key]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _invalid(key)
    if not isfinite(value) or value < minimum:
        _invalid(key)
    if maximum is not None and value > maximum:
        _invalid(key)


def _five_numbers(values, key, minimum, *, integers=False):
    items = values[key]
    if not isinstance(items, list) or len(items) != 5:
        _invalid(key)
    for item in items:
        if isinstance(item, bool) or not isinstance(item, (int, float)):
            _invalid(key)
        if integers and not isinstance(item, int):
            _invalid(key)
        if not isfinite(item) or item < minimum:
            _invalid(key)


def validate_d_config(values):
    if not isinstance(values, dict):
        _invalid("configuration")
    for key, minimum in {
        "region_size_bins": 1,
        "overlap_bins": 0,
        "num_part": 1,
        "size_threshold_bins": 1,
        "boundary_buffer_bins": 0,
    }.items():
        _integer(values, key, minimum)
    if values["overlap_bins"] >= values["region_size_bins"]:
        _invalid("overlap_bins")
    for key in BOOLEAN_FIELDS:
        if not isinstance(values[key], bool):
            _invalid(key)
    scale = values["scale"]
    if not isinstance(scale, str) or not scale.strip():
        _invalid("scale")
    for key in ("plateau", "max_gamma"):
        _number(values, key, 0)
        if values[key] <= 0:
            _invalid(key)
    for key in ("chaos_pct", "pctile_threshold", "pct_value"):
        _number(values, key, 0, 1)
    _five_numbers(values, "size_bins", 1, integers=True)
    _five_numbers(values, "variance_thresholds", 0)
    overrides = values["resolution_overrides"]
    if not isinstance(overrides, dict):
        _invalid("resolution_overrides")
    for resolution, override in overrides.items():
        if str(resolution) != "4000" or not isinstance(override, dict):
            _invalid("resolution_overrides")
        if set(override) != {"size_threshold_bins", "variance_thresholds"}:
            _invalid("resolution_overrides")
        _integer(override, "size_threshold_bins", 1)
        _five_numbers(override, "variance_thresholds", 0)
