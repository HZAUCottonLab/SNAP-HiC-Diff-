"""Settings-driven worker limits for vendored 3DNetMod stages."""

def worker_count(settings, work_items):
    """Clamp a multiprocessing pool to the configured processors and work."""
    try:
        processors = int(settings["processors"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("settings.processors must be a positive integer") from error
    if processors < 1:
        raise ValueError("settings.processors must be a positive integer")
    try:
        items = int(work_items)
    except (TypeError, ValueError) as error:
        raise ValueError("work_items must be an integer") from error
    return min(processors, max(1, items))
