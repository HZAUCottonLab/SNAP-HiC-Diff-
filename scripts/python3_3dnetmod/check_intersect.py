def check_intersect(a, b):
    """
    Checks to see if two features intersect.

    Parameters
    ---------
    a, b : Dict[str, Any]
        The two features to check for intersection.

    Returns
    -------
    bool
        True if the features intersect, False otherwise.

    Notes
    -----
    Features are represented as dicts with the following structure::

            {
                'chrom': str
                'start': int,
                'end'  : int,
            }

    See ``lib5c.parsers.bed.load_features()``.
    """
    if (
            a['chrom'] == b['chrom'] and
            a['end'] > b['start'] and
            b['end'] > a['start']
    ):
        return True
    return False
