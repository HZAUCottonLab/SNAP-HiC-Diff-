import numpy as np

def flatten_counts_single_region(regional_counts, discard_nan=False):
    """
    Flattens the counts for a single region into a flat, nonredundant list.

    Parameters
    ---------
    regional_counts : 2d numpy array
        The square, symmetric array of counts for one region.
    discard_nan : bool
        If True, nan's will not be present in the returned list.

    Returns
    -------
    list of floats
        A flat, nonredundant lists of counts. The ``(i, j)`` th element of the
        ``regional_counts`` array (for ``i >= j`` ) ends up at the
        ``(i*(i+1)/2 + j)`` th index of the flattened list. If ``discard_nan``
        was True, these indices will not necessarily match up and it will not be
        possible to unflatten the list.
    """
    flattened_list = [regional_counts[i, j]
                      for i in range(len(regional_counts))
                      for j in range(i + 1)]
    if discard_nan:
        flattened_list = list(filter(np.isfinite, flattened_list))
    return flattened_list
