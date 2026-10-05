
import numpy as np

# Thomas Gilgenast

def compute_regional_spatial_variance(regional_counts, window_width, ddof=1):
    n = len(regional_counts)  # the size of the region
    w = window_width  # the width of the window
    r = (w-1)/2  # a "radius" for the window

    variances = []  # list to store variances

    for i in range(3*r, n-r):
        for j in range(r, i-2*r+1):
            counts_slice = regional_counts[j-r:j+r+1, i-r:i+r+1]
            if not np.any(~np.isfinite(counts_slice)):
                variances.append(np.var(counts_slice, ddof=ddof))

    return np.mean(variances)
