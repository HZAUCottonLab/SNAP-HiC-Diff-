
import itertools
import numpy as np
"""
This module zeroes out communities that are less than an input size_threshold.
[1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5, 6, 6]
with a size threshold of 2 becomes:
[1, 1, 1, 0, 0, 3, 3, 3, 3, 0, 5, 5, 5, 0, 0]
"""
def size_threshold_array(array, size_threshold):
        grouped = [list(g) for k, g in itertools.groupby(array)]
        thresholded_grouped = []
        for item in grouped:
                if len(item)<=size_threshold:
                        item = [0]*len(item)
                thresholded_grouped.append(item)

        thresholded_array = [x for sublist in thresholded_grouped for x in sublist]

        return np.asarray(thresholded_array)
