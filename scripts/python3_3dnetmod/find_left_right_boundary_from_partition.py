import numpy as np
def find_left_right_boundary_from_partition(partition):

        #print 'in find_left_right_boundary_from_partition'
        idx = np.where(partition == 1)

        #print 'idx: '
        #print idx
        #print 'printed idx'
        min = np.min(idx)
        max = np.max(idx)

        return min, max
