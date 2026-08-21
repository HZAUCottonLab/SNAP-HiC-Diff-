
import numpy as np


"""
This module constructs a Newman Girvan null model for a given symmetric input nxn numpy array.
NGij = ki*kj/(2*m) where ki is the degree of node i, kj is the degree of node j, and m is the total sum of edge weights for the network
"""

def construct_NG_null(A):
        NG = np.outer(np.sum(A,1),np.sum(A,1))/np.sum(A)
        return NG
