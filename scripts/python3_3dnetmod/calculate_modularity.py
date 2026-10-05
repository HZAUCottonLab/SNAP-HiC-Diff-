"""
This module uses the  louvain algorithm to calculate modularity
Input:
A: a square symmetric matrix of 5C counts data that does not contain NaNs
P: a square symmetric matrix with the same dimensions as A that serves as the null model matrix
resolution: structural resolution parameter used to assign communities within the louvain algorithm
seed: an integer seed used to set the random state of the louvain algorithm
Output:
Q, a scalar value of the modularity of the network at the given structural resolution parameter

"""
import numpy as np
from genlouvain import genlouvain


def calculate_modularity(A, P, resolution, seed):
        s = np.sum(A)
        B = (A-resolution*P)/s # B is the modularity matrix
        #print 'mean of modularity matrix', np.mean(np.mean(B))
        #print 'trace of modularity matrix', np.trace(B)
        [Ci, Q] = genlouvain(B, seed)

        return Q, Ci
