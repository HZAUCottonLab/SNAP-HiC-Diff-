"""
Heidi Norton
This module uses the  louvain algorithm to assign communities
Input:
A: a square symmetric matrix of 5C counts data that does not contain NaNs
P: a square symmetric matrix with the same dimensions as A that serves as the null model matrix
num_part: the number of times the louvain algorithm should be applied
resolution: structural resolution parameter used to assign communities within the louvain algorithm
Output:
communities: a numpy matrix with num_part number of rows. Each row represents a community partition.
Q_array: a numpy array with num_part rows containing the Modularity value, Q, for each application of the louvain algorithm

The assign_sorted_communities module modifiesd community assignments such that they are in increasing order.
"""
import numpy as np
from genlouvain import genlouvain


def assign_communities_gen_louvain(A, P, num_part, resolution, seeds, seed_count):
        s = np.sum(A)
        B = (A-resolution*P)/s # B is the modularity matrix
        communities = np.zeros((num_part, len(B))) # each row represents a community partition
        Q_array = np.zeros((num_part))  # Modularity results of each partition will go here
        for i in range(num_part):
                seed = seeds[seed_count]
                [Ci, Q] = genlouvain(B, seed)
                #print Ci
                communities[i,:] = Ci
                Q_array[i] = Q
                if seed_count < len(seeds)-1: # condition to preven running out of seeds in the list
                        seed_count = seed_count + 1
                else:
                        seed_count = 0

        return communities, Q_array, seed_count




def assign_sorted_communities_gen_louvain(A, P, num_part, resolution, seeds, seed_count):
        s = np.sum(A)
        B = (A-resolution*P)/s # B is the modularity matrix
        communities_sorted = np.zeros((num_part, len(B))) # each row represents a community partition
        Q_array = np.zeros((num_part))  # Modularity results of each partition will go here
        for i in range(num_part):
                seed = seeds[seed_count]
                [Ci, Q] = genlouvain(B, seed)
                # Modify Ci such that communities are in increasing order
                value, loc, counts = np.unique(Ci, return_index = True, return_counts = True)
                counts_ordered = counts[np.argsort(loc)]
                Ci_new = np.repeat(list(range(1,len(value)+1)), counts_ordered)
                #print Ci_new
                communities_sorted[i,:] = Ci_new
                Q_array[i] = Q
                if seed_count < len(seeds)-1:
                        seed_count = seed_count + 1
                else:
                        seed_count = 0

        return communities_sorted, Q_array, seed_count




def assign_sorted_communities_gen_louvain_fixed_seed(A, P, num_part, resolution, seed):
        s = np.sum(A)
        B = (A-resolution*P)/s # B is the modularity matrix
        communities_sorted = np.zeros((num_part, len(B))) # each row represents a community partition
        Q_array = np.zeros((num_part))  # Modularity results of each partition will go here
        for i in range(num_part):
                [Ci, Q] = genlouvain(B, seed)
                # Modify Ci such that communities are in increasing order
                value, loc, counts = np.unique(Ci, return_index = True, return_counts = True)
                counts_ordered = counts[np.argsort(loc)]
                Ci_new = np.repeat(list(range(1,len(value)+1)), counts_ordered)
                #print Ci_new
                communities_sorted[i,:] = Ci_new
                Q_array[i] = Q

        return communities_sorted, Q_array
