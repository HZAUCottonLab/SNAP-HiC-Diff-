"""
05/21/16 Heidi Norton
This module finds a consensus partition from a number of community partitions by computing the adjusted RAND score between all possible pairs of
community partitions. The consensus partition is selected as the partition that has the highest similarity to all other consensus partitions.
Input:
communities: a (number of partitions) x (number of nodes) matrix containing community assignments
Output:
consensus: a consensus partition that is most similar to all other partitions
avg_simm: average similarity across all partitions.
"""
import numpy as np
import matplotlib.pyplot as pyplot
from sklearn.metrics.cluster import adjusted_rand_score

def get_similarity_consensus(communities):
        npart = communities.shape[0]
        #print 'npart', npart
        pairwise_simm = np.ones((npart,npart))

        # Calculate pairwise similarities between each partition
        for i in range(npart):
                for j in range(npart):
                        pairwise_simm[i,j] = adjusted_rand_score(communities[i,:], communities[j,:])

        # make matrix symmetric
        pairwise_simm = (pairwise_simm + np.transpose(pairwise_simm))/2
        #pyplot.imshow(pairwise_simm)
        #pyplot.show()
        #pyplot.clf()


        # compute average pairwise similarity
        avg_pairwise_simm = np.average(pairwise_simm, axis = 0)
        #print "average pairwise similarity",  avg_pairwise_simm


        # find the partition that has the maximal average similarity
        idx = np.argmax(avg_pairwise_simm)
        #print "idx of max average simm", idx

        consensus = communities[idx,:]

        # Compute global average similarity; this represents the 'quality' of the gamma value
        avg_simm = np.average(pairwise_simm)
        #print "global average simm", avg_simm

        return consensus, avg_simm
