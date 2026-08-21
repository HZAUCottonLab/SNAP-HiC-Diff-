import os
from outdir_helper import OUTDIR
import numpy as np
from get_similarity_consensus import get_similarity_consensus
import time

def consensus_partition(settings, level_blocks, cell_type, region, tag_MMCP):
        """
        Finds similarity consensus of levels and saves a text output

        INPUTS
        ------
        settings, dict of strings
                dictionary corresponding settings categories to their string values
        level_blocks, dict of ndarrays
                key = level, value = ndarray of num_part*(number of gammas in level) by num_node community partitions per level
        region, str
                the region that the levels came from

        OUTPUTS
        -------
        none

        consensus_dict, dict of arrays
                key = level, value = 1 x n array where n is number of nodes representing the consensus partition per level
        """

        print("\ngenerating similarity consensus partitions in region %s" %region)
        consensus_dict = {}
        for avg_num_comm in level_blocks:
                level_consensus = get_similarity_consensus(level_blocks[avg_num_comm])[0]
                if settings['plots'] == 'True':
                        dir = OUTDIR + '/MMCP/consensus_files/'
                        if not os.path.isdir(dir): os.makedirs(dir)
                        filename = 'Consensus_' + tag_MMCP + '_' + cell_type + '_' + region + '_' + str(avg_num_comm) + '.txt'
                        np.savetxt(os.path.join(dir,filename), level_consensus)
                consensus_dict[avg_num_comm] = level_consensus

        return consensus_dict, region
