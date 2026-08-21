import os
from outdir_helper import OUTDIR
from read_seeds import read_seeds
import numpy as np
from construct_nulls import construct_NG_null
from assign_communities_gen_louvain import assign_sorted_communities_gen_louvain
from plot_partition import plot_partition
from make_scatterplots import make_scatterplots
from get_partition_colormap import get_colormap
from get_gammas import get_gammas
from convert_tag_genomewide import convert_tag_genomewide

def modularity_maximization(settings, counts_array, cell_type, region,seed_count, tag_GPS,tag_MMCP, hybrid = 'none'):
        """
        Creates partition block for counts array at every gamma in a range
        Saves image and text files for partition blocks at each gamma
        Makes a scatterplot of avg. num of communities vs gamma

        INPUTS
        ------
        settings, dict of strings
                dictionary corresponding settings categories to their string values
        counts_array, ndarray
                a symmetrical array of counts
        region, str
                the region of the array
        seed_count, int
                the current counter of seed location

        OUTPUTS
        -------
        communities_dict, dict of ndarrays
                key = gamma, value = num_part x num_nodes 2D partition block
        seed_count, int
                an updated counter of seed location
        """
        print('in modularity maximization inner script for: ', cell_type, region)
        print('\npartioning %s region' %region)
        assert not np.any(np.isnan(counts_array)) # check for nans
        tag_GPS_test = convert_tag_genomewide(tag_GPS, settings)
        if 'single_gamma' in list(settings.keys()):
                gammas = [float(settings['single_gamma'])]
        else:
                gamma_file = OUTDIR + '/GPS/gamma_dynamic_range/gamma_dynamic_range_' + cell_type + '_' + region  + '_' + tag_GPS
                gamma_file2 = OUTDIR + '/GPS/gamma_dynamic_range/gamma_dynamic_range_' + cell_type + '_' + region  + '_' + tag_GPS_test

                if os.path.isfile(gamma_file):
                        gammas = get_gammas(gamma_file)
                elif os.path.isfile(gamma_file2):
                        gammas = get_gammas(gamma_file2)
                else:
                        print('could not find gamma file')
                        return None
        print('gammas: ', gammas)


        communities_dict = {}
        NG = construct_NG_null(counts_array) # null model
        num_comm = []
        if 'num_part' not in settings:
                settings['num_part'] = 20 #default being 20

        if 'seed_file' in list(settings.keys()):
                gen_louvain_seeds = read_seeds(settings['seed_file'], settings['num_part']) # a list of int seeds
        else:
                gen_louvain_seeds = read_seeds('gen_louvain_100_seeds.txt', settings['num_part'])
        for gamma in gammas:
                print('gamma %.2f being partitioned' %gamma)

                # Assign communities using the Louvain algorithm numpart times.
                # Ci is a num_part x num_nodes 2D numpy array of community assignments. Each row represents a different community partition.
                # Q is a list of modularity values from each application of the algorithm
                [Ci, Q, seed_count] = assign_sorted_communities_gen_louvain(counts_array, NG, int(settings['num_part']), gamma, gen_louvain_seeds, seed_count)
                print('ran genlouvain')

                # Add partition block to dictionary
                communities_dict[gamma] = Ci

                # add mean number of communities at current gamma
                mean_number_of_communities = round(np.mean(Ci.max(axis=1)))
                num_comm.append(mean_number_of_communities)

                identifier = cell_type + '_' + region + '_' + tag_MMCP
                # Save image of partition across num_part for a given gamma
                if hybrid == 'none':
                        if eval(settings["plots"]):
                                filename = 'Partitions_across_num_part_' + identifier + '_' + str(gamma) + '.png'
                                plot_partition(Ci, get_colormap(), filename, subdir='per_gamma')

        if hybrid == 'none':
                filename = 'Avg_num_comm_v_gamma_' + identifier + '.png'

        return communities_dict, region, seed_count
