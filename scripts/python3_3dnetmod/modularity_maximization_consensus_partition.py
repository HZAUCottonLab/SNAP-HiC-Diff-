import copy
from counts_to_array_zero_diagonal import counts_to_array_zero_diagonal
from modularity_maximization import modularity_maximization
from consensus_partition import consensus_partition
from trash_community_removal_avg_num_comm import trash_community_removal
from write_unique_communities_genomewide import write_unique_communities_genomewide
from plot_partition_adj_color import plot_partition
from get_partition_colormap import get_colormap

"""
Heidi Norton
Dan Emerson
Harvey Huang

This script runs the Louvain-like greedy algorithm for a range of resolution parameters (gamma) determined by determine_gamma_range_genome_wide.py
Num_part (typically 20) partitions are determined for each gamma value.
A consensus partition for each hierarchy level is determined.

All user-defined settings should be specified in a settings file

"""
def MMCP(chr,settings,sample, badregions,tag_GPS,tag_MMCP,tag_preprocess):
    #cell_type_input = str(chr) + '_' + settings['sample_1']
    print('in MMCP for : ', chr, sample)
    cell_type_input = str(chr) + '_' + sample
    counts = 'input/' +  cell_type_input + '_' + tag_preprocess + "finalpvalues.counts"
    bed = 'input/' +  cell_type_input + '_' + tag_preprocess + "final.bed"
    #cell_type = str(chr) + settings['sample_1']
    cell_type = str(chr) + sample
    counts_as_arrays = counts_to_array_zero_diagonal(counts,bed, 1)[0] # zero out diagonal
    seed_count = 0
    communities_dict = {}
    consensus_dict = {}

    keys = list(counts_as_arrays.keys())
    print('regions before deleting bad regions', keys)
    for region in badregions:
        if region[1] == chr:
             if region[0] in counts_as_arrays:
                 del counts_as_arrays[region[0]]

    keys = list(counts_as_arrays.keys())
    print('regions after deleting bad regions', keys)
    if 'single_region' in list(settings.keys()):
        if settings['single_region'] in list(counts_as_arrays.keys()):
            regions = [settings['single_region']]
        else:
            print('selected single region not in this chromosome, exiting')
            return 0
    else:
        regions = list(counts_as_arrays.keys())

    if len(list(counts_as_arrays.keys()))> 0:
        #for region in counts_as_arrays:
        for region in regions:
            print('in for loop')
            communities, region, seed_count = modularity_maximization(settings, counts_as_arrays[region], cell_type, region, seed_count,tag_GPS, tag_MMCP, hybrid = 'none')
            communities_dict[region] = communities
        #for region in counts_as_arrays:
        for region in regions:
            print('running consensus partition')
            consensus, region = consensus_partition(settings, communities_dict[region], cell_type, region, tag_MMCP)
            consensus_dict[region] = consensus
            for gamma in communities_dict[region]:
                if settings['plots'] == 'True':
                    filename = cell_type + '_' + tag_MMCP + '_' + region + '_' + str(gamma) + '.png'
                    plot_partition(communities_dict[region][gamma], consensus[gamma], get_colormap(), filename, subdir = 'per_gamma_adjusted_color')
        print('entering trash_communit_removal')
        continuation = trash_community_removal(settings, communities_dict, consensus_dict, counts_as_arrays, cell_type, bed, tag_MMCP, hybrid = 'none')

        if continuation == -1: #exit program, no variances
            print('continuation = -1, exiting')
            return 0
        prefinal = write_unique_communities_genomewide(settings, cell_type,cell_type_input, 'none',tag_MMCP,tag_preprocess)


if __name__ == '__main__':
        main()
