import os
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import glob
from get_partition_colormap import get_colormap
from create_bin_location_dict import create_bin_location_dict
from plot_consensus_partition import plot_consensus_partition
from create_results_output_genomewide import create_results_output_genomewide
from make_scatterplots import make_scatterplots
from get_variance import get_variance_about_boundaries
from get_variance import get_variance_about_boundaries_size_thresholded
from read_seeds import read_seeds
from plot_partition import plot_partition
#import pandas as pd
from size_threshold_array import size_threshold_array
from flatten_counts_single_region_geometric import flatten_counts_single_region
from create_results_output_genomewide_remove_trash import create_results_output_genomewide_remove_trash
from outdir_helper import OUTDIR, consensus_filtered_filename



def remove_low_count_communities(array, counts, pct_threshold):
        print('array pre low count removal: ')
        print(array)


        threshold = np.percentile(flatten_counts_single_region(counts, discard_nan=True),pct_threshold)
        grouped = [list(g) for k, g in itertools.groupby(array)]
        for item in grouped:
                if item[0] !=0:
                        idx = np.where(array==item[0])
                        min_bin = min(idx[0])
                        max_bin = max(idx[0]) + 1
                        comm_column = counts[min_bin:max_bin, max_bin]
                        pct_value = float(len(comm_column[comm_column>threshold]))/float(len(comm_column))
                        if pct_value <=0.20:
                                array[array==item[0]] = 0

        print('adjusted array: ')

        return array

def trash_community_removal(settings, level_blocks, level_consensus, counts_as_arrays, cell_type, bed_file, tag, hybrid = 'none'):
        """
        This module takes an input set of partitions and consensus partitions and writes a results file of communities that pass the size and variance threshold. Typically, for the genome-wide pipeline, the variance threshold is 100AUC (no variance thresholding) at this stage. Thus, when used with the recommended settings, this script writes a results file of all communities that pass the user-defined size threshold
        """
        print('in trash_community_removal')
        bin_location = create_bin_location_dict(bed_file) #contains input directory in name

        my_cmap = get_colormap()
        level_variances = {}
        level_variances_presizethreshold = {}
        level_distances_squared = {}
        level_boundaries = {}
        level_consensus_zeroed = {}
        level_consensus_size_thresholded = {}
        variance_list = {}
        dist_list_squared = {}
        for region in level_blocks:
                level_consensus_zeroed[region] = {} # dictionary with key = level, value = similarity consensus after zeroing out communities
                level_consensus_size_thresholded[region] = {} # dictionary with key = level, value = similarity consensus after zeroing out variance AND size thresholded communities
                level_variances[region] = {} # dictionary with key = level, value = variance corresponding to the boundaries in similarity consensus
                level_variances_presizethreshold[region] = {}
                level_boundaries[region] = {} # dictionary with key = level, value = list of boundaries
                level_distances_squared[region] = {} # dictionary with key = level



        for region in level_blocks:
                filename_qualifier = cell_type + '_' + region
                #print 'REGION', region
                # Stores boundary variances to be plotted in histogram
                variance_list[region] = []
                dist_list_squared[region] = []
                for avg_num_comm in level_blocks[region]:
                        num_nodes = len(level_consensus[region][avg_num_comm])
                        # First index of each community for this consensus in a list
                        boundaries_pre_size_threshold = np.where(level_consensus[region][avg_num_comm][1:] != level_consensus[region][avg_num_comm][:-1])[0]+1

                        bound_vars, distances, distances_squared, boundaries = get_variance_about_boundaries_size_thresholded(level_blocks[region][avg_num_comm], boundaries_pre_size_threshold, 0, type='regular')

                        bound_vars_nonsizethresholded, x, y  = get_variance_about_boundaries(level_blocks[region][avg_num_comm], boundaries_pre_size_threshold, type='regular')
                        # concatenate boundary variances to list to be plotted (ignore the first and last entries)
                        variance_list[region] += bound_vars[1:-1] # if return_edges = true
                        #variance_list[region]+=bound_vars
                        dist_list_squared[region]+=distances_squared
                        # Store corresponding boundaries and their variance values into dictionary by vertically stacking
                        # Add left and right edges to boundaries before concatenating
                        level_variances[region][avg_num_comm] = np.vstack((np.append(np.insert(boundaries, 0, 0), len(level_consensus[region][avg_num_comm])), np.array(bound_vars))) # if return_edges = True
                        level_variances_presizethreshold[region][avg_num_comm]=np.vstack((np.append(np.insert(boundaries_pre_size_threshold, 0, 0), len(level_consensus[region][avg_num_comm])), np.array(bound_vars_nonsizethresholded)))


                        level_boundaries[region][avg_num_comm] = boundaries
                        level_distances_squared[region][avg_num_comm] = distances_squared # used for making violin plots





        two_sided_variance_distribution_sum = []
        two_sided_variance_distribution_max = []
        for region in level_variances:
                for level in level_variances[region]:
                        for i in range(len(level_variances[region][level][1])-1):
                                maximum = max(level_variances[region][level][1,i], level_variances[region][level][1,i+1])
                                #print 'maximum', maximum
                                two_sided_variance_distribution_max.append(maximum)
        #print 'percentage non-zero'
        try:
                float(np.count_nonzero(two_sided_variance_distribution_max))/len(two_sided_variance_distribution_max)
        except ZeroDivisionError as e:
                z = e
                print(z)
                return -1
        count = 0
        for item in two_sided_variance_distribution_max:
                if item <1:
                        count +=1


        # Handle edge case of nan values in variance distribution
        two_sided_variance_distribution_max = np.array(two_sided_variance_distribution_max)
        max_variance = np.nanmax(two_sided_variance_distribution_max) # Find the maximal value in the distribution excluding nans
        idx = np.isnan(two_sided_variance_distribution_max)
        two_sided_variance_distribution_max[idx] = max_variance + 1 # Set nans to be the maximum + 1, thus ensuring they are the greatest value in the distribution

        for region in level_variances:
                for gamma in level_variances[region]:
                        idx = np.isnan(level_variances[region][gamma])
                        level_variances[region][gamma][idx] = max_variance + 1
                        idx = np.isnan(level_variances_presizethreshold[region][gamma])
                        level_variances_presizethreshold[region][gamma][idx] = max_variance + 1






        #variance_threshold = global_variance_threshold_max
        variance_threshold = np.percentile(two_sided_variance_distribution_max, 100) # at this stage, we keep all communities
        for region in level_blocks:
                community_coord_list = []
                filename_qualifier = cell_type + '_' + region + '_variance_threshold_two_sided_max_100' + '_AUC_' + tag

                for avg_num_comm in level_blocks[region]:
                        # Make a copy of consensus to zero out
                        level_consensus_zeroed[region][avg_num_comm] = level_consensus[region][avg_num_comm].copy()

                        # Iterate through each boundary
                        for i in range(level_variances[region][avg_num_comm].shape[1]-1):
                                if max(level_variances[region][avg_num_comm][1,i], level_variances[region][avg_num_comm][1,i+1]) > variance_threshold:

                                        level_consensus_zeroed[region][avg_num_comm][level_variances[region][avg_num_comm][0, i]:level_variances[region][avg_num_comm][0, i+1]] = 0
                        #print 'consensus post thresholding'
                        #print level_consensus_zeroed[region][avg_num_comm]
                        if hybrid == 'none':
                                # Plot the un-zeroed out consensus strip
                                if eval(settings["plots"]):
                                        fname = 'Consensus_pre_threshold_' + filename_qualifier + '_hierarchy_level_' + str(avg_num_comm)
                                        plot_consensus_partition(level_consensus[region][avg_num_comm], fname + '.png',  OUTDIR + '/variance_thresholding/consensus/prethreshold', my_cmap, vmax=max(level_consensus[region][avg_num_comm]))
                                        plt.clf()
                                        plt.close()

                                        filename3 = 'variance_barplot_%s_level_%3f.png' %(filename_qualifier, float(avg_num_comm))
                                        plt.figure(figsize = (6,2))
                                        strip_dir = OUTDIR + '/MMCP/variance_barplots/'
                                        if not os.path.isdir(strip_dir): os.makedirs(strip_dir)
                                        for i in range(1,len(level_variances[region][avg_num_comm][1])):
                                                width = 6
                                                pos = level_variances[region][avg_num_comm][0][i] - width/2
                                                plt.axvline(x = level_variances[region][avg_num_comm][0][i], ls = 'dashed', lw = '1', color = '0.75')
                                                plt.bar(pos, level_variances[region][avg_num_comm][1][i],  width = width, color = '0.25')
                                        plt.ylim((0,70))
                                        plt.xlim((0, num_nodes))
                                        plt.savefig(os.path.join(strip_dir, filename3), dpi = 900,  bbox_inches='tight')
                                        plt.clf()
                                        plt.close()

                                # save a .txt file of community assignment vectors for springforce diagrams of thresholded community assignments
                                fname = consensus_filtered_filename(region, filename_qualifier, avg_num_comm)
                                dir = OUTDIR + '/MMCP/consensus_filtered_files/'
                                if not os.path.isdir(dir): os.makedirs(dir)
                                np.savetxt(os.path.join(dir, fname), level_consensus_zeroed[region][avg_num_comm])



        if hybrid =='none':

                filename_qualifier = cell_type + '_' + tag

                create_results_output_genomewide_remove_trash(level_consensus, level_variances_presizethreshold, bin_location, filename_qualifier, settings, counts_as_arrays, 'none',tag)


        else:
                print('script does not handle hybrids!')


if __name__ == '__main__':
        main()
