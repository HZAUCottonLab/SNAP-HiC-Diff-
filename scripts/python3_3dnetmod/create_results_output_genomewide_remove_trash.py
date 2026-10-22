"""
05/23/16
Heidi Norton and Jennifer Phillips-Cremins

Writes a community results file of the format:
#Region\tGamma_Parameter\tCommunity_Assignment\tChr\tStart_Pos\tStop_Pos


Input:
method: method used to determine significant consensus partitions

Modified 6/12/16 to only write communities that have a number > 0

"""
import os
from outdir_helper import OUTDIR
import numpy as np
from flatten_counts_single_region_geometric import flatten_counts_single_region

def create_results_output_genomewide_remove_trash(communities_dict, level_variances, bin_location, corrected_counts_file_input, settings, counts_as_array, hybrid,tag):
        """
        Example communities_dict[region][gamma]: [1 1 1 2 2 2 3 3 4 4 5 5]
        Example boundaries, aka level_variances[region][gamma][0]: [0, 3, 6, 8, 10, 12]
        Note that the boundaries range from 0 to 12 although there are only 12 nodes in the sample network.


        """
        if hybrid == 'none':
                dir = OUTDIR + '/MMCP/communities_pre_filters/'
        else:
                dir = OUTDIR + '/MMCP/communities_pre_filters/hybrids/'
        if not os.path.isdir(dir): os.makedirs(dir)

        fname1 ='Unthresholded_communities_' + corrected_counts_file_input  + '.txt'

        community_output_unthresholded = open(dir + fname1, 'w')


        fname2  = 'Communities_' + corrected_counts_file_input + '.txt'



        community_output_thresholded = open(dir + fname2, 'w')
        print('community_output_thresholded in write communities remove trash: ')

        header = '#Region\tGamma\tCommunity_Assignment\tChr\tStart_Pos\tStop_Pos\tLeft_Variance\tRight_Variance'
        print(header, file=community_output_unthresholded)
        print(header, file=community_output_thresholded)

        for region in communities_dict:
                counts = counts_as_array[region]
                np.fill_diagonal(counts, np.max(counts))
                #counts.diagonal = max(counts)
                threshold = np.percentile(flatten_counts_single_region(counts, discard_nan=True),float(settings['pctile_threshold']))
                print('threshold: ', threshold)
                for gamma in communities_dict[region]:
                        for i in range(level_variances[region][gamma].shape[1]-1):
                                community_start = int(level_variances[region][gamma][0,i])
                                community_stop = int(level_variances[region][gamma][0, i+1] - 1) # community_start to community_stop includes the bin numbers contained within the community
                                comm_column = counts[int(community_start):int(community_stop+1), int(community_stop)] #
                                #print 'community size: ', community_stop - community_start + 1
                                #print 'comm_column: ', comm_column
                                #print 'threshold : ', threshold
                                #print 'len(comm_column[comm_column>threshold])', len(comm_column[comm_column>threshold])
                                #print 'float(len(comm_column[comm_column>threshold]))/float(len(comm_column))', float(len(comm_column[comm_column>threshold]))/float(len(comm_column))
                                pct_value = float(len(comm_column[comm_column>threshold]))/float(len(comm_column)) # pct value is the % of counts along the edge of the community that must be greater than a given threshold, defined by a pctile count within the region
                                print('pct_value: percent of counts along edge greater than threshold', pct_value)
                                community_number = communities_dict[region][gamma][community_start]
                                community_start_coord = bin_location[region][community_start][2]
                                community_stop_coord = bin_location[region][community_stop][3]
                                left_var = level_variances[region][gamma][1,i]
                                right_var = level_variances[region][gamma][1, i+1]
                                chr= bin_location[region][i][1]

                                temp = "%s\t%.2f\t%s%i\t%s\t%i\t%i\t%3e\t%3e" % (region, gamma, 'community_',int(community_number),chr,community_start_coord, community_stop_coord, left_var, right_var)
                                print(temp, file=community_output_unthresholded) # all communities, including trash communities and communities too small
                                if community_stop - community_start + 1 >=2: # ensure that all communities are at least 2 nodes
                                        if pct_value >= float(settings['pct_value']): # line to modify for pct
                                                print('=======')
                                                print('passed pct threshold: ', pct_value)
                                                print('community size: ', community_stop - community_start + 1)
                                                print('=======')
                                                print(temp, file=community_output_thresholded)
                                        else:
                                                print('========')
                                                print('did not pass pct threshold: ', pct_value)
                                                print('community size: ', community_stop - community_start + 1)
                                                print('=======')
                                                print(temp)
                                else:
                                        print('community too small!')
                                        print('community size: ', community_stop - community_start + 1)


        community_output_thresholded.close()
        community_output_unthresholded.close()
