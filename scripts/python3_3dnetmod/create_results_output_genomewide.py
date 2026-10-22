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

def create_results_output_genomewide(communities_dict, level_variances, variance_threshold, bin_location, corrected_counts_file_input, num_part, settings, hybrid = 'none'):
        """
        Example communities_dict[region][gamma]: [1 1 1 2 2 2 3 3 4 4 5 5]
        Example boundaries, aka level_variances[region][gamma][0]: [0, 3, 6, 8, 10, 12]
        Note that the boundaries range from 0 to 12 although there are only 12 nodes in the sample network.


        """
        if hybrid == 'none':
                dir1 = OUTDIR + '/variance_thresholding/results_files/consensus_communities_unthresholded/'
        else:
                dir1 = OUTDIR + '/variance_thresholding/results_files/hybrids/consensus_communities_unthresholded/'
        if not os.path.isdir(dir1): os.makedirs(dir1)

        fname1 ='Unthresholded_communities_' + str(num_part) + '_partitions_' + corrected_counts_file_input + '.txt'
        community_output_unthresholded = open(os.path.join(dir1,fname1), 'w')

        if hybrid == 'none':
                dir2 = OUTDIR + '/variance_thresholding/results_files/consensus_communities_thresholded/'
        else:
                dir2 = OUTDIR + '/variance_thresholding/results_files/hybrids/consensus_communities_thresholded/'
        if not os.path.isdir(dir2): os.makedirs(dir2)

        fname2 = 'Communities_' + str(num_part) + '_partitions_' + corrected_counts_file_input + '.txt'
        community_output_thresholded = open(os.path.join(dir2,fname2), 'w')

        header = '#Region\tGamma\tCommunity_Assignment\tChr\tStart_Pos\tStop_Pos\tLeft_Variance\tRight_Variance'
        print(header, file=community_output_unthresholded)
        print(header, file=community_output_thresholded)

        for region in communities_dict:
                for gamma in communities_dict[region]:

                        for i in range(level_variances[region][gamma].shape[1]-1):
                                community_start = level_variances[region][gamma][0,i]
                                community_stop = level_variances[region][gamma][0, i+1] - 1 # community_start to community_stop includes the bin numbers contained within the community
                                community_number = communities_dict[region][gamma][community_start]
                                community_start_coord = bin_location[region][community_start][2]
                                community_stop_coord = bin_location[region][community_stop][3]
                                left_var = level_variances[region][gamma][1,i]
                                right_var = level_variances[region][gamma][1, i+1]
                                #print 'left_var: ', left_var
                                #print 'right_var: ', right_var
                                #print 'community_start_coord: ', community_start_coord
                                #print 'community_stop_coord: ', community_stop_coord
                                chr= bin_location[region][i][1]
                                #print 'chr: ', chr
                                if community_stop - community_start + 1 > int(settings['size_threshold']):
                                        temp = "%s\t%.2f\t%s%i\t%s\t%i\t%i\t%3e\t%3e" % (region, gamma, 'community_',int(community_number),chr,community_start_coord, community_stop_coord, left_var, right_var)
                                        print(temp, file=community_output_unthresholded)
                                        if max(level_variances[region][gamma][1,i], level_variances[region][gamma][1,i+1]) <= variance_threshold: # if community is within variance threshold, write to results file
                                                print(temp, file=community_output_thresholded)
                                #else:
                                        #print 'found small community, size = ', community_stop - community_start + 1

        community_output_thresholded.close()
        community_output_unthresholded.close()
