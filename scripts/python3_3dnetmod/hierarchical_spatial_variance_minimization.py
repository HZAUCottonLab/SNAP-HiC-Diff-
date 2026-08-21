import sys
import glob
import os
from outdir_helper import OUTDIR
import copy
import numpy as np
from merge_communities import merge_communities
from get_variance_distribution_and_threshold import get_variance_distribution_and_threshold
from plot_variance_distribution import *
from plot_size_distributions_per_stratum import plot_size_distributions_per_stratum
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from load_region_start_end_coord import load_region_start_end_coord
from create_bin_location_dict import create_bin_location_dict
from find_left_right_boundary_from_partition import find_left_right_boundary_from_partition
import math
from scipy.stats import percentileofscore




#def HSVM(settings,community_dict,cell_type_input,cell_type,tag_HSVM,tag_preprocess):
def HSVM(settings,community_dict,sample, tag_HSVM, tag_preprocess, chr):
        """
        This module takes an input dictionary of communities (key is region, value is list of community dictionaries) and identifies communities that do not have either boundary at the edge of the region (network). The communities are then thresholded based on the variance specified in the settings file. The resulting communities that pass the size and variance threshold are written to a results file.
        """
        # cell_type needs to have major chromosome in the name
        bad_communities = []
        sizes = []
        variance_settings = []
        variance_order = []
        final_communities = {}
        final_communities_sizes = {} # for troubleshooting
        variance_distributions = {}

        for size_inst in settings:
                if size_inst[:6]=='size_s':
                        sizes.append(int(settings[size_inst]))
                        final_communities[int(settings[size_inst])] = []
                        final_communities_sizes[int(settings[size_inst])] = []
        sizes = sorted(sizes)

        variance_tag = ''
        size_tag = ''
        for var_inst in settings:
                if var_inst[:10] == 'var_thresh':
                        variance_settings.append(int(settings[var_inst]))
                        variance_order.append(var_inst)
        variance_settings = [x for (y,x) in sorted(zip(variance_order,variance_settings))]

        #variance_settings = sorted(variance_settings)

        print('sizes: ', sizes)
        print('variance_settings: ', variance_settings)
        for item in variance_settings:
                variance_tag = variance_tag + str(item) + '_'


        for item in sizes:
                size = str(int(item)/1000) + 'kb'
                size_tag = size_tag + size + '_'



        if len(sizes) != len(variance_settings):
                print('must have the same number of size stratifications as variance thresholds! Exiting!')
                return None

        for subchrom in community_dict:
                #have to work this out
                #print 'looking for bed file: ', 'input/' + cell_type_input + "_" + tag_preprocess + "final.bed"
                bed_file = 'input/' + subchrom + '_' + sample + '_' + tag_preprocess + "final.bed"
                #input/chr18.14_CKKO1Neuronalchr18_100_50_True_Falsefinal.bed
                print('looking for bed_file: ', bed_file)
                (region_idx, bin_size)  = load_region_start_end_coord(bed_file)
                bin_location = create_bin_location_dict(bed_file)
                for region in community_dict[subchrom]:
                        start_boundary = region_idx[region]['start']
                        end_boundary = region_idx[region]['end']
                        for community in community_dict[subchrom][region]:
                                if int(community['start']) - int(start_boundary)<=bin_size*4:
                                        bad_communities.append(community)
                                        print('lost community to edge')
                                        print('start community: ', community['start'])
                                        print('start region: ', start_boundary)
                                elif int(end_boundary) - int(community['end']) <=bin_size*4:
                                        print('lost community to edge')
                                        print('end community: ', community['end'])
                                        print('end region: ', end_boundary)
                                        bad_communities.append(community)
                                elif abs(int(community['end']) - int(community['start'])) <= bin_size*int(settings['size_threshold']):
                                        print('lost community to size threshold')
                                        bad_communities.append(community)
                                else:
                                        previous_size = 0
                                        for i in range(len(sizes)):
                                                if abs(int(community['end']) - int(community['start'])) > previous_size and abs(int(community['end']) - int(community['start'])) <= sizes[i]:
                                                        final_communities[sizes[i]].append(community)
                                                        final_communities_sizes[sizes[i]].append(abs(int(community['end']) - int(community['start'])))
                                                previous_size = sizes[i]


        if settings['plots'] == 'True':
                plot_size_distributions_per_stratum(final_communities_sizes, chr + sample, sizes, tag_HSVM, settings)
        variance_distributions, var_thresholds, thresholded_communities = get_variance_distribution_and_threshold(final_communities, sizes, variance_settings)


        dir = OUTDIR + '/HSVM/variance_thresholds/'
        if not os.path.isdir(dir): os.makedirs(dir)
        variance_filename = 'Variance_thresholds_' + chr + sample + '_' + tag_HSVM + '.txt' #HOW DOES CHROMOSOME FIT IN HERE?
        variance_out = open(dir + variance_filename, 'w')
        header = '#Min_size' + '\t' + 'max_size' + '\t' + 'threshold_AUC' + '\t' + 'threshold_value'
        print(header, file=variance_out)
        previous_size = 0
        for i in range(len(sizes)):
                size = sizes[i]
                var_thresh = variance_settings[i]
                if var_thresh == 0:
                        var_thresh  = percentileofscore(variance_distributions[size], 0, kind = 'weak')
                else:
                        var_thresh = variance_settings[i]
                print(str(previous_size) + '\t' + str(size) + '\t' + str(var_thresh) + '\t' + str(var_thresholds[size]), file=variance_out)
                previous_size = size

        variance_out.close()

        dir = OUTDIR + '/HSVM/variance_distributions/'
        if not os.path.isdir(dir): os.makedirs(dir)

        previous_size = 0
        colors=cm.rainbow(np.linspace(0,1,len(sizes)))
        for i in range(len(sizes)):
                size = sizes[i]
                var_thresh = variance_settings[i]
                if settings['plots'] == 'True':
                        histname = str(previous_size) + '_to_' + str(size) + chr + sample + '_'  + tag_HSVM + '.png' # HOW DOES CHROMOSOME FIT IN HERE?
                        plot_variance_distribution(variance_distributions[size], var_thresholds[size], dir + histname, color = colors[i])
                        plot_variance_distribution_zoom(variance_distributions[size], var_thresholds[size], dir + '60_zoom' + histname, 60, color = colors[i])
                        plot_variance_distribution_zoom(variance_distributions[size], var_thresholds[size], dir + '10_zoom' + histname, 10, color = colors[i])
                        plot_variance_distribution_zoom(variance_distributions[size], var_thresholds[size], dir + '250_zoom' + histname, 250, color = colors[i])
                        plot_variance_distribution_zoom(variance_distributions[size], var_thresholds[size], dir + '4_zoom' + histname, 4, color = colors[i])
                        plot_variance_distribution_zoom(variance_distributions[size], var_thresholds[size], dir + '5_zoom' + histname, 5, color = colors[i])
                        plot_variance_distribution_zoom(variance_distributions[size], var_thresholds[size], dir + '6_zoom' + histname, 6, color = colors[i])
                        plot_variance_distribution_zoom(variance_distributions[size], var_thresholds[size], dir + '7_zoom' + histname, 7, color = colors[i])
                        plot_variance_distribution_zoom(variance_distributions[size], var_thresholds[size], dir + '150_zoom' + histname, 150, color = colors[i])
                        #plot_variance_distribution_zoom(variance_distributions[size], var_thresholds[size], dir + '200_zoom' + histname, 200, color = colors[i])


                previous_size = size


        total_thresholded_communities = []
        for i in range(len(sizes)):
                size = sizes[i]
                total_thresholded_communities += thresholded_communities[size]


        print('number of total_thresholded_communities: ', len(total_thresholded_communities))


        dir = OUTDIR + '/HSVM/variance_thresholded_communities/unmerged/'
        if not os.path.isdir(dir): os.makedirs(dir)
        final_output = 'Communities_'+ chr + sample + '_' + tag_HSVM + '.txt' # HOW DOES CHROMOSOME FIT IN HERE?
        print('final output: ', dir + final_output)
        community_output_thresholded = open(dir + final_output, 'w')
        header = '#Region\tGamma\tCommunity_Assignment\tChr\tStart_Pos\tStop_Pos\tLeft_Variance\tRight_Variance'
        print(header, file=community_output_thresholded)
        for j in range(0,len(total_thresholded_communities)):
                print('total_thresholded_communities[j]', total_thresholded_communities[j])
                temp = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(total_thresholded_communities[j]['region'], total_thresholded_communities[j]['gamma'], total_thresholded_communities[j]['number'], total_thresholded_communities[j]['chr'], total_thresholded_communities[j]['start'], total_thresholded_communities[j]['end'], total_thresholded_communities[j]['left_var'], total_thresholded_communities[j]['right_var'])
                print(temp, file=community_output_thresholded)  #final list of communities with bounded removed and below threshold
        community_output_thresholded.close()


        total_merged_communities = merge_communities(total_thresholded_communities, bin_size, int(settings['boundary_buffer']))
        dir = OUTDIR + '/HSVM/variance_thresholded_communities/merged/'
        if not os.path.isdir(dir): os.makedirs(dir)
        final_output = 'Merged_' + 'Communities_' + chr + sample  + '_' + tag_HSVM + '.txt' # HOW DOES CHROMOSOME FIT IN HERE?
        community_output_thresholded_merged = open(dir + final_output, 'w')
        header = '#Chr\tUpstream_boundary_start\tUpstream_boundary_stop\tDownstream_boundary_start\tDownstream_boundary_stop'
        print(header, file=community_output_thresholded_merged)
        for community in total_merged_communities:
                temp = community['chr'] + "\t" + str(community['start1']) + "\t" + str(community['start2']) + "\t" + str(community['stop1']) + "\t" + str(community['stop2'])
                print(temp, file=community_output_thresholded_merged)
        community_output_thresholded_merged.close()
