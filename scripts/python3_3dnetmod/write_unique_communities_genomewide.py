import glob
import os
import numpy as np
from create_bin_location_dict import create_bin_location_dict
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from plot_consensus_partition import plot_consensus_partition
from get_partition_colormap import get_colormap
from load_region_start_end_coord import load_region_start_end_coord
from read_settings import read_settings
from outdir_helper import OUTDIR

def write_unique_communities_genomewide(settings, cell_type, cell_type_input,hybrid,tag_MMCP,tag_preprocess):
    """
    This script reads in a community results file and writes only the unique community coordinates. In cases of identically matched coordinates, the community with the lowest total boundary variance is selected.
    """
    if hybrid == 'none':

        file_name = 'Communities_' + cell_type + '_' + tag_MMCP + '.txt'

        results_file = OUTDIR + '/MMCP/communities_pre_filters/' +  file_name
        print('in write_unique_communities_genomewide, results file to write: ')
        print(results_file)
        print('results_file:' , results_file)
    else:
        print('script does not currently handle hybrids')

        #file_name = 'Communities_' + settings['num_part'] + '_partitions_hybrid_' + str(hybrid) + '_' + cell_type + '_gammastep_' + settings['gamma_step'] + '_variance_threshold_two_sided_max_' + settings['variance_threshold'] + '_' + settings['type_var_threshold']  + '_hierarchy_' + settings['hierarchies'] + '.txt'
        #results_file = 'output/variance_thresholding/results_files/hybrids/consensus_communities_thresholded/' + file_name
    bin_location = create_bin_location_dict('input/' + cell_type_input + '_' + tag_preprocess + 'final.bed' )


    #region_overlap_dict = load_region_overlap_coord('input/' + cell_type_input + 'final.bed')




    my_cmap = get_colormap()

    input = open(results_file, 'r')
    new_filename = file_name
    #new_filename = 'TEST_' + file_name
    communities = [] #list dictionaries of all communities, irrespective of region
    community_coords_dict = {} # dictionary where key = community coordinates, value = another dict with keys 'community_list' and 'variance_list'
    all_communities = [] # list of all community coordinates. Used to assess whether a given community is unique
    all_variance = []
    count = 0
    hierarchy_levels = {} # dictionary where key = region, value = list of hierarchy levels
    regions = []

    i = 0
    for line in input:
        if line.startswith('#'):
            continue
            #Region Gamma   Community_Assignment    Chr     Start_Pos       Stop_Pos        Left_Variance   Right_Variance
        else:
            community = {}
            community['region'] = line.strip().split('\t')[0]
            region = line.strip().split('\t')[0]
            community['gamma'] = line.strip().split('\t')[1]
            community_name = line.strip().split('\t')[2]
            community['number'] = int(community_name.split('_')[1])
            community['chr'] = line.strip().split('\t')[3]
            community['start'] = int(line.strip().split('\t')[4])
            community['stop'] = int(line.strip().split('\t')[5])
            community['left_var'] = float(line.strip().split('\t')[6])
            community['right_var'] = float(line.strip().split('\t')[7])
            community['hierarchy'] = line.strip().split('\t')[1]
            communities.append(community)
            if community['region'] not in regions:
                regions.append(community['region'])
                hierarchy_levels[community['region']] = []
            if community['hierarchy'] not in hierarchy_levels[community['region']]:
                hierarchy_levels[community['region']].append(community['hierarchy'])
            community_coordinates = ((line.strip().split('\t')[3], int(line.strip().split('\t')[4]), int(line.strip().split('\t')[5])))


            if community_coordinates not in community_coords_dict:
                community_coords_dict[community_coordinates] = {}
                community_coords_dict[community_coordinates]['community_list'] = []
                community_coords_dict[community_coordinates]['variance_list'] = []
                community_coords_dict[community_coordinates]['community_list'].append(community)
                community_coords_dict[community_coordinates]['variance_list'].append(community['left_var'] + community['right_var'])





    if hybrid == 'none':
            dir = OUTDIR + '/MMCP/unique_communities/results_files/'



    unique_communities = {}


    if  not os.path.isdir(dir): os.makedirs(dir)
    trimmed_communities = open(os.path.join(dir,new_filename), 'w')

    print('trimmed_communities: ', trimmed_communities)
    header = '#Region\tGamma_Parameter\tCommunity_Assignment\tChr\tStart_Pos\tStop_Pos\tLeft_var\tRight_Var'
    print(header, file=trimmed_communities)

    prefinalcommunity = []
    for community in community_coords_dict:
        #print "community_coords_dict[community]['variance_list']:", community_coords_dict[community]['variance_list']
        idx = np.argmin(community_coords_dict[community]['variance_list'])
        #print 'min variance:', community_coords_dict[community]['variance_list'][idx]
        community_to_write = community_coords_dict[community]['community_list'][idx]
        temp = '%s\t%s\t%s\t%s\t%i\t%i\t%3e\t%3e' %(community_to_write['region'], community_to_write['gamma'], community_to_write['number'], community_to_write['chr'], community_to_write['start'], community_to_write['stop'], community_to_write['left_var'], community_to_write['right_var'])
        print(temp, file=trimmed_communities)
        prefinalcommunity.append({'region':community_to_write['region'],'gamma':community_to_write['gamma'],'number':community_to_write['number'],'chr':community_to_write['chr'],'start':community_to_write['start'],'stop':community_to_write['stop'],'left_var':community_to_write['left_var'],'right_var':community_to_write['right_var']})
    trimmed_communities.close()

    regions = []
    community_dict = {}
    for community in prefinalcommunity:
        region = community['region']
        if region not in regions:
            regions.append(region)
            community_dict[region] = []
        community_dict[region].append(community)
    return community_dict


def main():
    settings = read_settings('settings2.txt')
    chr = 'chr6'
    cell_type_input = str(chr) + '_' + settings['sample_1']
    cell_type = str(chr) + settings['sample_1']
    write_unique_communities_genomewide(settings, cell_type, cell_type_input,'none')


if __name__ == '__main__':
    main()
