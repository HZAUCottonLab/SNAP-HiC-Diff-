import numpy as np
import math
from check_intersect_equality import check_intersect






def merge_communities(total_thresholded_communities, bin_size, d):
        comparison_communities = total_thresholded_communities[:]
        comparison_options = np.arange(0, len(total_thresholded_communities))

        tracking_dict = {}



        while len(comparison_options) > 1:
                community = comparison_communities[comparison_options[0]]
                                # take first element of comparison options, index into total thresholded communities
                community_name = community['chr'] + '_' + str(community['start']) + '_' + str(community['end'])
                tracking_dict[community_name] = {}
                tracking_dict[community_name]['chr'] =  community['chr']
                tracking_dict[community_name]['communities'] = []
                tracking_dict[community_name]['start_coord'] = {}
                tracking_dict[community_name]['end_coord'] = {}
                tracking_dict[community_name]['start_coord']['chrom'] = community['chr']
                tracking_dict[community_name]['start_coord']['start'] = int(community['start']) - d
                tracking_dict[community_name]['start_coord']['end'] = int(community['start']) + d
                tracking_dict[community_name]['start_coord']['list'] = []
                tracking_dict[community_name]['end_coord']['start'] = int(community['end']) - d
                tracking_dict[community_name]['end_coord']['end'] = int(community['end']) + d
                tracking_dict[community_name]['end_coord']['chrom'] = community['chr']
                tracking_dict[community_name]['community_name'] = community_name
                tracking_dict[community_name]['end_coord']['list'] = []
                to_delete = []
                for option in comparison_options:
                        #community2 = comparison_communities[comparison_options[option]]
                        community2 = comparison_communities[option]
                        if community2['chr'] == community['chr']:
                                boundary1 = {}
                                boundary1['chrom'] = community2['chr']
                                boundary1['start'] = int(community2['start'])
                                boundary1['end'] = int(community2['start'])
                                boundary2 = {}
                                boundary2['chrom'] = community2['chr']
                                boundary2['start'] = int(community2['end'])
                                boundary2['end'] = int(community2['end'])
                                if check_intersect(tracking_dict[community_name]['start_coord'], boundary1) and check_intersect(tracking_dict[community_name]['end_coord'], boundary2):
                                        tracking_dict[community_name]['communities'].append(community2)
                                        tracking_dict[community_name]['start_coord']['list'].append(boundary1['start'])
                                        tracking_dict[community_name]['end_coord']['list'].append(boundary2['end'])
                                        if boundary1['start'] == tracking_dict[community_name]['start_coord']['start']:
                                                #print 'expanding start boundary start coordinate. communities involved: '
                                                #print community
                                                #print community2
                                                tracking_dict[community_name]['start_coord']['start'] = tracking_dict[community_name]['start_coord']['start'] - d
                                        elif boundary1['start'] == tracking_dict[community_name]['start_coord']['end']:
                                                tracking_dict[community_name]['start_coord']['end'] = tracking_dict[community_name]['start_coord']['end'] + d
                                                #print 'expanding start boundary end coordinate. communities involved: '
                                                #print community
                                                #print community2
                                        else:
                                                print('not expanding start coordinates. communities involved: ')
                                                #print community
                                                #print community2

                                        if boundary2['start'] == tracking_dict[community_name]['end_coord']['start']:
                                                print('expanding end boundary strat coordinate. communities involved: ')
                                                print(community)
                                                print(community2)
                                                tracking_dict[community_name]['end_coord']['start'] = tracking_dict[community_name]['end_coord']['start'] - d
                                        elif boundary2['start'] == tracking_dict[community_name]['end_coord']['end']:
                                                tracking_dict[community_name]['end_coord']['end'] = tracking_dict[community_name]['end_coord']['end'] + d
                                                #print 'expanding end boundary end coordinate. communities involved: '
                                                #print community
                                                #print community2
                                        else:
                                                print('not expanding end coordinates. communities involved: ')
                                                #print community
                                                #print community2

                                        #to_delete.append(option)
                                        to_delete.append(comparison_options.tolist().index(option))
                #print 'deleting because they matched: ', to_delete
                comparison_options = np.delete(comparison_options, to_delete)

        tracking_dict_list = []
        for key in tracking_dict:
                tracking_dict_list.append(tracking_dict[key])


        comparison_options = np.arange(0, len(tracking_dict_list))


        tracking_dict2 = {}
        while len(comparison_options) > 1:
                #print 'tracking_dict_list[comparison_options[0]]'
                community = tracking_dict_list[comparison_options[0]]
                community_name = community['community_name']
                tracking_dict2[community_name] = {}
                tracking_dict2[community_name]['chr'] = community['chr']
                tracking_dict2[community_name]['communities'] = community['communities']
                tracking_dict2[community_name]['start_coord'] = community['start_coord']
                tracking_dict2[community_name]['end_coord'] = community['end_coord']
                tracking_dict2[community_name]['community_name'] = community_name
                to_delete = []
                for option in comparison_options:
                        community2 = tracking_dict_list[option]
                        #print 'community2: '
                        #print community2
                        if community2['chr'] == community['chr']:
                                boundary1 = {}
                                boundary1['chrom'] = community2['chr']
                                boundary1['start'] = min(community2['start_coord']['list']) # might want to make this community2['start_coord']['start']. this would extend the search window
                                boundary1['end'] = max(community2['start_coord']['list'])
                                boundary2 = {}
                                boundary2['chrom'] = community2['chr']
                                boundary2['start'] = min(community2['end_coord']['list'])
                                boundary2['end'] = max(community2['end_coord']['list'])
                                if check_intersect(tracking_dict2[community_name]['start_coord'], boundary1) and check_intersect(tracking_dict2[community_name]['end_coord'], boundary2):
                                        tracking_dict2[community_name]['communities'].append(community2)
                                        tracking_dict2[community_name]['start_coord']['list'].append(boundary1['start'])
                                        tracking_dict2[community_name]['end_coord']['list'].append(boundary2['end'])
                                        """
                                        if boundary1['start'] == tracking_dict2[community_name]['start_coord']['start']:
                                                #print 'expanding start boundary start coordinate. communities involved: '
                                                #print community
                                                #print community2
                                                tracking_dict2[community_name]['start_coord']['start'] = tracking_dict2[community_name]['start_coord']['start'] - d
                                        elif boundary1['start'] == tracking_dict2[community_name]['start_coord']['end']:
                                                tracking_dict2[community_name]['start_coord']['end'] = tracking_dict2[community_name]['start_coord']['end'] + d
                                                #print 'expanding start boundary end coordinate. communities involved: '
                                                #print community
                                                #print community2
                                        else:
                                                print 'not expanding start coordinates. communities involved: '
                                                #print community
                                                #print community2

                                        if boundary2['start'] == tracking_dict2[community_name]['end_coord']['start']:
                                                #print 'expanding end boundary strat coordinate. communities involved: '
                                                #print community
                                                #print community2
                                                tracking_dict2[community_name]['end_coord']['start'] = tracking_dict2[community_name]['end_coord']['start'] - d
                                        elif boundary2['start'] == tracking_dict2[community_name]['end_coord']['end']:
                                                tracking_dict2[community_name]['end_coord']['end'] = tracking_dict2[community_name]['end_coord']['end'] + d
                                                #print 'expanding end boundary end coordinate. communities involved: '
                                                #print community
                                                #print community2
                                        else:
                                                print 'not expanding end coordinates. communities involved: '
                                                #print community
                                                #print community2
                                        """

                                        to_delete.append(comparison_options.tolist().index(option))
                #print 'deleting because they matched: ', to_delete
                comparison_options = np.delete(comparison_options, to_delete)

                #print 'comparison options after deleting: '
                #print comparison_options



        merged_communities = []


        #print 'tracking dict communities: '

        for community in tracking_dict2:
                #print community
                merged_community = {}
                merged_community['chr'] = tracking_dict2[community]['chr']
                merged_community['start1'] = int(min(tracking_dict2[community]['start_coord']['list']))
                merged_community['start2'] = int(max(tracking_dict2[community]['start_coord']['list']))
                temp_start = int(np.mean(tracking_dict2[community]['start_coord']['list']))
                merged_community['stop1'] = int(min(tracking_dict2[community]['end_coord']['list']))
                merged_community['stop2'] = int(max(tracking_dict2[community]['end_coord']['list']))
                temp_stop = int(np.mean(tracking_dict2[community]['end_coord']['list']))

                merged_communities.append(merged_community)





        return merged_communities
