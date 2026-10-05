#from load_community_list import load_community_list
from load_community_list_abbreviated import load_community_list
import sys
from check_intersect import check_intersect

from read_settings import read_settings
from get_master_consistent_list import get_master_consistent_list


def find_master_consistent_communities(element):
        chr = element[0]
        sample1_list = element[1]
        sample2_list = element[2]
        size_overlap = element[3]

        community = sample1_list[0]
        if 'start1' in list(community.keys()):
                (chr,communities1_overlap, communities2_overlap) = find_consistent_communities_merged(element)
        elif 'start' in list(community.keys()):
                (chr,communities1_overlap, communities2_overlap) = find_consistent_communities_unmerged(element)
        else:
                print('community does not have start or start1 in keys!')
        master_communities = get_master_consistent_list(communities1_overlap, communities2_overlap)

        return master_communities
def find_consistent_communities(element):
        chr = element[0]
        sample1_list = element[1]
        sample2_list = element[2]
        size_overlap = element[3]

        community = sample1_list[0]
        if 'start1' in list(community.keys()):
                (chr,communities1_overlap, communities2_overlap) = find_consistent_communities_merged(element)
        elif 'start' in list(community.keys()):
                (chr,communities1_overlap, communities2_overlap) = find_consistent_communities_unmerged(element)
        else:
                print('community does not have start or start1 in keys!')
                return None
        return (chr,communities1_overlap, communities2_overlap)


def find_consistent_communities_unmerged(element):

    chr = element[0]
    sample1_list = element[1]
    sample2_list = element[2]
    size_overlap = element[3]

    s1_intersect_s2 = 0
    s2_intersect_s1 = 0


    community1_names = []
    community2_names = []

    communities1_overlap = []

    communities2_overlap = []
    for community1 in sample1_list:
        community_count = 0
        a1 = {}
        a1['chrom'] = community1['chr']
        a1['start'] = int(community1['start']) - size_overlap
        a1['end'] = int(community1['start']) + size_overlap
        b1 = {}
        b1['chrom'] = community1['chr']
        b1['start'] = int(community1['end']) - size_overlap
        b1['end'] = int(community1['end']) + size_overlap

        name = community1['chr'] + '_' + community1['start'] + '_' + community1['end']
        community1_names.append(name)

        for community2 in sample2_list:
            a2 = {}
            a2['chrom'] = community2['chr']
            a2['start'] = int(community2['start']) - 1
            a2['end'] = int(community2['start']) + 1
            b2 = {}
            b2['chrom'] = community2['chr']
            b2['start'] = int(community2['end']) - 1
            b2['end'] = int(community2['end']) + 1



            if check_intersect(a1, a2) and check_intersect(b1,b2):
                community_count +=1
        if community_count > 0:
            s1_intersect_s2 +=1
            communities1_overlap.append(community1)


    for community2 in sample2_list:
        community_count = 0
        a1 = {}
        a1['chrom'] = community2['chr']
        a1['start'] = int(community2['start']) - size_overlap
        a1['end'] = int(community2['start']) + size_overlap
        b1 = {}
        b1['chrom'] = community2['chr']
        b1['start'] = int(community2['end']) - size_overlap
        b1['end'] = int(community2['end']) + size_overlap

        name = community2['chr'] + '_' + community2['start'] + '_' + community2['end']
        community2_names.append(name)

        for community1 in sample1_list:
            a2 = {}
            a2['chrom'] = community1['chr']
            a2['start'] = int(community1['start']) - 1
            a2['end'] = int(community1['start']) + 1
            b2 = {}
            b2['chrom'] = community1['chr']
            b2['start'] = int(community1['end']) - 1
            b2['end'] = int(community1['end']) + 1

            if check_intersect(a1, a2) and check_intersect(b1,b2):
                community_count +=1
        if community_count > 0:
            s2_intersect_s1 +=1
            communities2_overlap.append(community2)



    c1_v_2_exact_match = 0
    for name in community1_names:
        if name in community2_names:
            c1_v_2_exact_match +=1

    print('number of communities sample1: ', len(sample1_list))
    print('number of consistent communities: ', len(communities1_overlap))

    return (chr,communities1_overlap, communities2_overlap)



def find_consistent_communities_merged(element):

    chr = element[0]
    sample1_list = element[1]
    sample2_list = element[2]
    size_overlap = element[3]

    s1_intersect_s2 = 0
    s2_intersect_s1 = 0


    community1_names = []
    community2_names = []

    communities1_overlap = []

    communities2_overlap = []

    for community1 in sample1_list:
        community_count = 0
        a1 = {}
        a1['chrom'] = community1['chr']
        a1['start'] = int(community1['start1']) - size_overlap
        a1['end'] = int(community1['start2']) + size_overlap
        b1 = {}
        b1['chrom'] = community1['chr']
        b1['start'] = int(community1['end1']) - size_overlap
        b1['end'] = int(community1['end1']) + size_overlap

        name = community1['chr'] + '_' + community1['start1'] + '_' + community1['end2']
        community1_names.append(name)

        for community2 in sample2_list:
            a2 = {}
            a2['chrom'] = community2['chr']
            a2['start'] = int(community2['start1']) - 1
            a2['end'] = int(community2['start2']) + 1
            b2 = {}
            b2['chrom'] = community2['chr']
            b2['start'] = int(community2['end1']) - 1
            b2['end'] = int(community2['end2']) + 1



            if check_intersect(a1, a2) and check_intersect(b1,b2):
                community_count +=1
        if community_count > 0:
            s1_intersect_s2 +=1
            communities1_overlap.append(community1)




    for community2 in sample2_list:
        community_count = 0
        a1 = {}
        a1['chrom'] = community2['chr']
        a1['start'] = int(community2['start1']) - size_overlap
        a1['end'] = int(community2['start2']) + size_overlap
        b1 = {}
        b1['chrom'] = community2['chr']
        b1['start'] = int(community2['end1']) - size_overlap
        b1['end'] = int(community2['end2']) + size_overlap

        name = community2['chr'] + '_' + community2['start1'] + '_' + community2['end2']
        community2_names.append(name)

        for community1 in sample1_list:
            a2 = {}
            a2['chrom'] = community1['chr']
            a2['start'] = int(community1['start1']) - 1
            a2['end'] = int(community1['start2']) + 1
            b2 = {}
            b2['chrom'] = community1['chr']
            b2['start'] = int(community1['end1']) - 1
            b2['end'] = int(community1['end2']) + 1

            if check_intersect(a1, a2) and check_intersect(b1,b2):
                community_count +=1
        if community_count > 0:
            s2_intersect_s1 +=1
            communities2_overlap.append(community2)



    c1_v_2_exact_match = 0
    for name in community1_names:
        if name in community2_names:
            c1_v_2_exact_match +=1

    print('number of communities sample1: ', len(sample1_list))
    print('number of consistent communities: ', len(communities1_overlap))

    return (chr,communities1_overlap, communities2_overlap)



if __name__ == '__main__':
    main()
