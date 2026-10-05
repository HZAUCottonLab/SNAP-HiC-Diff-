from check_intersect import check_intersect
import numpy as np

def load_community_dict(results_path):
        input = open(results_path, 'r')
        community_dict = {}
        regions = []
        for line in input:
                if line.startswith('#'):
                        continue
                else:
                        community = {}
                        region = line.strip().split('\t')[0]
                        if region not in regions:
                                regions.append(region)
                                community_dict[region] = []
                        community['region'] = line.strip().split('\t')[0]
                        community['gamma'] = line.strip().split('\t')[1]
                        comm_name = line.strip().split('\t')[2]
                        community['number'] = line.strip().split('\t')[2]
                        community['chr'] = line.strip().split('\t')[3]
                        community['start'] = line.strip().split('\t')[4]
                        community['end'] = line.strip().split('\t')[5]
                        community['left_var'] = float(line.strip().split('\t')[6])
                        community['right_var'] = float(line.strip().split('\t')[7])
                        community_dict[region].append(community)

        input.close()

        return community_dict


def load_community_dict_centered_bin(results_path, bin_size):
        input = open(results_path, 'r')
        community_dict = {}
        regions = []
        for line in input:
                if line.startswith('#'):
                        continue
                else:
                        community = {}
                        region = line.strip().split('\t')[0]
                        if region not in regions:
                                regions.append(region)
                                community_dict[region] = []
                        community['region'] = line.strip().split('\t')[0]
                        community['gamma'] = line.strip().split('\t')[1]
                        comm_name = line.strip().split('\t')[2]
                        community['number'] = line.strip().split('\t')[2]
                        community['chr'] = line.strip().split('\t')[3]
                        community['start'] = line.strip().split('\t')[4] + bin_size/2
                        community['end'] = line.strip().split('\t')[5] - bin_size/2
                        community_dict[region].append(community)

        input.close()

        return community_dict
