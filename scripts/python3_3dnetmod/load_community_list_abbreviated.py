

def load_community_list(results_path):
        input = open(results_path, 'r')
        communities = [] # this will be a list of dicts

        for line in input:
                if line.startswith('#'):
                        continue
                else:
                        community = {}
                        community['chr'] = line.strip().split('\t')[0]
                        community['start'] = line.strip().split('\t')[1]
                        community['end'] = line.strip().split('\t')[2]
                        communities.append(community)
        input.close()

        return communities





def load_community_list_from_community_dict(community_dict):
        community_list = []
        for region in community_dict:
                for community in community_dict[region]:
                        community_list.append(community)



        return community_list








def load_community_list_centered_bin(results_path, bin_size):
        input = open(results_path, 'r')
        communities = []
        for line in input:
                if line.startswith('#'):
                        continue
                else:
                        community = {}
                        community['region'] = line.strip().split('\t')[0]
                        community['level'] = line.strip().split('\t')[1]
                        comm_name = line.strip().split('\t')[2]
                        community['number'] = line.strip().split('\t')[2]
                        community['chr'] = line.strip().split('\t')[3]
                        community['start'] = int(line.strip().split('\t')[4]) + bin_size/2 # start is defined by midpoint of first bin
                        community['end'] = int(line.strip().split('\t')[5]) - bin_size/2 # end is defined by midpoint of last bin
                        communities.append(community)
        input.close()

        return communities
