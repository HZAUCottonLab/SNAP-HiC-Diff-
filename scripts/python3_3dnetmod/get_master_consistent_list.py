

def get_master_consistent_list(communities1, communities2):
        master_consistent_list = []
        master_consistent_names = []
        for community in communities1:
                community_id = community['chr'] + '_' + community['start'] + '_' + community['end']
                if community_id not in master_consistent_names:
                        master_consistent_names.append(community_id)
                        master_consistent_list.append(community)

        for community in communities2:
                community_id = community['chr'] + '_' + community['start'] + '_' + community['end']
                if community_id not in master_consistent_names:
                        master_consistent_names.append(community_id)
                        master_consistent_list.append(community)


        return master_consistent_list
