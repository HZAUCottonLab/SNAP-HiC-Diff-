import numpy as np


def get_plateaus(avg_num_comm_array, gammas, plateau_size):
        hierarchies = {}
        temp_list = []
        level_blocks = {}
        #print 'keys in communities_per_gamma: should be gammas'
        #for key in communities_per_gamma:
                #print key


        #previous_avg_num_comm = round(np.mean(communities_per_gamma[gammas[0]].max(axis = 1)))
        previous_avg_num_comm = avg_num_comm_array[0]
        for i in range(len(gammas)):
                gamma = gammas[i]
                #print 'gamma in get_hierarchies', gamma
                #avg_num_comm = round(np.mean(communities_per_gamma[gamma].max(axis=1)))
                avg_num_comm = avg_num_comm_array[i]
                #print 'avg_num_comm', avg_num_comm
                if avg_num_comm == previous_avg_num_comm:
                        temp_list.append(gamma)
                elif avg_num_comm > previous_avg_num_comm:
                        if len(temp_list)>=plateau_size:
                                #hierarchies[previous_avg_num_comm] = {}
                                hierarchies[previous_avg_num_comm] = temp_list
                                #for item in temp_list:
                                        #hierarchies[previous_avg_num_comm][item]=communities_per_gamma[item]
                        temp_list = []
                        temp_list.append(gamma)
                if i == len(gammas) - 1:
                        if len(temp_list) > 1:
                                hierarchies[avg_num_comm] = {}
                                #for item in temp_list:
                                        #hierarchies[avg_num_comm][item] = communities_per_gamma[item]
                                hierarchies[avg_num_comm] = temp_list

                previous_avg_num_comm = avg_num_comm

        """
        for avg_num_comm in hierarchies:
                if avg_num_comm == 1:
                        continue
                level_blocks[avg_num_comm] = np.concatenate([hierarchies[avg_num_comm][gamma] for gamma in sorted(hierarchies[avg_num_comm])])
        """
        print('hierarchies: ')
        print(hierarchies)

        return hierarchies
