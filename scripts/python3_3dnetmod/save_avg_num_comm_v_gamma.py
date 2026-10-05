import os
from outdir_helper import OUTDIR


def save_avg_num_comm_v_gamma(num_comm, gammas, settings, cell_type,tag):
        for region in gammas:
                filename = 'Avg_num_comm_v_gamma_' + cell_type + '_' + region  + '_' + tag
                dir = OUTDIR + '/GPS/gamma_dynamic_range/'
                dynamic_gammas = open(dir + filename, 'w')
                print('#gamma' + '\t' + 'avg_num_comm' + '\t', file=dynamic_gammas)
                for i in range(len(num_comm[region])):
                        gamma = gammas[region][i]
                        num = num_comm[region][i]
                        temp = str(gamma) + '\t' + str(num)
                        print(temp, file=dynamic_gammas)
                dynamic_gammas.close()
