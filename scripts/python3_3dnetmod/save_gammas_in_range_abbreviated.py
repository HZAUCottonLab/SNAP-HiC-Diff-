import os
from outdir_helper import OUTDIR


def save_gammas_in_range(settings, cell_type, new_gammas,tag):
        dir = OUTDIR + '/GPS/gamma_dynamic_range/'
        if not os.path.isdir(dir): os.makedirs(dir)
        for region in new_gammas:
                filename = 'gamma_dynamic_range_' + cell_type + '_' + region  + '_' + tag
                final_gammas = open(dir + filename, 'w')
                print('#gamma', file=final_gammas)
                for i in range(len(new_gammas[region])):
                        print(new_gammas[region][i], file=final_gammas)
                final_gammas.close()
