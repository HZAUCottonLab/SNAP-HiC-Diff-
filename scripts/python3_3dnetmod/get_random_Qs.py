import numpy as np
from outdir_helper import OUTDIR
from construct_nulls import construct_NG_null
from null_model_und_sign import null_model_und_sign
from calculate_modularity import calculate_modularity
from scipy.stats import t
import os
import glob
import csv
import sys

def get_random_Qs(settings, counts_as_arrays, cell_type, region,tag, chr):
    seed = 3443363518
    random_Qs_dict={}
    random_Qconfs_dict={}
    random_networks = {}
    random_networks_NG = {}
    dir = OUTDIR + '/GPS/random_networks/'
    if not os.path.isdir(dir): os.makedirs(dir)
    existing_random_networks = glob.glob(OUTDIR + '/GPS/random_networks/' + cell_type + '_' + chr + str(region) + '_random_network_*_' + tag + '.csv')
    print('existing_random_networks:')
    print(existing_random_networks)
    random_Qs_dict[region] = []
    random_Qconfs_dict[region] = []
    if 'num_rand' not in settings:
        settings['num_rand'] = 1  #default setting
    if 'switch' not in settings:
        settings['switch'] = 'True' #default setting

    for i in range(int(settings['num_rand'])):
        fname = dir + cell_type + '_' + chr + region +'_random_network_' + str(i) + '_' + tag + '.csv'
        if fname in existing_random_networks:
            print('already have random network', i)
            print('filename:')
            print(fname)
            R=np.array(list(csv.reader(open(fname, "rb"),delimiter=','))).astype('float')
        else:

            print('counts_as_arrays.keys(): ')
            print(list(counts_as_arrays.keys()))
            passin = counts_as_arrays[region]
            R,rpos,rneg = null_model_und_sign(passin,settings['switch'],1,wei_freq=0.1)  #calculate null, True = swap bins with sign change,wei_freq = freq of random # generation
            np.savetxt(fname, R, delimiter=',', fmt='%10.5f')
        random_networks[i] = R

    return random_networks
