from counts_to_array import counts_to_array
from outdir_helper import OUTDIR
from construct_nulls import construct_NG_null
from get_random_Qs import get_random_Qs
from real_random_Qs import real_random_Qs
from counts_to_array import counts_to_array
import numpy as np
import os
import glob
import sys

def determine_max_gamma_per_region(elements):
        settings = elements[0]
        chr = elements[1]
        region = elements[2]
        tags_preprocess = elements[3]
        print('tags_preprocess: ', tags_preprocess)
        tags_GPS = elements[4]
        cell_type = elements[5]
        print('cell_type: ', cell_type)
        counts = 'input/' + chr + '_' + cell_type + '_' + tags_preprocess + "finalpvalues.counts"
        bed = 'input/' + chr + '_' + cell_type + '_' + tags_preprocess + "final.bed"
        print('counts: ', counts)
        print('bed: ', bed)
        counts_as_arrays, _  = counts_to_array(counts,bed)
        inst = region

        random_networks_NG = {}
        real_Qs_dict = {}
        random_Qs_dict = {}
        random_Qconfs_dict = {}
        mean_diffs_dict = {}
        gammas = {}
        max_gamma_list = []


        random_networks = get_random_Qs(settings, counts_as_arrays, cell_type, region, tags_GPS, chr)
        for i in range(0,len(random_networks)):
                random_networks_NG[i] = []
                random_networks_NG[i] =  construct_NG_null(random_networks[i]) #find null matrix
        A = counts_as_arrays[region]
        np.fill_diagonal(A, 0) # zero out the diagonal of the network to match random rewiring scenario
        NG = construct_NG_null(A) # construct the Newman  Girvan null model
        real_Qs_dict[inst] = []
        random_Qs_dict[inst] = []
        random_Qconfs_dict[inst] = []
        mean_diffs_dict[inst] = []
        gammas[inst] = []

        gammas[inst] = []
        gamma = 0.0
        counter = 0
        real_Q = 1
        avgdiff_change = 10
        counterlist = []
        historylist = []
        size = 10
        particular_gamma = 0
        maxcounter = 1
        max_change = -1

        if 'rate' not in settings:
                settings['rate']  = 0.0048 #default rate

        if 'max_gamma' in list(settings.keys()):
                while gamma <= float(settings['max_gamma']):
                        real_Q, random_Qs, conf_Q, num_communities = real_random_Qs(settings,A,NG,random_networks,random_networks_NG,gamma)
                        diff = real_Q - random_Qs
                        counterlist.append(counter)
                        historylist.append(diff)
                        history = np.abs(historylist[-1*int(size):]) #capture last 10 diffs
                        if gamma >= 1.5: # assume that convergence has to occur atleat past 1.5 gamma
                                if diff < 0:
                                        max_gamma_list.append(gamma) #just pick gamma as real and random cross paths after 1.5 (to stabilize must go past 1.5)
                                        particular_gamma = gamma
                                        break
                                maxtest = max(historylist) #retrieve largest positive (careful not to use abs)
                                if maxtest > max_change:
                                        max_change = maxtest
                                        for i, j in enumerate(historylist):
                                                 if j == maxtest:
                                                        maxcounter = i
                                diffwmax = max_change - historylist[-1] #diff between largest encountered and current beyond 1.5 gamma
                                adjustmentpergamma = (float(counter) - float(maxcounter))
                                if int(counter) > int(maxcounter):
                                        avgdiff_change = float(diffwmax)/adjustmentpergamma
                                        if avgdiff_change <= float(settings['rate']):
                                                max_gamma_list.append(gamma)
                                                particular_gamma = gamma
                                if np.sum(np.abs(np.diff(history)))/float((size-1)) < float(settings['rate']): # find mean of change in difference over history, converge if below 0.005 change in difference per 0.1 gamma (one less diff than size of window (size -1))
                                        avgdiff_change = float(settings['rate'])/10 #ensure convergence
                                        max_gamma_list.append(gamma)
                                        particular_gamma = gamma
                        if real_Q < -0.1: # stop if real modularity below -0.1
                                max_gamma_list.append(gamma)
                                particular_gamma = gamma
                                break
                        counter = counter + 1
                        real_Qs_dict[inst].append(real_Q)
                        random_Qs_dict[inst].append(random_Qs)
                        random_Qconfs_dict[inst].append(conf_Q)
                        gammas[inst].append(gamma)
                        mean_diffs_dict[inst].append(diff)
                        gamma = gamma + 0.1

        else:
                while (((avgdiff_change > (float(settings['rate'])) or avgdiff_change < 0))): #cutoff at steadystate difference or random crosses paths with  real, removed requirement	 to be below 0.05
                        real_Q, random_Qs, conf_Q, num_communities = real_random_Qs(settings,A,NG,random_networks,random_networks_NG,gamma)
                        diff = real_Q - random_Qs
                        counterlist.append(counter)
                        historylist.append(diff)
                        history = np.abs(historylist[-1*int(size):]) #capture last 10 diffs
                        if gamma >= 1.5: # assume that convergence has to occur atleat past 1.5 gamma
                                if diff < 0:
                                        max_gamma_list.append(gamma) #just pick gamma as real and random cross paths after 1.5 (to stabilize must go past 1.5)
                                        particular_gamma = gamma
                                        break
                                maxtest = max(historylist) #retrieve largest positive (careful not to use abs)
                                if maxtest > max_change:
                                        max_change = maxtest
                                        for i, j in enumerate(historylist):
                                                 if j == maxtest:
                                                        maxcounter = i
                                diffwmax = max_change - historylist[-1] #diff between largest encountered and current beyond 1.5 gamma
                                adjustmentpergamma = (float(counter) - float(maxcounter))
                                if int(counter) > int(maxcounter):
                                        avgdiff_change = float(diffwmax)/adjustmentpergamma
                                        if avgdiff_change <= float(settings['rate']):
                                                max_gamma_list.append(gamma)
                                                particular_gamma = gamma
                                if np.sum(np.abs(np.diff(history)))/float((size-1)) < float(settings['rate']): # find mean of change in difference over history, converge if below 0.005 change in difference per 0.1 gamma (one less diff than size of window (size -1))
                                        avgdiff_change = float(settings['rate'])/10 #ensure convergence
                                        max_gamma_list.append(gamma)
                                        particular_gamma = gamma
                        if real_Q < -0.1: # stop if real modularity below -0.1
                                max_gamma_list.append(gamma)
                                particular_gamma = gamma
                                break
                        counter = counter + 1
                        real_Qs_dict[inst].append(real_Q)
                        random_Qs_dict[inst].append(random_Qs)
                        random_Qconfs_dict[inst].append(conf_Q)
                        gammas[inst].append(gamma)
                        mean_diffs_dict[inst].append(diff)
                        gamma = gamma + 0.1

        # Write a text file of real modularity and random modularity across gamma values
        dir = OUTDIR + '/GPS/gamma_dynamic_range/'
        fname = 'Real_random_modularity_v_gamma_' + chr + '_' + cell_type + '_' + region + '_' + tags_GPS + '.txt'
        out_file = open(dir + fname, 'w')
        header = '#Gamma\tReal_Q\tRandom_Q\tMean_diffs'
        print(header, file=out_file)
        for i in range(0,len(real_Qs_dict[region])):
                temp = '%f\t%f\t%f\t%f' %(gammas[region][i],real_Qs_dict[region][i],random_Qs_dict[region][i],mean_diffs_dict[region][i])
                print(temp, file=out_file)
        out_file.close()


        print("max_gamma_list ", max_gamma_list)
        if len(max_gamma_list) == 0:
                if gamma > 0:
                        max_gamma_list.append(gamma - 0.1)
                else:
                        max_gamma_list.append(gamma)
        elif particular_gamma == 0: #still overlooked
                max_gamma_list.append(gamma - 0.1)

        max_gamma_dict = {}
        max_gamma_dict[cell_type] = max(max_gamma_list)

        return max_gamma_dict
