import glob
from outdir_helper import OUTDIR
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as pyplot
from counts_to_array import counts_to_array
#from determine_max_gamma_genomewide import determine_max_gamma
from construct_nulls import construct_NG_null
from make_scatterplots import make_scatterplots
from make_scatterplots import make_scatterplots_zoom
from real_Qs_avg import real_Qs_avg
from get_plateaus import get_plateaus
from save_avg_num_comm_v_gamma import save_avg_num_comm_v_gamma
from save_gammas_in_range_abbreviated import save_gammas_in_range
from convert_tag_genomewide import convert_tag_genomewide

"""
Heidi Norton, Daniel Emerson

Determine the useful range of the structural resolution parameter, gamma, by computing modularity for a real and randomly rewired networks across a range of structural resolution parameter values.

Final useful gamma is determined as the gamma at which the modularity of the real network converges with modularity of the random network.

Assumes convergence has to occur after 1.5 gamma.

User defined settings should be specified in settings.txt

"""

def GPS(chr, settings, sample, badregions,tag_GPS,tag_MMCP,tag_preprocess, max_gamma):

    if 'single_gamma' in list(settings.keys()):
        return 0
    avg_num_comm_files = glob.glob(OUTDIR + '/GPS/gamma_dynamic_range/Avg_num_comm_v_gamma_' + str(chr) + sample + '_*' + tag_GPS)

    tag_GPS_test = convert_tag_genomewide(tag_GPS, settings)
    tag_MMCP_test = convert_tag_genomewide(tag_MMCP, settings)
    avg_num_comm_files2 = glob.glob(OUTDIR + '/GPS/gamma_dynamic_range/Avg_num_comm_v_gamma_' + str(chr) + sample + '_*' + tag_GPS_test)

    if len(avg_num_comm_files) > 0:
        print('already found gamma plateau sweep files')
        unique_communities = glob.glob(OUTDIR + '/MMCP/unique_communities/results_files/Communities_' + str(chr) + sample + '_' + tag_MMCP + '.txt')
        print('unique_communities: ', unique_communities)
        if len(unique_communities) > 0:  #if unique communities already found in addition to same Avg_num_comm_v_gamma
            print("exiting MMCP: already found community output with identical tag and sample name")
            return 1
        else: # just skip max gamma determination with Avg_num_comm_v_gamma present (still computes unique community)
            return 0
    elif len(avg_num_comm_files2) > 0:
        print('already found gamma plateau sweep files')
        unique_communities = glob.glob(OUTDIR + '/MMCP/unique_communities/results_files/Communities_' + str(chr) + sample + '_' + tag_MMCP_test + '.txt')
        print('unique_communities: ', unique_communities)
        if len(unique_communities) > 0:
            print("exiting MMCP: already found community output with identical tag and sample name")
            return 1
        else:
            return 0
    else:

        print('looked for: ', OUTDIR + '/GPS/gamma_dynamic_range/Avg_num_comm_v_gamma_' + chr + sample + '_' + tag_GPS)
        print('and: ')
        print(OUTDIR + '/GPS/gamma_dynamic_range/Avg_num_comm_v_gamma_' + chr + sample + '_' + tag_GPS_test)

    if 'num_rand' not in settings:
        settings['num_rand'] = 1

    num_rand = int(settings['num_rand'])
    cell_type_input = str(chr) + '_' + sample
    counts = 'input/' +  cell_type_input + '_' + tag_preprocess + "finalpvalues.counts"
    bed = 'input/' +  cell_type_input + '_' + tag_preprocess + "final.bed"
    cell_type = str(chr) + sample
    counts_as_arrays,bin_location = counts_to_array(counts, bed)
    for region in badregions:
        if region[1] == chr:
             if region[0] in counts_as_arrays:
                  del counts_as_arrays[region[0]]
    gammas_finer = {}
    num_comm = {}
    real_Q_finer = {}
    random_Qs_finer = {}
    plateau_dict = {}

    if 'gamma_step' not in settings:
        settings['gamma_step'] = 0.01 #set default to 0.01 if missing

    gamma_step = float(settings['gamma_step'])

    try:
        x = 1/int(1/gamma_step)
    except ZeroDivisionError:
         print("Cannot exceed step size of 1 gamma")

    if max_gamma == None:
         return None
    else:

         for region in counts_as_arrays:
              real_Q = 1
              gamma = 0.0
              gammas_finer[region] = []
              num_comm[region] = []
              real_Q_finer[region] = []
              random_Qs_finer[region] = []

              # Perform a finer sampling of gamma to get plateaus in avg # community v gamma
              gamma = 0
              mod = []
              #while gamma <= max(gammas[region]):
              A = counts_as_arrays[region]
              np.fill_diagonal(A, 0) # zero out the diagonal of the network to match random rewiring scenario
              NG = construct_NG_null(A) # construct the Newman  Girvan null model
              while gamma <= max_gamma:
                  real_Q,num_communities  = real_Qs_avg(settings,A,NG,gamma)
                  gammas_finer[region].append(gamma)
                  num_comm[region].append(num_communities)
                  real_Q_finer[region].append(real_Q)
                  gamma = gamma + float(settings['gamma_step'])


              identifier = cell_type + '_' + chr + region + '_' + tag_GPS
              # Make plots
              if settings['plots'] == 'True':

                  filename = 'Avg_num_comm_v_gamma_' + identifier + '.png'
                  make_scatterplots(gammas_finer[region], num_comm[region], filename, OUTDIR + '/GPS/scatterplots/')

                  filename = 'Avg_num_comm_v_gamma_zoom1' + identifier + '.png'
                  make_scatterplots_zoom(gammas_finer[region], num_comm[region], filename, OUTDIR + '/GPS/scatterplots/', (2,3))


                  filename = 'Avg_num_comm_v_gamma_zoom2' + identifier + '.png'
                  make_scatterplots_zoom(gammas_finer[region], num_comm[region], filename, OUTDIR + '/GPS/scatterplots/', (3,4))

                  filename = 'Avg_num_comm_v_Q_' + identifier + '.png'
                  make_scatterplots(real_Q_finer[region], num_comm[region], filename, OUTDIR + '/GPS/scatterplots/')

                  filename = 'Q_v_Avg_num_comm_' + identifier + '.png'
                  make_scatterplots(num_comm[region], real_Q_finer[region], filename, OUTDIR + '/GPS/scatterplots/')

                  filename = 'Q_v_gamma_' + identifier + '.png'
                  make_scatterplots(gammas_finer[region], real_Q_finer[region], filename, OUTDIR + '/GPS/scatterplots/')



         save_avg_num_comm_v_gamma(num_comm, gammas_finer, settings, cell_type,tag_GPS)

         plateau_gammas = {}
         for region in counts_as_arrays:
              plateau_gammas[region] = []
              plateau = get_plateaus(num_comm[region], gammas_finer[region], int(settings['plateau']))
              for avg_num_comm in plateau:
                   if avg_num_comm > 1:   # only want to keep gammas that lead to at least 2 communities
                        plateau_gammas[region].append(np.mean(plateau[avg_num_comm]))

              print('final gammas: ')
              print(plateau_gammas[region])
              plateau_gammas[region].sort()
         save_gammas_in_range(settings, cell_type, plateau_gammas,tag_GPS)
         return 0




if __name__ == '__main__':
        main()
