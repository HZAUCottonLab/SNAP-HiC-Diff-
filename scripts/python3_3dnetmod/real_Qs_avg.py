import numpy as np
from calculate_modularity import calculate_modularity
from assign_communities_gen_louvain import assign_sorted_communities_gen_louvain
from scipy.stats import t
from read_seeds import read_seeds


def real_Qs_avg(settings,A,NG,gamma):

    if 'gamma_part' not in settings:
        settings['gamma_part'] = 1
    num = int(settings['gamma_part'])
    real_Qs_dict = {}
    #Q_real, Ci = calculate_modularity(A, NG, gamma, seed)
    seed = 3443363518
    seed_count = 0
    if 'seed_file' in list(settings.keys()):
        gen_louvain_seeds = read_seeds(settings['seed_file'], num)
    else:
        gen_louvain_seeds = read_seeds('gen_louvain_100_seeds.txt', num)
    [Ci, Q, seed_count] = assign_sorted_communities_gen_louvain(A, NG, num, gamma, gen_louvain_seeds, seed_count)
    Q_real = np.mean(Q)
    print('Q_real: ', Q_real)


    num_communities = round(np.mean(Ci.max(axis=1)))

    return Q_real, num_communities
