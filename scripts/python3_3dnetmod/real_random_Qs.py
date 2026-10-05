import numpy as np
from calculate_modularity import calculate_modularity
from scipy.stats import t


def real_random_Qs(settings,A,NG,random_networks,random_networks_NG,gamma):
    seed = 3443363518
    real_Qs_dict = {}
    Q_real, Ci = calculate_modularity(A, NG, gamma, seed)
    avg_random_Qs_list = np.zeros(len(random_networks))

    for i in range(0,len(random_networks)):
         Q_random, _ = calculate_modularity(random_networks[i], random_networks_NG[i], gamma, seed)
         avg_random_Qs_list[i] = Q_random

    avg_random_Qs = np.mean(avg_random_Qs_list)
    std_Q= np.std(avg_random_Qs_list)
    t_bounds = t.interval(0.95, len(avg_random_Qs_list) - 1) #student ts
    conf_Q = t_bounds[1]*(std_Q/np.sqrt(len(avg_random_Qs_list))) #95% confidence
    num_communities = max(Ci)

    return Q_real,avg_random_Qs,conf_Q, num_communities
