import numpy as np



def get_variance_distribution_and_threshold(final_communities, sizes, variance_settings):

        variance_distribution = {}
        variance_thresholds = {}
        thresholded_communities = {}
        for i in range(len(sizes)):
                size = sizes[i]
                var_thresh = variance_settings[i]
                thresholded_communities[size] = []
                if len(final_communities[size]) == 0:
                        #variance_distribution[size] = np.zeros(0)
                        variance_distribution[size] = []
                        variance_thresholds[size] = 0
                        #thresholded_communities[size] = np.zeros(0)
                        thresholded_communities[size] = []

                else:
                        #variance_distribution[size] = np.zeros(len(final_communities[size]))
                        variance_distribution[size] = []
                        print('variance_distribution:', variance_distribution)
                        for j in range(0, len(final_communities[size])):
                                #variance_distribution[size][j] = max(float(final_communities[i]['left_var']),float(final_communities[size][i]['right_var']))
                                variance_distribution[size].append(max(float(final_communities[size][j]['left_var']),float(final_communities[size][j]['right_var'])))
                        if var_thresh > 0:
                                variance_threshold = np.percentile(variance_distribution[size], var_thresh)
                        else:
                                variance_threshold = var_thresh
                        for k in range(0, len(final_communities[size])):
                                if max(float(final_communities[size][k]['left_var']), float(final_communities[size][k]['right_var'])) <= variance_threshold:
                                        thresholded_communities[size].append(final_communities[size][k])



                        variance_thresholds[size] = variance_threshold
        return variance_distribution, variance_thresholds, thresholded_communities
