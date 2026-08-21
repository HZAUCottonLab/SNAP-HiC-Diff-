import os
from outdir_helper import OUTDIR
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as pyplot




def plot_size_distributions_per_stratum(final_communities_sizes, cell_type, sizes, HSVM_tag, settings):
        dir = OUTDIR + '/HSVM/size_hist/'
        if not os.path.isdir(dir): os.makedirs(dir)
        for size in sizes:
                size_dist = final_communities_sizes[size]
                print('size_dist: ')
                print(size_dist)
                filename = str(size) + HSVM_tag + '.png'
                results, edges = np.histogram(size_dist, bins = 50, range = (0,3000000), density=True)
                binWidth = edges[1] - edges[0]
                #adjusted_edges = edges - binWidth/2
                pyplot.bar(edges[:-1], results*binWidth, binWidth, color = '0.95')
                pyplot.xlim([0,3000000])
                pyplot.savefig(os.path.join(dir,filename), dpi = 300, bbox_inches = 'tight')
                pyplot.clf()
                pyplot.close()
