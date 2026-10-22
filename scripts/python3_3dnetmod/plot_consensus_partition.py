"""
05/21/16 Heidi Norton
This module plots consensus partitions and saves a .csv file of community assignments.

Does not clear figure before consensus. Requires debugging
"""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np



def plot_consensus_partition(consensus, figname, directory, my_cmap, vmax=None):

        length = consensus.size
        to_view = np.tile(consensus, (length,1))

        if vmax is None:
                fig1 = pyplot.imshow(to_view, extent = [0, to_view.shape[1], 0, to_view.shape[1]], interpolation = 'none', aspect = 0.03, cmap = my_cmap, vmin=0)
        else:
                fig1 = pyplot.imshow(to_view, extent = [0, to_view.shape[1], 0, to_view.shape[1]], interpolation = 'none', aspect = 0.03, cmap = my_cmap, vmin=0, vmax=vmax)

        # Create mean boundary lines
        # Locations where the community value changes (can't use np.unique because 0 appears many times)
        boundaries = np.where(consensus[1:] != consensus[:-1])[0]+1 #0 selects the numpy array out of a nested array; +1 to obtain 1st "new element" instead of last "old element"
        for boundary in boundaries:
                pyplot.axvline(boundary, linewidth = 1, color = 'k')

        fig1.axes.get_yaxis().set_visible(False)
        fig1.axes.get_xaxis().set_visible(False)
        dir = directory
        if not os.path.isdir(dir): os.makedirs(dir)
        pyplot.savefig(os.path.join(dir,figname), dpi =300, bbox_inches = 'tight')
