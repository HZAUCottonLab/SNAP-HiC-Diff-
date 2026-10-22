import os
from outdir_helper import OUTDIR
import matplotlib.pyplot as plt

def plot_partition(Ci, cmap, filename, vmin=0, vmax=None, subdir='', vlines=None):
        """ Plot a partition block Ci with a given cmap, and optional vmin, vmax kwargs for cmap

        Inputs
                Ci: ndarray, n_part by n_nodes
                cmap: matplotlib colormap
                filename: str, filename to save partition block as
                vmin: int/float, min normalization value for cmap
                vmax: int/float, max normalization value for cmap
                subdir: str, subfolder after output/communities_across_numpart/ to save file in
                vlines: iterable ints, x-coordinates of consensus boundary vertical lines

        Outputs
                none
        """
        plt.clf()

        # Optionally specify the upper range of colormap
        if vmax is None:
                plt.imshow(Ci, interpolation = 'none', cmap=cmap, vmin=vmin)
                plt.axis('off')
        else:
                plt.imshow(Ci, interpolation = 'none', cmap=cmap, vmin=vmin, vmax=vmax)
                plt.axis('off')

        # plot consensus lines if they are inputted
        if not vlines is None:
                for xcoord in vlines:
                        plt.axvline(x=xcoord, linewidth='2', color='k')

        # Create directory
        dir = OUTDIR + '/MMCP/communities_across_numpart/' + subdir
        if not os.path.isdir(dir): os.makedirs(dir)

        # Save figure
        plt.savefig(os.path.join(dir,filename), dpi = 900, bbox_inches = 'tight')
        plt.clf()
        plt.close()
