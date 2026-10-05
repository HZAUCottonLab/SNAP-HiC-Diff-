from matplotlib import cm
import numpy as np
import matplotlib as mpl
mpl.use('Agg')

def get_colormap():
        # Custom jet colormap with 0 = black
        my_cmap = cm.get_cmap('jet', 100)
        my_cmap_values = my_cmap(np.arange(100))
        my_cmap_values = np.insert(my_cmap_values, 0, [0,0,0,0], axis=0)
        my_cmap = mpl.colors.LinearSegmentedColormap.from_list("my_cmap", my_cmap_values)
        return my_cmap
