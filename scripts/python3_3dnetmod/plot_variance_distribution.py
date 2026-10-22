import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt



def plot_variance_distribution(variance_distribution, variance_threshold, output, color = '0.87'):
        mpl.rcdefaults()
        results, edges = np.histogram(variance_distribution, bins = 50, density = True)
        binWidth = edges[1] - edges[0]
        plt.bar(edges[:-1], results*binWidth, binWidth, color = color)
        plt.axvline(x = variance_threshold, color = 'red')
        plt.savefig(output, dpi = 300, bbox_inches = 'tight')
        plt.clf()
        plt.close()


def plot_variance_distribution_zoom(variance_distribution, variance_threshold, output, xmax, color = '0.87'):
        mpl.rcdefaults()
        #results, edges = np.histogram(variance_distribution, bins = 50, range = (0,50), density=True)
        variance_distribution = np.asarray(variance_distribution)
        to_plot = variance_distribution[variance_distribution<xmax]
        to_plot.tolist()
        if len(to_plot) > 0:
                #results, edges = np.histogram(to_plot, bins = 50, density = True)
                #binWidth = 0.3
                binWidth = (xmax*0.6)/15.0 # so that different zooms have the same bin width
                results, edges = np.histogram(to_plot, bins= np.arange(min(to_plot), max(to_plot) + binWidth, binWidth), density = True)
                #binWidth = edges[1] - edges[0]
                plt.bar(edges[:-1], results*binWidth, binWidth, color = color)
                plt.axvline(x = variance_threshold, color = 'red')
                plt.xlim((0,xmax))
                plt.savefig(output, dpi = 300, bbox_inches = 'tight')
                plt.clf()
                plt.close()

def plot_variance_distribution_zoomout(variance_distribution, variance_threshold, output, xmax, color = '0.87'):
        mpl.rcdefaults()
        #results, edges = np.histogram(variance_distribution, bins = 50, range = (0,50), density=True)
        variance_distribution = np.asarray(variance_distribution)
        to_plot = variance_distribution[variance_distribution<xmax]
        to_plot.tolist()
        if len(to_plot)>0:
                #results, edges = np.histogram(to_plot, bins = 50, density = True)
                #binWidth = 5
                binWidth = (xmax*0.3)/15.0
                results, edges = np.histogram(to_plot, bins= np.arange(min(to_plot), max(to_plot) + binWidth, binWidth), density = True)
                #binWidth = edges[1] - edges[0]
                plt.bar(edges[:-1], results*binWidth, binWidth, color = color)
                plt.axvline(x = variance_threshold, color = 'red')
                plt.xlim((0,xmax))
                plt.savefig(output, dpi = 300, bbox_inches = 'tight')
                plt.clf()
                plt.close()
