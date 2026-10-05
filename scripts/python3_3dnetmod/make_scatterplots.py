import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import os

def make_scatterplots(x, y, filename, dir):
        matplotlib.rcdefaults()
        if not os.path.isdir(dir): os.makedirs(dir)
        fig = pyplot.figure(figsize = (6, 4))
        ax1 = pyplot.subplot(111)
        l=ax1.plot(x, y, '.', color = 'g', markersize = 5)
        l = ax1.fill_between(x, y)
        ax1.set_xlim(np.min(x), np.max(x))
        ax1.set_ylim(0, np.max(y) + 5)
        ax1.grid('off')
        ax1.tick_params(axis = 'both', which = 'major', labelsize = 6)

        # Change the fill into a blueish color
        l.set_facecolors([[.5,.8,.5,0]])

        # change th edge color and thicknesss
        l.set_edgecolors([[0, .5, 0, .3]])
        l.set_linewidths([1])
        pyplot.savefig(os.path.join(dir,filename), dpi = 300)
        pyplot.clf()


def make_scatterplots_zoom(x, y, filename, dir, xrange):
        matplotlib.rcdefaults()
        if not os.path.isdir(dir): os.makedirs(dir)
        fig = pyplot.figure(figsize = (6, 4))
        ax1 = pyplot.subplot(111)
        l=ax1.plot(x, y, '.', color = 'g', markersize = 10)
        l = ax1.fill_between(x, y)
        #ax1.set_xlim(np.min(x), np.max(x))
        ax1.set_xlim(xrange)
        #ax1.set_ylim(0, np.max(y) + 5)
        ax1.set_ylim(0, 12)
        ax1.grid('off')
        ax1.tick_params(axis = 'both', which = 'major', labelsize = 6)

        # Change the fill into a blueish color
        l.set_facecolors([[.5,.8,.5,0]])

        # change th edge color and thicknesss
        l.set_edgecolors([[0, .5, 0, .3]])
        l.set_linewidths([1])
        pyplot.savefig(os.path.join(dir,filename), dpi = 300)
        pyplot.clf()


def make_scatterplots_no_line(x, y, filename, dir):

        if not os.path.isdir(dir): os.makedirs(dir)
        fig = pyplot.figure(figsize = (6, 4))
        ax1 = pyplot.subplot(111)
        l=ax1.plot(x, y, '.', color = 'g', markersize = 5)
        #l = ax1.fill_between(x, y)
        ax1.set_xlim(np.min(x), np.max(x))
        ax1.set_ylim(0, np.max(y) + 5)
        ax1.grid('on')
        ax1.tick_params(axis = 'both', which = 'major', labelsize = 6)

        # Change the fill into a blueish color
        #l.set_facecolors([[.5,.8,.5,0]])

        # change th edge color and thicknesss
        #l.set_edgecolors([[0, .5, 0, .3]])
        #l.set_linewidths([1])
        pyplot.savefig(os.path.join(dir,filename), dpi = 300)
        pyplot.clf()


def make_scatterplots_threshold(x, y, threshold, filename, dir):

        if not os.path.isdir(dir): os.makedirs(dir)
        fig = pyplot.figure(figsize = (6, 4))
        ax1 = pyplot.subplot(111)
        l=ax1.plot(x, y, '.', color = 'g', markersize = 5)
        l = ax1.fill_between(x, y)
        ax1.set_xlim(np.min(x), np.max(x))
        ax1.set_ylim(0, np.max(y) + 5)
        ax1.grid('on')
        ax1.tick_params(axis = 'both', which = 'major', labelsize = 6)

        # Change the fill into a blueish color
        l.set_facecolors([[.5,.8,.5,0]])

        # change th edge color and thicknesss
        l.set_edgecolors([[0, .5, 0, .3]])
        l.set_linewidths([1])
        pyplot.axhline(y=threshold, color = 'b')
        pyplot.savefig(os.path.join(dir,filename), dpi = 300)
        pyplot.clf()
