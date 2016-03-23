# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:18:08 2015

@author: sweel
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import islice, izip_longest
from matplotlib.gridspec import GridSpec as gs


def estimated_autocorrelation(vec):
    """Estimate autocorrelation formula according to
    https://en.wikipedia.org/wiki/Autocorrelation#Estimation
    convolution method using np.correlate

    Parameters
    ----------
    x :
        array of points, from on edge
    """
    X = np.array(vec)
    n = len(X)
    variance = X.var()
    X = X-X.mean()
    ret = np.correlate(X, X, mode='full')[-n:]
    result = ret/(variance*(np.arange(n, 0, -1)))
    return (result, n)


def plotauto(autodata, savename):
    """Plot the autocorrelation coefficient curves for a population
    with various minimum edge lengths given by minL and maxL

    Parameters
    ----------
    autoData :
        Dictionary containing list of outputs from autoCorOut:

        autoData[media_type][cell_name] = autoCorOut(cell,edgeData)
    saveName :
        file name
    """

    minl = 10
    maxl = 45
    celltemp = {}
    colors = {
        'YPD': 'b',
        'YPE': 'g',
        'YPR': 'm',
        'YPL': 'r'}

    media = sorted(autodata.keys())
    gs1 = gridspec.GridSpec(2, 2)
    ax1 = plt.figure(savename, figsize=(11, 8.25))

    for idx, mem in enumerate(media):
        ax2 = plt.subplot(gs1[idx])
        plt.ylim(-.5, 1.)
        ax2.text(7.5, .5, mem)
        ax2.axhline(color='k', alpha=0.5)

        for a, b in enumerate(range(minl, maxl, 5)):
            automedia = autodata[mem].values()
            X = np.ravel(
                [el for lis in automedia for el in lis[0]
                 if len(el) > b])
            celltemp.setdefault(mem, [])
            temp = []
            for i in islice(izip_longest(*X, fillvalue=.0), 0, 15):
                temp.append(np.mean(i))
            celltemp[mem].append(temp)

            plt.plot(celltemp[mem][a],
                     alpha=.35 * np.log(a+2), color=colors[mem])

    plt.suptitle(savename +
                 "Pop. AutoCor Coef for various minimum "
                 "edge lengths of %4.2f - %4.2f $\mu$m in "
                 " increments of %4.3f $\mu$m" %
                 (minl*.055, maxl*.055, .055), fontsize=13)
    plt.savefig(savename+'.png')
    plt.show()
    return ax1


def plotcelledgeauto(labs, autodata, savename):
    """Plot the autocorrelation coefficient curves for each cell and each
    edge  various minimum edge lengths given by minL and maxL

    Parameters
    ----------
    labs:
        media type label
    autoData :
        Dictionary containing list of outputs from autoCorOut:
        autoData[media_type][cell_name] = autoCorOut(cell,edgeData)
    saveName :
        file name to save as
    """
    n = int(np.ceil(np.sqrt(len(autodata[labs].keys()))))
    grid = gs(n, n)
    plt.figure(figsize=(17, 11))

    for ind, name in enumerate(autodata[labs].keys()):
        ax1 = plt.subplot(grid[ind])
        X = autodata[labs][name]
        for i in X[0]:
            ax1.plot([j for j in i])
    plt.suptitle(labs)
    plt.savefig('%s_%s.png' % (labs, savename))
    plt.close()