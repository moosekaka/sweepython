# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 13:01:49 2015
Functions for vionlinplots and stat sig bars
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
import seaborn as sns
import statsmodels.sandbox.stats.multicomp as mp
import config
import pandas as pd
import scipy.stats as sp


def mkdataset(datavals, datakeys):
    """Makes a flattened dataset with datavals in the first column and
    dataKeys in the second column for passing to mltComp function.
    The dimensions of dataVals must be compatible with dataKeys
    """
    temp = {}
    for k, j in enumerate(datavals):
        x = [j for j in datavals[k]]
        temp.update({y: datakeys[k] for y in x})
    return temp


def printtab(data):
    """Print every line of data with a newline for significant diffs
    """
    for line in data:
        if line[2] == True:
            print line


def mltcomp(dataset):
    """Performs multicomparison testing with stat function (default ranksums)
    and multiple testing correction method (default Holms)
    """
    mod = mp.MultiComparison(
        np.array(dataset.keys()),
        np.array(dataset.values()))
    res = mod.allpairtest(
        sp.ranksums, method='hs')
    restab = [
        (i[0], i[1], i[5], res[1][2][h]) for h, i
        in enumerate(res[2])]
    return restab


def label_diff(axes, labels, data, text):
    """label the significant samples with some pretty with a bar and '**'

    Parameters
    ----------
    ax :
        handle to plot to be labelled
    labels :
        labels for the samples
    data :
        data from the output of mltComp
    text :
        the text above the bar, default '**'
    """
    temp = []
    X = axes.get_xticks()
    Y = axes.get_ylim()
    for i in data:
        if i[2] == True:
            temp.append(
                (labels.index(i[0]), labels.index(i[1])))

    for k, j in enumerate(temp):
        y = np.log10(20*k+30)*.45*Y[1]   # vert. positions for labels
        props = {
            'connectionstyle': 'arc3', 'arrowstyle': '|-|',
            'shrinkA': 2, 'shrinkB': 2, 'lw': 1.25,
            'edgecolor': 'm'}
        axes.annotate(
            text, xy=(X[j[0]]+.4, y), zorder=10,
            size='large', weight='bold')
        axes.annotate(
            '', xy=(X[j[0]], y), xytext=(X[j[1]], y),
            arrowprops=props, ha='center', alpha=.3)


def count_funccall():
    """return a global index counter
    """
    config.counter += 1
    return config.counter


def pltvln(grphtitle, savename, vartoplot,
           labels, restab, **kwargs):
    """Plot using seaborn violin plot and pretty tables of stat sig.

    Parameters
    ----------
    grphTitle :
        Title for violinplot
    saveName :
        Filename to save as
    vartoPlot :
        The data to plot after applying user-defined functions
    labels :
        Labels for the x -axis
    resTab :
        Output of mltComp which uses statsmodels multiple testing
    """
    count_funccall()
    print'\n%02d. %s' % (config.counter, savename)+"\n"+"="*79
    printtab(restab)
    plt.figure(savename, figsize=(11, 8.5))

    frames = [
        pd.Series(vartoplot[i], name=labels[i])
        for i in range(len(labels))]
    srs = pd.concat(frames, axis=1)
    with sns.plotting_context('talk', font_scale=1.4):
        ax1 = sns.violinplot(srs)

        sns.stripplot(data=srs,
                      jitter=.05,
                      ax=ax1,
                      size=2,
                      alpha=.4,
                      color="gray")

        plt.title(grphtitle)
        pltlims = kwargs.pop('pltLims', None)  # for customizing ylims
        if pltlims is not None:
            plt.ylim(pltlims)
        label_diff(ax1, labels, restab, '**')
        plt.savefig('%02d_%s.png' % (config.counter, savename))
        return ax1
