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


def mkDataset(dataVals, dataKeys):
    """Makes a flattened dataset with datavals in the first column and
    dataKeys in the second column for passing to mltComp function.
    The dimensions of dataVals must be compatible with dataKeys
    """
    temp = {}
    for k, j in enumerate(dataVals):
        x = [j for j in dataVals[k]]
        temp.update({y: dataKeys[k] for y in x})
    return temp


def printTab(data):
    """Print every line of data with a newline for significant diffs
    """
    for line in data:
        if line[2] == True:
            print(line)


def mltComp(dataSet):
    """Performs multicomparison testing with stat function (default ranksums)
    and multiple testing correction method (default Holms)
    """
    mod = mp.MultiComparison(
        np.array(dataSet.keys()),
        np.array(dataSet.values()))
    res = mod.allpairtest(
        sp.ranksums, method='hs')
    resTab = [
        (i[0], i[1], i[5], res[1][2][h]) for h, i
        in enumerate(res[2])]
    return resTab


def label_diff(ax, labels, data, text):
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
    X = ax.get_xticks()
    Y = ax.get_ylim()
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
        ax.annotate(
            text, xy=(X[j[0]]+.4, y), zorder=10,
            size='large', weight='bold')
        ax.annotate(
            '', xy=(X[j[0]], y), xytext=(X[j[1]], y),
            arrowprops=props, ha='center', alpha=.3)


def countFuncCall():
    config.counter += 1
    return config.counter


def pltVln(grphTitle, saveName, vartoPlot,
           labels, resTab, *args, **kwargs):
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
    countFuncCall()
    print('\n%02d. %s' % (config.counter, saveName)+"\n"+"="*79)
    printTab(resTab)
    plt.figure(saveName)
    ax1 = sns.violinplot(vartoPlot,
                         cut=kwargs.pop('cut', .5),  # changed default value
                         names=labels, **kwargs)
    plt.title(grphTitle)
    pltLims = kwargs.pop('pltLims', None)  # for customizing ylims
    if pltLims is not None:
        plt.ylim(pltLims)
    label_diff(ax1, labels, resTab, '**')
    plt.savefig('%02d_%s.png' % (config.counter, saveName))
    return(ax1)
