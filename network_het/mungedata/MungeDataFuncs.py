# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 18:10:24 2015
Functions for munging dataset and plotting in sns and pandas
@author: sweel
"""
from tvtk.api import tvtk
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as sp
import seaborn as sns
import statsmodels.sandbox.stats.multicomp as mp


def adjust_axes(fighandle, numticks):
    """adjust number of ticks
    """
    ax_loc = fighandle.axes.shape
    for loc in range(ax_loc[0]):
        x = fighandle.axes[ax_loc[0]-1, loc]
        xbounds = x.get_xbound()
        x.set_xticks(np.round(
            np.linspace(0., xbounds[1], numticks),
            2))

        y = fighandle.axes[loc, 0]
        ybounds = y.get_ybound()
        y.set_yticks(np.round(
            np.linspace(0., ybounds[1], numticks),
            2))


def bpts_inten(vtkdata, bptscoord, radinf=.3):
    """return list of branchpoints ptIDs in nodes with
    averaged intensity values within a radius of influence from data
    """
    big = {}
    polydata = tvtk.PolyData()  # vtk object as input for Locator
    polydata.points = vtkdata.points

    for point in bptscoord:
        result = tvtk.IdList()
        loc = tvtk.PointLocator(data_set=polydata)
        loc.build_locator()
        loc.find_points_within_radius(
            radinf, bptscoord[point], result)
        big[point] = result
    return big


def boxviol(datf, vals, group, **kwargs):
    """plot violin with boxplot bounds and stripplots
    """
    with sns.plotting_context('talk', font_scale=1.25):
        _, ax1 = plt.subplots(1, 1)

        #  VIOLINPLOT
        sns.violinplot(x=group,
                       ax=ax1,
                       y=vals,
                       data=datf)

#        xpos = [tic for tic in ax1.get_xticks()]
#
#        #  BOXPLOT
#        dic = datf.boxplot(vals,
#                           by=group,
#                           ax=ax1,
#                           showbox=False,
#                           showmeans=True,
#                           showfliers=False,
#                           whiskerprops={'linewidth': 0},
##                           medianprops={'linewidth': 3},
#                           capprops={'linewidth': 2,
#                                     'markersize': 1,
#                                     'color': '#555555'},
#                           positions=xpos,
#                           return_type='dict')

#        for dline in dic[dic.keys()[0]]['medians']:
#            dline.set_color('#FFFFFF')
#            x1, x2 = dline.get_xdata()
#            dline.set_xdata([x1 + .1, x2 - .1])

        ax1.set_title('')
        plt.suptitle('')
        #  STRIPPLOT
        sns.stripplot(x=group,
                      y=vals,
                      ax=ax1,
                      data=datf,
                      jitter=.05,
                      size=3,
                      alpha=.4)

        #  LABELS AND LIMS
#        pltlims = []
#        for cap in dic[dic.keys()[0]]['caps']:
#            pltlims.append(cap.get_ydata()[0])
#        lowerb = min(pltlims)
#        upperb = max(pltlims)
#        plt.ylim(max(0, (lowerb - .05 * lowerb)), 1.05 * upperb)
#        pltlims = kwargs.pop('pltLims', None)  # for customizing ylims
#        if pltlims is not None:
#            plt.ylim(pltlims)
#        else:
#            plt.ylim(0.)
        plt.show()
        return ax1


def multiple_test(dataset, col1, col2):
    """Performs multicomparison testing with stat function (default ranksums)
    and multiple testing correction method (default Holms)

    Parameters
    ----------
        col1:
            the data vals in group to be tested
        col2:
            the group labels for testing
    Returns
    -------
        res:
            statsmodel data object, use res[0] to get the simpletable class
    """
    mod = mp.MultiComparison(dataset[col1], dataset[col2])
    res = mod.allpairtest(
        sp.mannwhitneyu, method='hs')
    return res


def result_to_excel(res, sheetname, writer):
    """Write out result of multiple_test to excel worksheet
    Parameters
    ----------
        res:
            the data vals in group to be tested
        sheetname:
            the sheetname which the result will be inserted into
        writer:
            an instance of Pandas ExcelWriter
    """
    obj = res[0]
    txt = obj.title.split('\n')
    txt2 = [i for i in txt[2].split(',')]
    txt3 = [txt[1], txt2[0], txt2[1]]

    ftemp = obj.as_html()
    ftemp = pd.read_html(ftemp)
    gtemp = pd.DataFrame(ftemp[0])
    gtemp.columns = gtemp.ix[0, :]
    gtemp = gtemp.drop(gtemp.index[0])
    gtemp = gtemp.reset_index()

    for ind, line in enumerate(txt3):
        gtemp.ix[ind, 'stat params.'] = line

    gtemp = gtemp.convert_objects(convert_numeric=True)
#    gtemp.to_excel(writer,
#                   sheetname,
#                   columns=['group1',
#                            'group2',
#                            'pval',
#                            'pval_corr',
#                            'reject',
#                            'stat params.'])
    with open('output.txt','a') as f:
        f.write(u'\\multicolumn{6}{c}{'+sheetname+'}\n')
        f.write(res[0].as_latex_tabular().encode('utf8'))
#    print(gtemp.to_string(columns=['group1',
#                            'group2',
#                            'pval',
#                            'pval_corr',
#                            'reject',
#                            'stat params.'],
#                            index=False))
#    ypetrue = []
#    for i in res[2]:
#        if ((i[0] == 'YPE') | (i[1] == 'YPE')) & i[5] == True:
#            ypetrue.append(i)
    return (sheetname)


def convert_longform(datedgs, colname):
    """ for converting a dataframe of cell data whose values are lists of edges
    to long form with a media type column label
    """
    dftemp = pd.DataFrame({i: pd.Series(j) for i, j
                           in datedgs.ix[:, colname].iteritems()})
    dftemp = dftemp.stack().reset_index(1)
    dftemp.columns = ['cell', colname]
#    dftemp['cat'] = dftemp['cell'].apply(lambda x: x[:3])
    dfout = dftemp.reset_index(drop=True)
    return dfout


def pairplotter(dataset):
    """plot pairplot for single cell
    """
    for idx in dataset.index:
        dlist = {c: d for c, d in dataset.ix[idx].iteritems()}
        cell_conn = pd.DataFrame(dlist)
        cell_conn['cell'] = idx
        dfpair = cell_conn.ix[:, 2:8]
        gr1 = sns.PairGrid(dfpair)
        gr1 = gr1.map(sns.regplot)
        gr1.savefig('%s_pairplot.png' % idx)
        print '%s done' % idx
        plt.close()
    gr1 = gr1.map_lower(plt.scatter)
    gr1 = gr1.map_diag(sns.kdeplot, lw=3, legend=False)
    gr1 = gr1.map_upper(sns.regplot)
