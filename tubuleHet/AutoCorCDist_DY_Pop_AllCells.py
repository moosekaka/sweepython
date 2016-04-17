# -*- coding: utf-8 -*-
"""
Created on Fri Jul 03 15:41:34 2015
Script to plot the autocorrelation coefficients of the various actual and
fitted DY distributions. Run lags\\MakeInputForLags.py in order to get the
fitted distributions pickle file ('*lagsunscaled)
@author: sweel
"""
import os
import os.path as op
import cPickle as pickle
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import wrappers as wr
from tubuleHet.autoCor.AutoPopFunc import autocorout
sns.set_context("talk")
sns.set(style="darkgrid")
sns.set(rc={"legend.markerscale": 3})


def binedges(coldata, interval):
    """
    Bins the edge lists according to length of edge array
    """
    coldata = coldata.dropna()
    edgelen = coldata.apply(len)
    edgs = pd.DataFrame({'auto_cor': temp, 'len': edgelen})
    edgs['cut'] = pd.cut(edgs.len,
                         interval,
                         labels=interval[:-1],
                         include_lowest=True,
                         right=False)
    edgs = edgs.dropna()  # drops edges < min len
    return edgs


def lagdist(edgebinned, thresh, label):
    """
    returns an array of lag distances (k) with the autocor coef at that k
    """
    binned_edges = edgebinned.loc[edgebinned.cut == thresh]
    dftemp = pd.DataFrame({i: pd.Series(j) for i, j
                           in binned_edges.ix[:, 'auto_cor'].iteritems()})

    dftemp = dftemp.stack().reset_index(0)
    dftemp.columns = ['lag', 'auto_cor']
    dftemp['thresh'] = thresh
    dftemp['type'] = label
    return dftemp


# =============================================================================
#           Data initialization
# =============================================================================
plt.close('all')
datadir = op.join(os.getcwd(), 'data')
rawdir = op.join(os.getcwd(), 'output')
vtkF = wr.ddwalk(op.join(rawdir, 'normalizedVTK'),
                 '*skeleton.vtk', start=5, stop=-13)

ACU = defaultdict(dict)  # uniform dist autocors
ACN = defaultdict(dict)  # Normal dist autocors
ACS = defaultdict(dict)  # Shuffled dist autocors
ACDY = defaultdict(dict)  # DY_scaled  autocors


# =============================================================================
#           Main Function block
# =============================================================================
for mtype in sorted(vtkF.keys())[:]:
    for cell in vtkF[mtype].keys():
        with open(op.join(rawdir,
                          'fitted_data',
                          '%s.pkl' % cell), 'rb') as inpt:
            (Norm, NormPermute, randNDY, randUDY, lineId) = pickle.load(inpt)
        ACU[mtype][cell] = []
        ACN[mtype][cell] = []
        ACS[mtype][cell] = []
        ACDY[mtype][cell] = []

        ACDY[mtype][cell].append(autocorout(Norm[0]))
        ACU[mtype][cell].append(autocorout(randUDY[0]))
        ACN[mtype][cell].append(autocorout(randNDY[0]))
        ACS[mtype][cell].append(autocorout(NormPermute[0]))
        print "done calculating autocor for %s" % cell

# =============================================================================
# Calculate population autocorr coef of YPE for real vs random
# =============================================================================
autocor_type = {'actual YPE': ACDY['YPE'],
                'normal': ACN['YPE'],
                'shuffled': ACS['YPE'],
                'uniform': ACU['YPE']}
colors = ["medium green",
          "greyish blue",
          "yellowy brown",
          "reddish grey"]
real_rand_lags = pd.DataFrame()


for dist_type in sorted(autocor_type.keys()):
    autodata_cell = autocor_type[dist_type]
    autodata_pop = pd.concat([pd.Series(np.squeeze(autodata_cell[k]), name=k)
                              for k in autodata_cell.keys()], axis=1)
    real_rand_edges = pd.DataFrame()  # holds the edges binned by length
    minl = 10
    maxl = 40
    # set the open right bound to be longer then the longest edge
    bins = np.r_[np.linspace(minl, maxl, num=4, endpoint=True), 200]

    # bin the edges according to edgelen
    for col in autodata_pop.columns:
        temp = autodata_pop.ix[:, col]
        real_rand_edges = real_rand_edges.append(binedges(temp, bins))
        real_rand_edges = real_rand_edges.reset_index(drop=True)

    # group by lag distances
    for threshold_len in bins[:-1]:
        lags_k = lagdist(real_rand_edges,
                         threshold_len,
                         dist_type)
        real_rand_lags = real_rand_lags.append(lags_k,
                                               ignore_index=True)

# =============================================================================
# Calculate population autocorr coef by carbon type
# =============================================================================
ferm_resp_lags = pd.DataFrame()

for carbon in sorted(ACDY.keys()):
    autodata_cell = ACDY[carbon]
    autodata_pop = pd.concat([pd.Series(np.squeeze(autodata_cell[k]), name=k)
                              for k in autodata_cell.keys()], axis=1)
    ferm_resp_edges = pd.DataFrame()  # holds the edges binned by length
    minl = 10
    maxl = 40
    # set the open right bound to be longer then the longest edge
    bins = np.r_[np.linspace(minl, maxl, num=4, endpoint=True), 200]

    # bin the edges according to edgelen
    for col in autodata_pop.columns:
        temp = autodata_pop.ix[:, col]
        ferm_resp_edges = ferm_resp_edges.append(binedges(temp, bins))
        ferm_resp_edges = ferm_resp_edges.reset_index(drop=True)

    # group by lag distances
    for threshold_len in bins[:-1]:
        lags_k = lagdist(ferm_resp_edges,
                         threshold_len,
                         carbon)
        ferm_resp_lags = ferm_resp_lags.append(lags_k,
                                               ignore_index=True)

#    ================================================================
    # Plots
#    ================================================================
# real vs random
with sns.plotting_context('talk', font_scale=1.4):
    real_rand_lags = real_rand_lags.loc[real_rand_lags.lag <= 15]
    FIG_T = sns.factorplot(x='lag',
                           y='auto_cor',
                           col='thresh',
                           hue='type',
                           col_wrap=2,
                           palette=sns.xkcd_palette(colors),
                           scale=.5,
                           data=real_rand_lags)
    plt.show()
    FIG_T.despine(left=True)
    FIG_T.set_ylabels('Autocorr. Coeff.')
    for subp in FIG_T.axes:
        subp.set_yticks(np.arange(-.25, 1.25, .25))
        subp.set_yticks(np.arange(-.25, 1.25, .25))
        subp.set_xticks(np.arange(0, 15, 2))
        subp.set_xticklabels(np.arange(0, 15, 2))

# carbon type
with sns.plotting_context('talk', font_scale=1.4):
    ferm_resp_lags = ferm_resp_lags.loc[ferm_resp_lags.lag <= 15]
    FIGM = sns.factorplot(x='lag',
                          y='auto_cor',
                          col='thresh',
                          hue='type',
                          col_wrap=2,
                          scale=.5,
                          data=ferm_resp_lags)
    plt.show()
    FIGM.despine(left=True)
    FIGM.set_ylabels('Autocorr. Coeff.')
    for subp in FIGM.axes:
        subp.set_xticks(np.arange(0, 15, 2))
        subp.set_xticklabels(np.arange(0, 15, 2))

# =============================================================================
# curve fitting exponential
# =============================================================================

#def func(x, b):
#    ''' fit an exponential function to vect x
#    '''
#    return np.exp(-b * x)
#
#for mem in sorted(ACDY.keys()):
#    data = ferm_resp_edges.loc[(ferm_resp_edges.type == mem) & (ferm_resp_edges.thresh == 40)]
#    xdata = data.lag
#    ydata = data.auto_cor
#    popt, pcov = curve_fit(func, xdata, ydata)
#    perr = np.sqrt(np.diag(pcov))
#
#    print 'alpha for %s : %6.4f with std=%6.4f' % (mem, popt, perr)
#
#popt, pcov = curve_fit(func, xdata, ydata)

#with open('autocorNorm.pkl', 'wb') as output:
#    pickle.dump(ferm_resp_edges, output)
