# -*- coding: utf-8 -*-
"""
Created on Fri Jul 03 15:41:34 2015
Script to plot the autocorrelation coefficients of the various actual and
fitted DY distributions. Run lags\\MakeInputForLags.py in order to get the
fitted distributions pickle file ('*lagsunscaled)
@author: sweel
"""
import matplotlib.pyplot as plt
import os
import cPickle as pickle
import seaborn as sns
from tubuleHet.autoCor.AutoPopFunc import autocorout
from collections import defaultdict
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
sns.set_context("talk")
sns.set(style="whitegrid")
sns.set(rc={"legend.markerscale": 3})
# =============================================================================
#           Data initialization
# =============================================================================
plt.close('all')
# pylint: disable=C0103
dirlist = []
# pylint: enable=C0103
for root, dirs, files in os.walk(os.getcwd()):
    for f in dirs:
        if f.startswith('YP'):
            dirlist.append(
                os.path.join(root, f))

ACU = defaultdict(dict)  # uniform dist autocors
ACN = defaultdict(dict)  # Normal dist autocors
ACS = defaultdict(dict)  # Shuffled dist autocors
ACDY = defaultdict(dict)  # DY_scaled  autocors

for media in dirlist:
    labs = media[-3:]
    print'\nNow on %s' % labs+"\n"+"="*79
    with open('%s_lagsRFP.pkl' % labs, 'rb') as inpt:  # for RFP width
        (randNDY, randUDY, Norm, NormPermute, data) = pickle.load(inpt)

# =============================================================================
#           Main Function block
# =============================================================================
    for cell in data.keys():
        ACU[labs][cell] = []  # uniform dist autocors
        ACN[labs][cell] = []  # Normal dist autocors
        ACS[labs][cell] = []  # Shuffled dist autocors
        ACDY[labs][cell] = []  # DY_scaled  autocors

        # for autocor coeff
        ACDY[labs][cell].append(autocorout(cell, Norm))
        ACU[labs][cell].append(autocorout(cell, randUDY))
        ACN[labs][cell].append(autocorout(cell, randNDY))
        ACS[labs][cell].append(autocorout(cell, NormPermute))
# =============================================================================
# plot by type of distribution for YPE
# =============================================================================
f_lags = pd.DataFrame()
TYPE = {'actual YPE': ACDY['YPE'],
        'normal': ACN['YPE'],
        'shuffled': ACS['YPE'],
        'uniform': ACU['YPE']}
colors = ["medium green",
          "greyish blue",
          "yellowy brown",
          "reddish grey"]
# =============================================================================
# for autocor coeff of random vs real dist.
# =============================================================================
for mem in sorted(TYPE.keys()):
    autodata = TYPE[mem]
    temp = pd.concat([pd.Series(autodata[k][0], name=k) for k
                      in autodata.keys()], axis=1)
    autocorlags = []
    f_edges = pd.DataFrame()
    minl = 10
    maxl = 40
    for cll in temp.columns:
        temp2 = temp.ix[:, cll]
        temp2 = temp2.dropna()
        arrlen = temp2.apply(len)
        temp3 = pd.DataFrame({'auto_cor': temp2, 'len': arrlen})
        for thresh in np.linspace(minl, maxl, 4, endpoint=True):
            mask = (temp3['len'] >= thresh) & (temp3['len'] < thresh+10)
            temp3.loc[mask, ['len']] = thresh
        temp3.loc[temp3['len'] > maxl, ['len']] = maxl
        f_edges = f_edges.append(
            temp3.loc[temp3['len'] >= minl], ignore_index=True)

    for thresh in np.linspace(minl, maxl, 4, endpoint=True):
        bigf = f_edges.loc[f_edges.len == thresh]
        dftemp = pd.DataFrame({i: pd.Series(j) for i, j
                               in bigf.ix[:, 'auto_cor'].iteritems()})
        dftemp = dftemp.stack().reset_index(0)
        dftemp.columns = ['lag', 'auto_cor']
        dftemp['thresh'] = thresh
        dftemp['type'] = mem
        f_lags = f_lags.append(dftemp, ignore_index=True)

    f_lags.loc[(f_lags.lag > 11) & (f_lags.thresh == 10), ['thresh']] = None
    f_lags = f_lags[pd.notnull(f_lags['thresh'])]
    f_lags = f_lags.loc[f_lags.lag <= 15]

with sns.plotting_context('talk', font_scale=1.4):
    FIG_T = sns.factorplot(x='lag',
                           y='auto_cor',
                           col='thresh',
                           hue='type',
                           col_wrap=2,
                           palette=sns.xkcd_palette(colors),
                           scale=.5,
                           data=f_lags)
    plt.show()
    FIG_T.despine(left=True)
    FIG_T.set_ylabels('Autocorr. Coeff.')
    for subp in FIG_T.axes:
        subp.set_yticks(np.arange(-.25, 1.25, .25))
        subp.set_yticks(np.arange(-.25, 1.25, .25))
        subp.set_xticks(np.arange(0, 15, 2))
        subp.set_xticklabels(np.arange(0, 15, 2))

# =============================================================================
# for autocor corr for real dist by media and thresh length
# =============================================================================
f_lags = pd.DataFrame()
for mem in sorted(ACDY.keys()):
    autodata = ACDY[mem]
    temp = pd.concat([pd.Series(autodata[k][0], name=k) for k
                      in autodata.keys()], axis=1)
    autocorlags = []
    f_edges = pd.DataFrame()
    minl = 10
    maxl = 40
    for cll in temp.columns:
        temp2 = temp.ix[:, cll]
        temp2 = temp2.dropna()
        arrlen = temp2.apply(len)
        temp3 = pd.DataFrame({'auto_cor': temp2, 'len': arrlen})
        for thresh in np.linspace(minl, maxl, 4, endpoint=True):
            mask = (temp3['len'] >= thresh) & (temp3['len'] < thresh+10)
            temp3.loc[mask, ['len']] = thresh
        temp3.loc[temp3['len'] > maxl, ['len']] = maxl
        f_edges = f_edges.append(temp3.loc[temp3['len'] >= minl],
                                 ignore_index=True)

    for thresh in np.linspace(minl, maxl, 4, endpoint=True):
        bigf = f_edges.loc[f_edges.len == thresh]
        dftemp = pd.DataFrame({i: pd.Series(j) for i, j
                               in bigf.ix[:, 'auto_cor'].iteritems()})
        dftemp = dftemp.stack().reset_index(0)
        dftemp.columns = ['lag', 'auto_cor']
        dftemp['thresh'] = thresh
        dftemp['type'] = mem
        f_lags = f_lags.append(dftemp, ignore_index=True)

    f_lags.loc[(f_lags.lag > 11) & (f_lags.thresh == 10), ['thresh']] = None
    f_lags = f_lags[pd.notnull(f_lags['thresh'])]
    f_lags = f_lags.loc[f_lags.lag <= 15]

with sns.plotting_context('talk', font_scale=1.4):
    FIGM = sns.factorplot(x='lag',
                          y='auto_cor',
                          col='thresh',
                          hue='type',
                          col_wrap=2,
                          scale=.5,
                          data=f_lags)
    plt.show()
    FIGM.despine(left=True)
    FIGM.set_ylabels('Autocorr. Coeff.')
    for subp in FIGM.axes:
        subp.set_xticks(np.arange(0, 15, 2))
        subp.set_xticklabels(np.arange(0, 15, 2))

with sns.plotting_context('talk', font_scale=1.4):
    _, axes0 = plt.subplots(1, 1)
    sns.pointplot(x='lag',
                  y='auto_cor',
                  hue='type',
                  scale=.75,
                  data=f_lags.loc[f_lags.thresh == 40],
                  ax=axes0)
    plt.show()
    axes0.set_ylabel('Autocorr. Coeff.')
    axes0.set_xticks(np.arange(0, 15, 2))
    axes0.set_xticklabels(np.arange(0, 15, 2))

# =============================================================================
# curve fitting exponential
# =============================================================================


def func(x, b):
    return np.exp(-b * x)

for mem in sorted(ACDY.keys()):
    data = f_lags.loc[(f_lags.type == mem) & (f_lags.thresh == 40)]
    xdata = data.lag
    ydata = data.auto_cor
    popt, pcov = curve_fit(func, xdata, ydata)
    perr = np.sqrt(np.diag(pcov))

    print 'alpha for %s : %6.4f with std=%6.4f' % (mem, popt, perr)

popt, pcov = curve_fit(func, xdata, ydata)
with open('autocorRFP.pkl', 'wb') as output:
    pickle.dump(f_lags, output)

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
