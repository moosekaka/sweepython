# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 15:55:26 2015
CAlculates and plots POWER SPEC DENSTIY
@author: sweel
"""
import matplotlib.pyplot as plt
import os
import cPickle as pickle
import seaborn as sns
from tubuleHet.autoCor.AutoPopFunc import psd
from collections import defaultdict
import pandas as pd
import numpy as np
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

PSDY = defaultdict(dict)
PSDP = defaultdict(dict)
PSDN = defaultdict(dict)
PSDU = defaultdict(dict)

for media in dirlist:
    labs = media[-3:]
    print'\nNow on %s' % labs + "\n" + "=" * 79
    with open('%s_lagsUnscaled.pkl' % labs, 'rb') as inpt:
        (randNDY, randUDY, Norm, NormPermute, data) = pickle.load(inpt)[:-1]

# =============================================================================
#           Main Function block
# =============================================================================
    for cell in data.keys():
        PSDY[labs][cell] = []  # DY
        PSDP[labs][cell] = []  # Shuffled
        PSDN[labs][cell] = []  # Normal
        PSDU[labs][cell] = []  # uniform

        # for psd values
        PSDY[labs][cell].append(psd(cell, Norm, 40))
        PSDP[labs][cell].append(psd(cell, NormPermute, 40))
        PSDN[labs][cell].append(psd(cell, randNDY, 40))
        PSDU[labs][cell].append(psd(cell, randUDY, 40))
# =============================================================================
# Power spectrum density actual (YPE) vs random
# =============================================================================
colors = ["medium green",
          "greyish blue",
          "yellowy brown",
          "reddish grey"]
TYPE2 = {'actual YPE': PSDY['YPE'],
         'normal': PSDN['YPE'],
         'shuffled': PSDP['YPE'],
         'uniform': PSDU['YPE']}

PSD = pd.DataFrame()
for mem in sorted(TYPE2.keys()):

    autodata = TYPE2[mem]
    x1 = pd.concat([pd.Series(autodata[k][0][0], name=k) for k
                    in autodata.keys()], axis=1)
    y1 = pd.concat([pd.Series(autodata[k][0][1], name=k) for k
                    in autodata.keys()], axis=1)
    p_edges = pd.DataFrame()

    for cll in x1.columns:
        ptemp = y1.ix[:, cll]
        xtemp = x1.ix[:, cll]
        ptemp = ptemp.dropna()
        xtemp = xtemp.dropna()
        psdtemp = pd.DataFrame({'psd': ptemp, 'u': xtemp})
        p_edges = p_edges.append(psdtemp, ignore_index=True)
        psdx = pd.DataFrame({i: pd.Series(j) for i, j
                             in p_edges.ix[:, 'u'].iteritems()})
        psdy = pd.DataFrame({i: pd.Series(j) for i, j
                             in p_edges.ix[:, 'psd'].iteritems()})

    G = []
    bins = np.linspace(0, .5, 22)
    for edge in psdx.columns:
        A = psdy.ix[:, edge].groupby(pd.cut(psdx.ix[:, edge], bins)).mean()
        G.append(A.values)
    psd = pd.concat([pd.Series(i) for i in G], axis=1)
    psd['u'] = np.round(bins[1:], 2)
    psd2 = pd.melt(psd,
                   id_vars='u',
                   var_name='lineid',
                   value_name='psd')
    psd2['type'] = mem
    PSD = PSD.append(psd2, ignore_index=True)


with sns.plotting_context('talk', font_scale=1.4):
    _, axes5 = plt.subplots(1, 1)
    sns.pointplot(x='u',
                  y='psd',
                  hue='type',
                  palette=sns.xkcd_palette(colors),
                  scale=.95,
                  data=PSD,
                  ax=axes5)
    axes5.get_legend().set_visible(False)

# =============================================================================
# Power spectrum density actual by media
# =============================================================================
PSD2 = pd.DataFrame()
for mem in sorted(PSDY.keys()):

    autodata = PSDY[mem]
    x1 = pd.concat([pd.Series(autodata[k][0][0], name=k) for k
                    in autodata.keys()], axis=1)
    y1 = pd.concat([pd.Series(autodata[k][0][1], name=k) for k
                    in autodata.keys()], axis=1)
    p_edges = pd.DataFrame()

    for cll in x1.columns:
        ptemp = y1.ix[:, cll]
        xtemp = x1.ix[:, cll]
        ptemp = ptemp.dropna()
        xtemp = xtemp.dropna()
        psdtemp = pd.DataFrame({'psd': ptemp, 'u': xtemp})
        p_edges = p_edges.append(psdtemp, ignore_index=True)
        psdx = pd.DataFrame({i: pd.Series(j) for i, j
                             in p_edges.ix[:, 'u'].iteritems()})
        psdy = pd.DataFrame({i: pd.Series(j) for i, j
                             in p_edges.ix[:, 'psd'].iteritems()})

    G = []
    bins = np.linspace(0, .5, 22)
    for edge in psdx.columns:
        A = psdy.ix[:, edge].groupby(pd.cut(psdx.ix[:, edge], bins)).mean()
        G.append(A.values)
    psd = pd.concat([pd.Series(i) for i in G], axis=1)
    psd['u'] = np.round(bins[1:], 2)
    psd2 = pd.melt(psd,
                   id_vars='u',
                   var_name='lineid',
                   value_name='psd')
    psd2['type'] = mem
    PSD2 = PSD2.append(psd2, ignore_index=True)
sns.set(style="whitegrid")
with sns.plotting_context('talk', font_scale=1.4):
    _, axes6 = plt.subplots(1, 1)
    sns.pointplot(x='u',
                  y='psd',
                  hue='type',
                  scale=.95,
                  data=PSD2,
                  ax=axes6)
    axes6.get_legend().set_visible(False)
