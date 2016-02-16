# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 12:08:52 2015
O2 data graphs
@author: sweel
"""

import pandas as pd
import cPickle as pickle
import scipy.stats as sp
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mungedata import MungeDataFuncs as md
sns.plotting_context('talk', font_scale=1.4)
# pylint: disable=C0103
plt.close('all')

with open('o2data.pkl', 'rb') as INPUT:
    dfo2 = pickle.load(INPUT)
with open('munged_dataframe.pkl', 'rb') as INPUT:
    df = pickle.load(INPUT)
with open('lagedges.pkl', 'rb') as INPT:
    dflags = pickle.load(INPT)
df['lags_1'] = dflags
with open('cellVolume.pkl', 'rb') as INPT:
    dfsize = pickle.load(INPT)
    dfsize.index = dfsize.cell

dfvol = pd.DataFrame({'Vol': dfsize.loc[:, 'Vol'],
                      'mitolen': df.loc[:, 'mito_totlen'],
                      'media': df.loc[:, 'media']})

dfvol['mitovol'] = np.pi * (.15)**2 * dfvol.mitolen
dfvol['Vol Ratio'] = dfvol.mitovol / dfvol.Vol

#   Topology  (by cell)
dfcellcon = df.ix[:, ['mito_cell_avedyr',
                      'cell_coefvar_r',
                      'mito_beta_geo',
                      'mito_beta_top',
                      'mito_phi',
                      'mito_pk3',
                      'mito_avgdeg',
                      'mito_edgenum']]
dfcellcon = pd.concat([dfcellcon,
                       dfvol.loc[:, 'Vol Ratio'],
                       df.loc[:, 'media']], axis=1)

dfcellcon = dfcellcon[dfcellcon.mito_cell_avedyr <= 2000]  # exclude hi YPE's
dfcellcon[r'$\Delta \Psi$ Unscaled'] = dfcellcon["mito_cell_avedyr"]

# =============================================================================
#   O2 data here
# =============================================================================
dfcellcon['Number of Edges'] = dfcellcon.mito_edgenum
dfcellcon['Average Degree'] = dfcellcon.mito_avgdeg
dfcellcon['O2 per mito vol'] = ''
dfcellcon['OCR per cell'] = ''

#with sns.plotting_context('talk', font_scale=1.5):
#    g = sns.lmplot(x="Vol Ratio",
#                   y='Average Degree',
#                   col="media",
#                   hue='media',
#                   data=dfcellcon,
#                   col_wrap=2,
#                   size=5,
#                   scatter_kws={"s": 25, "alpha": 1})
#    g.set(ylim=(1, 3.5))


with sns.plotting_context('talk', font_scale=1.5):
    f, (ax1, ax0) = plt.subplots(2, 1,
                                 figsize=(11, 8.5),
                                 sharex=True)
    sns.violinplot('media',
                   'Vol Ratio',
                   data=dfcellcon,
                   hue='media',
                   ax=ax0)
    ax0.set_xlabel('')
    ax0.set_ylabel('Volume Ratio')
    ax0.get_legend().set_visible(False)

    A=dfo2.groupby('type').quantile(0.025).values.flatten()
    B=dfo2.groupby('type').quantile(0.975).values.flatten()
    C=dfo2.groupby('type').quantile(0.5).values.flatten()
    g = sns.barplot(x='type',
                    y='OCRmito',
                    estimator=np.median,
                    ci=None,
                    ecolor=[.25, .25, .25],
                    data=dfo2,
                    yerr=[C-A,B-C],
                    ax=ax1)
    ax1.set_ylabel(r'OCR per mito vol /$\mu m^{3}$')
    ax1.set_xlabel('')
#    ax1.get_legend().set_visible(False)


with sns.plotting_context('talk', font_scale=1.5):
    f, (ax1, ax3, ax2) = plt.subplots(3, 1,
                                      figsize=(11, 8.5),
                                      sharex=True)
    g = sns.barplot(x='type',
                    y='OCRmito',
                    estimator=np.median,
                    ci=None,
                    ecolor=[.25, .25, .25],
                    data=dfo2,
                    yerr=[C-A,B-C],
                    ax=ax1)

    ax1.set_ylabel(r'OCR per mito vol /$\mu m^{3}$')
    ax1.set_xlabel('')
    ax1.set_yticks(np.arange(0,.015,.005))
#    ax1.get_legend().set_visible(False)

    sns.violinplot('media',
                   'mito_cell_avedyr',
                   data=dfcellcon,
                   hue='media',
                   ax=ax2)
    ax2.set_ylabel(r'$\Delta \Psi$ ')
    ax2.set_ylim(0, 2000)
    ax2.set_xlabel('')
    ax2.get_legend().set_visible(False)

    sns.violinplot('media',
                   'Vol Ratio',
                   data=dfcellcon,
                   hue='media',
                   estimator=np.median,
                   ax=ax3)
    ax3.set_ylabel('Vol Ratio')
    ax3.set_ylim(0)
    ax3.set_xlabel('')
    ax3.get_legend().set_visible(False)
    sns.despine(bottom=True)
    plt.tight_layout(h_pad=3)
