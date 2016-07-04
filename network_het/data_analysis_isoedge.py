# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 17:40:33 2015

@author: sweel
"""
import pandas as pd
import cPickle as pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as sp
sns.plotting_context('talk', font_scale=1.4)
# pylint: disable=C0103
plt.close('all')


with open('munged_dataframe.pkl', 'rb') as INPUT:
    df = pickle.load(INPUT)
with open('lagedges.pkl', 'rb') as INPT:
    dflags = pickle.load(INPT)
df['lags_1'] = dflags
with open('cellVolume.pkl', 'rb') as INPT:
    dfsize = pickle.load(INPT)
    dfsize.index = dfsize.cell
with open('dic_names.pkl', 'rb') as inpt:
    dic = pickle.load(inpt)

dfscals = df.loc[:, ['mito_avgdeg',
                     'mito_cell_avedy',
                     'mito_cell_avedyr',
                     'mito_cell_stddy',
                     'mito_cell_stddyr',
                     'mito_charpl_uw',
                     'mito_charpl_w',
                     'mito_edgenum',
                     'mito_beta_geo',
                     'mito_beta_top',
                     'mito_phi',
                     'mito_pk3',
                     'mito_totlen',
                     'charpl_norm_len',
                     'charpl_norm_numedge',
                     'cell_coefvar',
                     'cell_coefvar_r',
                     'mito_cell_ave_gfp',
                     'mito_iso_dyr']]

dfvol = pd.DataFrame({'Vol': dfsize.loc[:, 'Vol'],
                      'mitolen': df.loc[:, 'mito_totlen'],
                      'media': df.loc[:, 'media']})
dfvol['mitovol'] = np.pi * (.15)**2 * dfvol.mitolen
dfvol['Vol Ratio'] = dfvol.mitovol / dfvol.Vol


dfscals = pd.concat([dfscals,
                     dfvol.loc[:, 'Vol Ratio'],
                     df.loc[:, 'media']], axis=1)
temp3 = []
for cell in dfscals.index:
    if pd.notnull(dfscals.mito_iso_dyr.loc[cell]):
        temp = dfscals.mito_iso_dyr.loc[cell]
        temp2 = pd.DataFrame({k: pd.Series(val) for k, val
                              in temp.items()})
        temp2.columns = [cell for i in temp2.columns]
        temp3.append(temp2)
temp4 = pd.concat(temp3, axis=1)

# =============================================================================
# isolated chunks analysis
# =============================================================================
dfchunks2 = pd.concat([temp4.mean(), temp4.count()], axis=1)
dfchunks2.columns = ['mean_isody', 'len_isody']
dfchunks2 = pd.concat([dfchunks2, dfscals.mito_cell_avedyr],
                      join='inner',
                      axis=1)
dfchunks2.reset_index(inplace=True)
dfchunks2['media'] = dfchunks2['index'].apply(lambda x: x[:3])
dfchunks2['thresh'] = None
dfchunks2.loc[dfchunks2.len_isody >= 9.5, ['thresh']] = 'long'
dfchunks2.loc[dfchunks2.len_isody < 9.5, ['thresh']] = 'short'
dfm2 = pd.melt(dfchunks2,
               value_vars=['mean_isody',
                           'len_isody',
                           'mito_cell_avedyr'],
               id_vars=['media', 'thresh'],
               value_name='value')

with sns.plotting_context('talk', font_scale=1.25):
    fig1, ax1 = plt.subplots(1, 1)
    sns.violinplot(x='media',
                   y='value',
                   hue='thresh',
                   data=dfm2[dfm2.variable == 'mean_isody'],
                   ax=ax1)
    sns.stripplot(x='media',
                  y='value',
                  hue='thresh',
                  data=dfm2[dfm2.variable == 'mean_isody'],
                  ax=ax1,
                  jitter=.2,
                  size=3,
                  alpha=.5)
    ax1.set_ylim(0, 2000)
    ax1.get_legend().set_visible(False)
    h, l = ax1.get_legend_handles_labels()
    h = h[:2]
    l = l[:2]
    ax1.legend(h, l, loc=2)
    fig1.suptitle(r'Long vs short isolated edges $\Delta \Psi$ distribution',
                  fontsize=22)

dfm3 = dfm2.loc[(dfm2.variable == 'mean_isody') |
                (dfm2.variable == 'mito_cell_avedyr')]
with sns.plotting_context('talk', font_scale=1.25):
    fig2, ax2 = plt.subplots(1, 1)
    sns.violinplot(x='media',
                   y='value',
                   hue='variable',
                   data=dfm3,
                   ax=ax2)
    ax2.set_ylim(0, 2000)
    ax2.get_legend().set_visible(False)
    h, l = ax2.get_legend_handles_labels()
    h = h[:2]
    l = l[:2]
    ax2.legend(h, l, loc=2)
    fig2.suptitle(r'Cell vs isolated edges $\Delta \Psi$ distribution',
                  fontsize=22)


dfchunks2['isol.egde len ($\mu m$)'] = dfchunks2.len_isody*.055
with sns.plotting_context('talk', font_scale=1.25):
    g = sns.lmplot(x='isol.egde len ($\mu m$)',
                   y='mean_isody',
                   col="media",
                   hue='media',
                   data=dfchunks2,
                   col_wrap=2,
                   size=5,
                   scatter_kws={"s": 25, "alpha": 1})
    g.set(xlim=(0.))
    g.set(ylim=(0, 2000))
    for j in g.axes:
        j.get_lines()[0].set_visible(False)


for mem in sorted(pd.unique(dfchunks2.media)):
    mask = dfchunks2.media == mem
    dftest = dfchunks2.loc[mask]
    for mea in ['mean_isody', 'mito_cell_avedyr']:
        res = sp.pearsonr(dftest.loc[mask, 'isol.egde len ($\mu m$)'],
                          dftest.loc[mask, mea])
        print '%-4s:%-25s:%8.3f %8.3f' % (mem, mea, res[0], res[1])

for mem in sorted(pd.unique(dfchunks2.media)):
    mask = dfchunks2.media == mem
    dftest = dfchunks2.loc[mask]
    res = sp.pearsonr(dftest.loc[mask, 'isol.egde len ($\mu m$)'],
                    dftest.loc[mask, 'mean_isody'])
    print '%-4s:%8.3f %8.3f' % (mem, res[0], res[1])