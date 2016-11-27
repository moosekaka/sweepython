# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 10:55:07 2015

@author: sweel
"""
import os.path as op
import pandas as pd
import cPickle as pickle
import scipy.stats as sp
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from network_het.mungedata import MungeDataFuncs as md
import statsmodels.formula.api as sm

sns.plotting_context('talk', font_scale=1.4)
sns.set(style="whitegrid")
# pylint: disable=C0103b
plt.close('all')
olddir = op.join('deprecated', 'OLD_NormFiles_021416')

with open('munged_dataframe_2016.pkl', 'rb') as INPUT:
    df = pickle.load(INPUT)
with open(op.join('cellVolume.pkl'), 'rb') as INPT:
    dfsize = pickle.load(INPT)
    dfsize.index = dfsize.cell
with open(op.join(olddir, 'dic_names.pkl'), 'rb') as inpt:
    dic = pickle.load(inpt)
with open('list_wtmb.pkl', 'rb') as inpt:
    wt = pickle.load(inpt)

wtall = df[df.media == 'WT'].index
wtfilt = df.index.intersection(wt)
exclude = wtall.difference(wtfilt)
df = df.loc[df.index.difference(exclude)]
dfsize = dfsize.loc[df.index.difference(exclude)]
newlabels = {'YPR': u'WT_YPR', 'YPL': u'WT_YPL', 'YPE': u'WT_YPE',
             'WT': u'WT_YPE', 'NUM1': u'ΔNUM1', 'MFB1': u'ΔMFB1',
             'YPT11': u'ΔYPT11', 'YPD': 'WT_YPD'}
wtypes = ['WT_YPD', 'WT_YPE', 'WT_YPL', 'WT_YPR']

df['media_new'] = df.media.map(newlabels)

newcollabels = {'mito_avgdeg': 'Average Degree',
                'mito_edgenum': 'Number of Edges',
                'mito_beta_geo': u'β geo.',
                'mito_beta_top': u'β top.',
                'mito_pk3': u'pk₃',
                'mito_phi': u'ϕ',
                'mito_bpts_dyraw': u'branchpoints ΔΨ',
                'mito_knn_uw': 'Nearest Neighbor Degree',
                'mito_cell_avedyr': u'ΔΨ Unscaled',
                'mito_btwcntr_uw': 'Betweeness Centr.',
                'mito_clscntr_uw': 'Closeness Centr.',
                'mito_clstcf_uw': 'Clustering Coeff.'}
df.rename(columns=newcollabels, inplace=True)

# =============================================================================
# setup subset of dataframes for stat tests
# =============================================================================
# label/filter lists for global and local connectivity
globcon = [u'ΔΨ Unscaled',
           'Amount density',
           u'β geo.',
           u'β top.',
           u'ϕ',
           u'pk₃',
           'Average Degree',
           'Number of Edges',
           'mito_tubew']

localcon = ['Amount density',
            u'branchpoints ΔΨ',
            'Betweeness Centr.',
            'Nearest Neighbor Degree',
            'Closeness Centr.',
            'Clustering Coeff.']

dfvol = pd.DataFrame({'Vol': dfsize.loc[:, 'Vol'],
                      'Surf': dfsize.loc[:, 'Surf'],
                      'mitolen': df.loc[:, 'mito_totlen'],
                      'media_new': df.loc[:, 'media_new']})

dfvol['mitovol'] = np.pi * (.15)**2 * dfvol.mitolen
dfvol['Vol Ratio'] = dfvol.mitovol / dfvol.Vol
dfvol['Amount density'] = dfvol.mitolen / dfvol.Surf

with open('cellVolume_complete.pkl', 'wb') as out:
    pickle.dump(dfvol, out)

# Network Connectivity Dataframe (by branchpoints)
dfconn = df.ix[:, localcon[1:]]
frames = []
for i in dfconn.columns:
    temp = md.convert_longform(dfconn, i)
    temp.index = temp.cell
    temp = temp.drop('cell', 1)
    frames.append(temp)
BIG = pd.concat(frames, axis=1)
BIG = BIG.groupby(BIG.index).mean()

dfconn = pd.concat([BIG,
                    dfvol.loc[:, ['Vol Ratio', 'Amount density']],
                    df.loc[:, 'media_new']], axis=1)

# Topology  (by cell)
dfcellcon = df.ix[:, [i for i in globcon if i != 'Amount density']]
dfcellcon = pd.concat([dfcellcon,
                       dfvol.loc[:, ['Vol Ratio', 'Amount density']],
                       df.loc[:, 'media_new']], axis=1)

dfcellcon = dfcellcon[dfcellcon[u'ΔΨ Unscaled'] <= 2000]  # exclude hi YPE's

# subset by mito amount density
subset = ((dfcellcon['Amount density'] > .5) &
          (dfcellcon['Amount density'] < .9))
dfsubset = dfcellcon[subset]
a = dfsubset.groupby('media_new').count()
b = dfcellcon.groupby('media_new').count()
c = 1. * a['Amount density'] / b['Amount density']

# =============================================================================
# table 5.3
# =============================================================================
print ("5.3 Test for local connectivity vs "
       "mito. amount density\n{}".format('*' * 79))
data = []
for mem in wtypes:
    mask = dfconn.media_new == mem
    dftest = dfconn.loc[mask]
    for mea in localcon:
        res = sp.pearsonr(dftest.loc[mask, 'Amount density'],
                          dftest.loc[mask, mea])
        data.append((mem, mea, res[0], res[1]))
dfl = pd.DataFrame(data)
dfl.columns = ['media_new', 'type', 'r', 'p']
dfl = dfl.pivot(index='type', columns='media_new', values='r')
dfl.to_clipboard(float_format='%6.4f')

# =============================================================================
# table 5.5
# =============================================================================
data = []
print u"5.5 Test for local conn. vs ΔΨ branchpoints\n{}".format('*' * 79)
for mem in wtypes:
    mask = dfconn.media_new == mem
    dftest = dfconn.loc[mask]
    for mea in localcon:
        res = sp.pearsonr(dftest.loc[mask, u'branchpoints ΔΨ'],
                          dftest.loc[mask, mea])
        data.append((mem, mea, res[0], res[1]))
dfl2 = pd.DataFrame(data)
dfl2.columns = ['media_new', 'type', 'r', 'p']
dfl2 = dfl2.pivot(index='type', columns='media_new', values='r')
dfl2.to_clipboard(float_format='%6.4f')

# =============================================================================
# table 5.2
# =============================================================================
print ("5.2 Test for global connectivity vs "
       "mito. amount density\n{}".format('*' * 79))
data1 = []
for mem in wtypes:
    mask = dfcellcon.media_new == mem
    dftest = dfcellcon.loc[mask]
    for mea in globcon:
        res = sp.pearsonr(dftest.loc[mask, 'Amount density'],
                          dftest.loc[mask, mea])
        data1.append((mem, mea, res[0], res[1]))

dfg = pd.DataFrame(data1)
dfg.columns = ['media_new', 'type', 'r', 'p']
dfg = dfg.pivot(index='type', columns='media_new', values='r')
dfg.to_clipboard(float_format='%6.4f')

# =============================================================================
# Table 5.4
# =============================================================================
print u"5.4 Test for global connectivity vs ΔΨ Unscaled\n{}".format('*' * 79)
data1 = []
for mem in wtypes:
    mask = dfcellcon.media_new == mem
    dftest = dfcellcon.loc[mask]
    for mea in globcon:
        res = sp.pearsonr(dftest.loc[mask, u'ΔΨ Unscaled'],
                          dftest.loc[mask, mea])
        data1.append((mem, mea, res[0], res[1]))

dfg2 = pd.DataFrame(data1)
dfg2.columns = ['media_new', 'type', 'r', 'p']
dfg2 = dfg2.pivot(index='type', columns='media_new', values='r')
dfg2.to_clipboard(float_format='%6.4f')

# print (u"\npercentage of population in subset [{:4.2f}-{:4.2f}]:"
#       .format(dfsubset.Amount density.min(), dfsubset.Amount density.max()))
# print "{}".format(c)
# =============================================================================
# Table 5.6
# =============================================================================
print ("5.6 Test for subset [{:4.2f}-{:4.2f}] of "
       "glob conn. vs Amount density \n{}"
       .format(dfsubset['Amount density'].min(),
               dfsubset['Amount density'].max(), '*' * 79))
data2 = []
for mea in globcon:
    for mem in wtypes[1:]:
        mask = dfcellcon.media_new == mem
        res = sp.pearsonr(dfsubset.loc[mask, 'Amount density'],
                          dfsubset.loc[mask, mea])
        data2.append((mem, mea, res[0], res[1]))
dfs = pd.DataFrame(data2)
dfs.columns = ['media_new', 'type', 'r', 'p']
dfs = dfs.pivot_table(index='type', columns='media_new', values='r')
dfs.to_clipboard(float_format='%6.4f')

# =============================================================================
# Table 5.7
# =============================================================================
print ("5.7 Test for subset [{:4.2f}-{:4.2f}] of "
       "local conn. vs Amount density \n{}"
       .format(dfsubset['Amount density'].min(),
               dfsubset['Amount density'].max(), '*' * 79))
subset2 = (dfconn['Amount density'] > .5) & (dfconn['Amount density'] < .9)
dfsubset2 = dfconn[subset2]
data3 = []
for mea in localcon:
    for mem in wtypes[1:]:
        mask = dfconn.media_new == mem
        res = sp.pearsonr(dfsubset2.loc[mask, 'Amount density'],
                          dfsubset2.loc[mask, mea])
        data3.append((mem, mea, res[0], res[1]))
dfs2 = pd.DataFrame(data3)
dfs2.columns = ['media_new', 'type', 'r', 'p']
dfs2 = dfs2.pivot_table(index='type', columns='media_new', values='r')
dfs2.to_clipboard(float_format='%6.4f')

# ============================================================================
# # PLOTS
# =============================================================================
g = md.boxviol(dfvol, 'Amount density', 'media_new')
g.set_title('Quasi Density', fontsize=24)
g.set_ylabel('Quasi Density')
plt.ylim([0, dfvol['Amount density'].quantile(0.995)])
h = md.boxviol(dfvol, 'Surf', 'media_new')
h.set_title(u'Cell Surface Area /µm²', fontsize=24)
h.set_ylabel('')
plt.ylim([0, dfvol.Surf.quantile(0.9995)])

# correlations for Amount density vs connectivity
with sns.plotting_context('talk', font_scale=1.25):
    g = sns.lmplot(x="Amount density",
                   y='Average Degree',
                   hue='media_new',
                   data=dfcellcon,
                   size=5,
                   scatter_kws={"s": 25, "alpha": 1})
    g.set(ylim=(0, 3.5))

with sns.plotting_context('talk', font_scale=1.25):
    g = sns.lmplot(x="Amount density",
                   y=u'ΔΨ Unscaled',
                   col='media_new',
                   hue='media_new',
                   data=dfcellcon,
                   col_wrap=2,
                   size=5,
                   scatter_kws={"s": 25, "alpha": 1})
    g.set(ylim=(0, 2000), xlim=(0))
    g.fig.suptitle(u'ΔΨ vs Mito. amount density', fontsize=22)
    g.fig.subplots_adjust(top=.9)
    for i in g.axes.flat:
        i.get_lines()[0].set_visible(False)

with sns.plotting_context('talk', font_scale=1.25):
    g = sns.lmplot(x="Amount density",
                   y='Average Degree',
                   col='media_new',
                   data=dfcellcon,
                   hue='media_new',
                   col_wrap=2,
                   size=5,
                   scatter_kws={"s": 25, "alpha": 1})
    g.set(ylim=(0, 4.))
    g.fig.suptitle('Mito. amount density vs average degr.', fontsize=22)
    g.fig.subplots_adjust(top=.9)

with sns.plotting_context('talk', font_scale=1.25):
    g = sns.lmplot(x="Amount density",
                   y='Average Degree',
                   col='media_new',
                   hue='media_new',
                   data=dfcellcon,
                   col_wrap=2,
                   size=5,
                   scatter_kws={"s": 25, "alpha": 1})
    g.set(ylim=(1, 3.5))
    g.fig.suptitle(u'Average degree vs mito. amount density', fontsize=22)
    g.fig.subplots_adjust(top=.9)

# subset of volume ratios plots
with sns.plotting_context('talk', font_scale=1.25):
    g = sns.lmplot(x="Amount density",
                   y='Average Degree',
                   col='media_new',
                   data=dfsubset,
                   hue='media_new',
                   col_wrap=2,
                   size=5,
                   scatter_kws={"s": 25, "alpha": 1})
    g.set(ylim=(1, 3.5))
    g.fig.suptitle('Mito. amount density vs average degr. (subset)',
                   fontsize=22)
    g.fig.subplots_adjust(top=.9)

    f = sns.lmplot(x="Amount density",
                   y=u'ΔΨ Unscaled',
                   col='media_new',
                   data=dfsubset,
                   col_wrap=2,
                   hue='media_new',
                   size=5,
                   scatter_kws={"s": 25, "alpha": 1})
    f.fig.suptitle(u'Mito. amount density vs ΔΨ (subset)', fontsize=22)
    f.fig.subplots_adjust(top=.9)
    for ax in f.axes:
        ax.get_lines()[0].set_visible(False)

#  ΔΨ vs tube width
dfwidth = pd.DataFrame({'tube width': df.loc[:, 'mito_tubew'],
                        'media_new': df.loc[:, 'media_new']})
g, axt = plt.subplots(1, 1, figsize=(11, 8.5))
g = sns.violinplot('media_new',
                   'tube width',
                   data=dfwidth, ax=axt)

x = pd.DataFrame({'width_cor_rfp': [i[0] for i in df.mito_widcoef],
                  'pval': [i[1] for i in df.mito_widcoef]},
                 index=df.index)

dftubew = pd.concat([x, df.loc[:, 'media_new']], axis=1)
with sns.plotting_context('talk', font_scale=1.25):
    f, ax0 = plt.subplots(1, 1, figsize=(11, 8.5), sharex=True)
    sns.barplot('media_new',
                'width_cor_rfp',
                data=dftubew.loc[dftubew.pval < .05],
                ax=ax0)
    ax0.set_xlabel('')
    ax0.set_ylabel('R value per cell')
    ax0.set_title('Distribution of R values for '
                  'correlation of cell tubewidth vs cell RFP intensity')

y = pd.DataFrame({'width_cor_dy': [i[0] for i in df.mito_widcoefDY],
                  'pval': [i[1] for i in df.mito_widcoefDY]},
                 index=df.index)

dftubew2 = pd.concat([y, df.loc[:, 'media_new']], axis=1)

with sns.plotting_context('talk', font_scale=1.25):
    f, ax0 = plt.subplots(1, 1, figsize=(11, 8.5), sharex=True)
    sns.violinplot('media_new',
                   'width_cor_dy',
                   data=dftubew2.loc[dftubew.pval < .05],
                   ax=ax0)
    ax0.set_xlabel('')
    ax0.set_ylabel('R value per cell')
    ax0.set_title(u'Distribution of R values for '
                  'correlation of cell tubewidth vs cell DY')


with sns.plotting_context('talk', font_scale=1.25):
    f, ax0 = plt.subplots(1, 1, figsize=(11, 8.5), sharex=True)
    sns.violinplot('media_new',
                   'mito_cell_ave_rfp',
                   data=df,
                   ax=ax0)
    ax0.set_xlabel('')
    ax0.set_ylim(0, df.mito_cell_ave_rfp.quantile([0, 0.95]).tolist()[1])

# =============================================================================
# O2 data
# =============================================================================
o2mito = {'YPR': .0180,
          'YPL': 0.0176,
          'YPE': 0.0113,
          'YPD': 0.0056}
o2cell = {'YPR': .540,
          'YPL': 0.321,
          'YPE': 0.238,
          'YPD': 0.066}
dfcellcon['O2 per mito vol'] = df.media.map(o2mito)
dfcellcon['OCR per cell'] = df.media.map(o2cell)

with sns.plotting_context('talk', font_scale=1.25):
    f, (ax1, ax3, ax2) = plt.subplots(3, 1,
                                      figsize=(11, 8.5))
    sns.barplot('media_new',
                'O2 per mito vol',
                data=dfcellcon[~(dfcellcon['OCR per cell'].isnull())],
                ax=ax1)
    ax1.set_ylabel(u'OCR per mito vol /µm³')
    ax1.set_yticks(np.arange(0.0, 0.022, 0.005))
    ax1.set_xlabel('')
    plt.setp(ax1.get_xticklabels(), visible=False)

    sns.boxplot('media_new',
                u'ΔΨ Unscaled',
                data=dfcellcon[~(dfcellcon['OCR per cell'].isnull())],
                notch=True,
                ax=ax2)
    ax2.set_ylabel(u'ΔΨ Unscaled')
    ax2.set_ylim(0, 1200)
    ax2.set_yticks(np.arange(0, 1500, 500))
    ax2.set_xlabel('')

    sns.boxplot('media_new',
                'Amount density',
                data=dfcellcon[~(dfcellcon['OCR per cell'].isnull())],
                notch=True,
                ax=ax3)
    ax3.set_ylabel('Amount Density')
    ax3.set_yticks(np.arange(0., 1., .2))
    ax3.set_ylim(0.2, .95)
    ax3.set_xlabel('')
    plt.setp(ax3.get_xticklabels(), visible=False)
