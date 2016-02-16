# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 10:55:07 2015
Data analysis  for local connectivity and function (including an boostrap comp
for branchpoints vs fake bpoints)
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

dflists = df.loc[:, ['mito_bpts_dy',
                     'mito_bpts_dyraw',
                     'mito_bptcoefvar_raw',
                     'mito_btwcntr_uw',
                     'mito_btwcntr_w',
                     'mito_clscntr_uw',
                     'mito_clscntr_w',
                     'mito_clstcf_uw',
                     'mito_clstcf_w',
                     'mito_edge_avedy',
                     'mito_edge_avedyr',
                     'mito_edge_coefvar',
                     'mito_edge_coefvarr',
                     'mito_edge_stddy',
                     'mito_edge_stddyr',
                     'mito_edgelen',
                     'mito_knn_uw',
                     'mito_knn_w',
                     'mito_bootbpts_dyraw',
                     'lags_1',
                     'media']]

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

# =============================================================================
#   Subdatasets
# =============================================================================
#   Network Connectivity (by branchpoints)
dfconn = df.ix[:, ['mito_bpts_dyraw',
                   'mito_btwcntr_uw',
                   'mito_knn_uw',
                   'mito_clscntr_uw',
                   'mito_clstcf_uw']]
frames = []
for i in dfconn.columns:
    temp = md.convert_longform(dfconn, i)
    temp.index = temp.cell
    temp = temp.drop('cell', 1)
    frames.append(temp)
BIG = pd.concat(frames, axis=1)
BIG = BIG.groupby(BIG.index).mean()

dfconn = pd.concat([BIG,
                    df.loc[:, 'charpl_norm_len'],
                    dfvol.loc[:, 'Vol Ratio'],
                    df.loc[:, 'media']], axis=1)

cols = ['mito_bpts_dyraw',
        'mito_btwcntr_uw',
        'mito_knn_uw',
        'mito_clscntr_uw',
        'mito_clstcf_uw',
        'charpl_norm_len']

dfconn['Nearest Neighbor Degree'] = dfconn.mito_knn_uw

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

cols = ['mito_beta_geo',
        'mito_beta_top',
        'mito_phi',
        'mito_pk3',
        'mito_avgdeg',
        'mito_edgenum']

#==============================================================================
# local conn corr with dyraw
#==============================================================================
pdlists = dflists.ix[:, ['mito_bpts_dyraw',
                         'mito_knn_uw',
                         'mito_btwcntr_uw',
                         'mito_clscntr_uw',
                         'mito_clstcf_uw',
                         'mito_bptcoefvar_raw',
                         'mito_bootbpts_dyraw']]

bpts_dyraw = pd.DataFrame({i: pd.Series(pdlists.ix[i, 0]) for i
                           in pdlists.index})

bpts_knn = pd.DataFrame({i: pd.Series(pdlists.ix[i, 1]) for i
                         in pdlists.index})
bpts_btw = pd.DataFrame({i: pd.Series(pdlists.ix[i, 2]) for i
                         in pdlists.index})
bpts_clsns = pd.DataFrame({i: pd.Series(pdlists.ix[i, 3]) for i
                           in pdlists.index})
bpts_clust = pd.DataFrame({i: pd.Series(pdlists.ix[i, 4]) for i
                           in pdlists.index})
bpts_coefvar = pd.DataFrame({i: pd.Series(pdlists.ix[i, 5]) for i
                             in pdlists.index})
boot_dyraw = pd.DataFrame({i: pd.Series(pdlists.ix[i, 6]) for i
                           in pdlists.index})

# =============================================================================
# giant cell analysis for local conn measures
# =============================================================================
bpts_dyrawG = bpts_dyraw.stack().reset_index()
bpts_knnG = bpts_knn.stack().reset_index()
bpts_btwG = bpts_btw.stack().reset_index()
bpts_clsnsG = bpts_clsns.stack().reset_index()
bpts_clustG = bpts_clust.stack().reset_index()


giant = pd.DataFrame({r'$\Delta \Psi$': bpts_dyrawG.ix[:, 2],
                      'near neigh deg': bpts_knnG.ix[:, 2],
                      'btw cent.': bpts_btwG.ix[:, 2],
                      'close cent.': bpts_clsnsG.ix[:, 2],
                      'clustering': bpts_clustG.ix[:, 2],
                      'zcell': bpts_dyrawG.ix[:, 1]},
                       index= bpts_dyrawG.index)
giant['media'] = giant['zcell'].apply(lambda x: x[:3])
giant = giant.ix[giant[r'$\Delta \Psi$'] <= 2000]
with sns.plotting_context('talk', font_scale=1.3):
    g= sns.lmplot(x='close cent.',
               y=r'$\Delta \Psi$',
               hue='media',
               data=giant,
               col_wrap=2,
               col='media')
    g.axes[0].get_lines()[0].set_visible(False)

cols =giant.columns[1:-2]
#for mem in sorted(pd.unique(giant.media)):
print '\nGiant cell local connectivity vs $\Delta \Psi$ correlations\n' + 79 * '-'
for mea in cols:
    for mem in sorted(pd.unique(giant.media)):
        mask = giant.media == mem
        dftest = giant.loc[mask]
        res = sp.pearsonr(dftest.loc[mask, '$\Delta \Psi$'],
                              dftest.loc[mask, mea])
        print '%-15s:%-5s:%6.3f:%8.4f' % (mea, mem, res[0], res[1])


def corrsp(df1, df2):
    """correlation of every column with every column using sp pearssonr
    """
    out = {}
    for cell in df1.columns:
        x1 = df1.ix[:, cell].dropna()
        y1 = df2.ix[:, cell].dropna()
        out[cell] = sp.pearsonr(x1, y1)
    output = pd.DataFrame({'Rval': [v[0] for v in out.values()],
                           'pval': [v[1] for v in out.values()]},
                          index=[k for k in out.keys()])
    return output


def cordyloc(df1, df2, label, axes):
    """corr local conn with bpts DY
    """
    cordf = corrsp(df1, df2)
    cordf = pd.concat([cordf,
                       df.loc[:, 'media']], axis=1)
    cordf.columns = [label, 'pval', 'type']
    cordf.reset_index(drop=True, inplace=True)
    tot = cordf.groupby('type').count().pval
    sigprop = {}
    for mem in cordf.type.unique():
        sigprop[mem] = 1. * cordf.ix[(cordf.type == mem) &
                                     (cordf.pval < .05)].pval.count()/tot[mem]
        print '%-18s pror. sig. : %5.3f %6s' % (label, sigprop[mem], mem)

    ghand = sns.violinplot(x='type',
                           y=label,
                           data=cordf.ix[cordf.pval < .05],
                           hue='type',
                           ax=axes)
    sns.stripplot(x='type',
                  y=label,
                  data=cordf.ix[cordf.pval < .05],
                  hue='type',
                  ax=axes,
                  size=3,
                  jitter=.15,
                  alpha=.5)

    ghand.set_ylim(-1, 1)
    ghand.set_xlabel('')
    ghand.get_legend().set_visible(False)
    sns.despine(bottom=True)
    return cordf

with sns.plotting_context('talk', font_scale=1.25):
    f, (ax0, ax1, ax2, ax3) = plt.subplots(4, 1,
                                           figsize=(11, 8.5),
                                           sharex=True)
    print '\nbranchpoints vs local connectivity\n' + 79 * '-'
    cordyloc(bpts_dyraw, bpts_knn, 'nearest neigh deg', ax0)
    cordyloc(bpts_dyraw, bpts_btw, 'btwn centr', ax1)
    cordyloc(bpts_dyraw, bpts_clsns, 'closeness centr', ax2)
    cordyloc(bpts_dyraw, bpts_clust, 'cluster coeff', ax3)
f.suptitle(
    r'Distributions of local connectivity correlations with branchpoints $\Delta \Psi$ Raw',
    fontsize=22)

with sns.plotting_context('talk', font_scale=1.25):
    f, (ax0, ax1, ax2, ax3) = plt.subplots(4, 1,
                                           figsize=(11, 8.5),
                                           sharex=True)
    print '\ncoef of var vs local connectivity\n' + 79 * '-'
    cordyloc(bpts_coefvar, bpts_knn, 'nearest neigh deg', ax0)
    cordyloc(bpts_coefvar, bpts_btw, 'btwn centr', ax1)
    cordyloc(bpts_coefvar, bpts_clsns, 'closeness centr', ax2)
    cordyloc(bpts_coefvar, bpts_clust, 'cluster coeff', ax3)


f.suptitle(
    r"Distributions of local connectivity correlations with branchpoints\
    $\Delta \Psi$ Coef of Var Raw",
    fontsize=22)

# =============================================================================
# bpoints bootstrap analy.
# =============================================================================
bootmean = boot_dyraw.apply(np.mean)
bptsmean = bpts_dyraw.apply(np.mean)
df = pd.concat([bptsmean, bootmean], axis=1)
x = pd.Series([i[:3] for i in bootmean.index])
df.columns = ['real', 'bootstrap.']
df.reset_index(drop=True, inplace=True)
df['type'] = x

dfm = pd.melt(df,
              id_vars=['type'],
              value_vars=['real', 'bootstrap.'],
              value_name=r'$\Delta \Psi$ raw')

with sns.plotting_context('talk', font_scale=1.25):
    fig2, ax2 = plt.subplots(1, 1)
    sns.stripplot(x='type',
                  y=r'$\Delta \Psi$ raw',
                  hue='variable',
                  data=dfm,
                  split=True,
                  jitter=.1,
                  size=5,
                  alpha=.5,
                  ax=ax2)
    sns.violinplot(x='type',
                   y=r'$\Delta \Psi$ raw',
                   hue='variable',
                   data=dfm,
                   split=True,
                   inner="quart",
                   ax=ax2)
    ax2.set_ylim(0, 2000)
    ax2.get_legend().set_visible(False)

    h, l = ax2.get_legend_handles_labels()
    h = h[:2]
    l = l[:2]
    ax2.legend(h, l, loc=2)
    fig2.suptitle(r'Distributions of branchpoints vs bootstrapped branchpoints mean $\Delta \Psi$',
                  fontsize=22)
