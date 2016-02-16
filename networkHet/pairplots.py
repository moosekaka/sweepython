# -*- coding: utf-8 -*-
"""
Created on Mon Aug 03 23:14:14 2015
Pair plots
@author: sweel
"""
import pandas as pd
import cPickle as pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mungedata import MungeDataFuncs as md
sns.plotting_context('talk', font_scale=1.4)
# pylint: disable=C0330
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
dfvol = pd.DataFrame({'Vol': dfsize.loc[:, 'Vol'],
                      'mitolen': df.loc[:, 'mito_totlen'],
                      'media': df.loc[:, 'media']})

dfvol['mitovol'] = np.pi * (.15)**2 * dfvol.mitolen
dfvol['Vol Ratio'] = dfvol.mitovol / dfvol.Vol

# =============================================================================
#   Subdatasets
# =============================================================================
#   Connectivity (by branchpoints)
dfconn = df.ix[:, ['mito_bpts_dy',
                   'mito_bpts_dyraw',
                   'mito_btwcntr_uw',
                   'mito_knn_uw',
                   'mito_clscntr_uw',
                   'mito_clstcf_uw',
                   'mito_clstcf_w']]
dfconn = pd.concat([dfconn,
                    dfvol.loc[:, 'Vol Ratio'],
                    df.loc[:, 'media']], axis=1)

#   Topology
ave_edgelen = df.ix[:, 'mito_edgelen'].apply(np.mean)
ave_lag = df.ix[:, 'lags_1'].apply(np.mean)
dftopo = df.ix[:, ['mito_cell_avedy',
                   'mito_cell_avedyr',
                   'cell_coefvar',
                   'cell_coefvar_r',
                   'lags_1',
                   'mito_avgdeg',
                   'mito_totlen',
                   'mito_edgenum',
                   'charpl_norm_len']]

dftopo['ave_edgelen'] = ave_edgelen
dftopo = pd.concat([dftopo,
                    dfvol.loc[:, 'Vol Ratio'],
                    df.loc[:, 'media']], axis=1)

#   Connectivity (by cell)
dfcellcon = df.ix[:, ['mito_cell_avedy',
                      'mito_cell_avedyr',
                      'cell_coefvar',
                      'cell_coefvar_r',
                      'mito_beta_geo',
                      'mito_beta_top',
                      'mito_phi',
                      'mito_pk3']]
dfcellcon = pd.concat([dfcellcon,
                       dfvol.loc[:, 'Vol Ratio'],
                       df.loc[:, 'media']], axis=1)

# =============================================================================
    #  for one condition all pairplots on  bpts connectivity
# =============================================================================
frames = []
for i in dfconn.columns[:-2]:
    temp = md.convert_longform(dfconn, i)
    temp.index = temp.cell
    temp = temp.drop('cell', 1)
    frames.append(temp)
BIG = pd.concat(frames, axis=1)
BIG['cat'] = ''
BIG.ix[:, 'cat'] = BIG.index.map(lambda x: x[:3])
filt = BIG.loc[:, 'mito_knn_uw'] > 3.334  # boolean indexer!
BIG.loc[filt, 'mito_knn_uw'] = BIG.loc[filt, 'mito_knn_uw']\
                               .apply(lambda x: min(x, 3.667))

for mem in sorted(pd.unique(dfconn['media'].values)):
    bpts_conn = BIG.loc[BIG.cat == mem]
    pg1 = sns.PairGrid(bpts_conn)
    pg1 = pg1.map_upper(sns.regplot,
                        scatter_kws={"s": 5},
                        line_kws={'lw': 1.5})
    pg1 = pg1.map_diag(sns.kdeplot, c='r')
    pg1 = pg1.map_lower(sns.regplot,
                        color='g',
                        scatter_kws={"s": 5},

                        line_kws={'lw': 1.5})
    md.adjust_axes(pg1, 5)
    plt.savefig('%s_bptsconn_pair.png' % mem)

# =============================================================================
    #  for one condition all pairplots on  topology
# =============================================================================

for mem in sorted(pd.unique(dftopo['media'].values)):
    cell_top = dftopo.loc[dftopo.media == mem]
    pg1 = sns.PairGrid(cell_top)
    pg1 = pg1.map_upper(sns.regplot,
                        scatter_kws={"s": 5},
                        line_kws={'lw': 1.5})
    pg1 = pg1.map_diag(sns.kdeplot, c='r')
    pg1 = pg1.map_lower(sns.regplot,
                        color='g',
                        scatter_kws={"s": 5},
                        line_kws={'lw': 1.5})
    md.adjust_axes(pg1, 4)
    plt.savefig('%s_topo_pair.png' % mem)

# =============================================================================
    #  for one condition all pairplots on cell connectivity
# =============================================================================

for mem in sorted(pd.unique(dfcellcon['media'].values)):
    cell_con = dfcellcon.loc[dfcellcon.media == mem]
    pg1 = sns.PairGrid(cell_con)
    pg1 = pg1.map_upper(sns.regplot,
                        scatter_kws={"s": 5},
                        line_kws={'lw': 1.5})
    pg1 = pg1.map_diag(sns.kdeplot, c='r')
    pg1 = pg1.map_lower(sns.regplot,
                        color='g',
                        scatter_kws={"s": 5},
                        line_kws={'lw': 1.5})
    md.adjust_axes(pg1, 4)
    plt.savefig('%s_cellconn_pair.png' % mem)
