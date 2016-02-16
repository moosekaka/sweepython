# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 12:07:17 2015
Stat test writeout and graphs
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
                     'mito_cell_ave_gfp']]

dfvol = pd.DataFrame({'Vol': dfsize.loc[:, 'Vol'],
                      'mitolen': df.loc[:, 'mito_totlen'],
                      'media': df.loc[:, 'media']})

dfvol['mitovol'] = np.pi * (.15)**2 * dfvol.mitolen
dfvol['Vol Ratio'] = dfvol.mitovol / dfvol.Vol

g = md.boxviol(dfvol, 'Vol Ratio', 'media')
g.set_title('Volume Ratio', fontsize=24)
g.set_ylabel('')
plt.ylim(0)
h = md.boxviol(dfvol, 'Vol', 'media')
h.set_title(r'Cell Volume $\mu m^{3}$', fontsize=24)
h.set_ylabel('')
plt.ylim(0)

dfscals = pd.concat([dfscals,
                     dfvol.loc[:, 'Vol Ratio'],
                     df.loc[:, 'media']], axis=1)
# =============================================================================
#    #  Stats test (multiple test holms sidak correctiong)
# =============================================================================

#writer = pd.ExcelWriter('output_lists.xlsx')
#for col in dflists.columns[:-1]:
#    DF = df.loc[:, [col, 'media']]
#    newstring = 'mean_edge_%s' % col
#
#    DF[newstring] = DF[col].apply(np.mean)
#    res = md.multiple_test(DF, newstring, 'media')
#    md.result_to_excel(res, newstring, writer)
#writer.save()
writer2 = pd.ExcelWriter('output_scalars.xlsx')
for col in dfscals.columns[:-2]:
    DF = dfscals.loc[:, [col, 'media']]
    newstring = '%s' % col

    res = md.multiple_test(DF, col, 'media')
    md.result_to_excel(res, newstring, writer2)
#writer2.save()
##   violinboxplots
#for col in dflists.columns[:-1]:
#    DF = df.loc[:, [col, 'media']]
#    DF[col] = DF[col].apply(np.mean)
#    axg = md.boxviol(DF, col, 'media')
#    axg.set_title(dic[col], fontsize=20)
#    axg.set_ylabel('')
#    plt.savefig('_%s.png' % col)
#    plt.close()
#for col in dfscals.columns[:-1]:
#    DF = df.loc[:, [col, 'media']]
#    axg = md.boxviol(DF, col, 'media')
#    axg.set_title(dic[col], fontsize=20)
#    axg.set_ylabel('')
#    plt.savefig('_%s.png' % col)
#    plt.close()
#

