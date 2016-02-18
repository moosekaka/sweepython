# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 10:55:07 2015

@author: sweel
"""
import pandas as pd
import cPickle as pickle
import scipy.stats as sp
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from networkHet.mungedata import MungeDataFuncs as md
import statsmodels.formula.api as sm

sns.plotting_context('talk', font_scale=1.4)
sns.set(style="whitegrid")
# pylint: disable=C0103b
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
                      'Surf': dfsize.loc[:, 'Surf'],
                      'mitolen': df.loc[:, 'mito_totlen'],
                      'media': df.loc[:, 'media']})

dfvol['mitovol'] = np.pi * (.15)**2 * dfvol.mitolen
dfvol['Vol Ratio'] = dfvol.mitovol / dfvol.Vol
dfvol['Quasik'] = dfvol.mitolen / dfvol.Surf

g = md.boxviol(dfvol, 'Quasik', 'media')
g.set_title('Quasi Density', fontsize=24)
g.set_ylabel('Quasi Density')
plt.ylim(0)
h = md.boxviol(dfvol, 'Surf', 'media')
h.set_title(r'Cell Surface Area $\mu m^{2}$', fontsize=24)
h.set_ylabel('')
plt.ylim(0)

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
                    dfvol.loc[:, ['Vol Ratio', 'Quasik']],
                    df.loc[:, 'media']], axis=1)

cols = ['mito_bpts_dyraw',
        'mito_btwcntr_uw',
        'mito_knn_uw',
        'mito_clscntr_uw',
        'mito_clstcf_uw',
        'charpl_norm_len']

print("\nstat test for local connectivity vs quasidensity\n" + 79 * '=')
dfl = pd.DataFrame()

data = []
for mem in sorted(pd.unique(dfconn.media)):
    mask = dfconn.media == mem
    dftest = dfconn.loc[mask]
    for mea in cols:
        res = sp.pearsonr(dftest.loc[mask, 'Quasik'],
                          dftest.loc[mask, mea])
        data.append((mem, mea, res[0], res[1]))
dfl = pd.DataFrame(data)
dfl.columns = ['media', 'type', 'r', 'p']
dfl = dfl.pivot(index='type', columns='media', values='r')
dfl.to_clipboard(float_format='%6.4f')
#        print '%-4s:%-15s:%6.3f:%15.13f' % (mem, mea, res[0], res[1])

dfconn['Nearest Neighbor Degree'] = dfconn.mito_knn_uw
with sns.plotting_context('talk', font_scale=1.5):
    g = sns.lmplot(x="Quasik",
                   y='Nearest Neighbor Degree',
                   hue='media',
                   data=dfconn,
                   #                   col_wrap=2,
                   size=5,
                   scatter_kws={"s": 25, "alpha": 1})
    g.set(ylim=(0, 4.5))

#   Topology  (by cell)
dfcellcon = df.ix[:, ['mito_cell_avedyr',
                      'cell_coefvar_r',
                      'mito_beta_geo',
                      'mito_beta_top',
                      'mito_phi',
                      'mito_pk3',
                      'mito_avgdeg',
                      'mito_edgenum',
                      'mito_tubew']]
dfcellcon = pd.concat([dfcellcon,
                       dfvol.loc[:, ['Vol Ratio', 'Quasik']],
                       df.loc[:, 'media']], axis=1)
#==============================================================================
#         OLS fits of conn with surface density with the interaction term test
#==============================================================================


def ferment(t):
    if t == 'YPD':
        return 'ferment'
    else:
        return 'resp'
dfcellcon['ferment'] = dfcellcon.media.apply(lambda x: ferment(x))

result = sm.ols(formula='mito_avgdeg~ Quasik',
                data=dfcellcon[dfcellcon.ferment == 'resp']).fit()
print(result.summary())


dfcellcon = dfcellcon[dfcellcon.mito_cell_avedyr <= 2000]  # exclude hi YPE's
dfcellcon[r'$\Delta \Psi$ Unscaled'] = dfcellcon["mito_cell_avedyr"]

with sns.plotting_context('talk', font_scale=1.5):
    g = sns.lmplot(x="Quasik",
                   y='mito_avgdeg',
                   hue='media',
                   data=dfcellcon,
                   #                   col_wrap=2,
                   size=5,
                   scatter_kws={"s": 25, "alpha": 1})
    g.set(ylim=(0, 3.5))

#==============================================================================
#
#==============================================================================
dfwidth = pd.DataFrame({'tube width': df.loc[:, 'mito_tubew'],
                        'media': df.loc[:, 'media']})
g = md.boxviol(dfwidth, 'tube width', 'media')
writer = pd.ExcelWriter('output_tube.xlsx')
newstring = 'mean_%s' % dfwidth.columns[1]
res = md.multiple_test(dfwidth, 'tube width', 'media')
md.result_to_excel(res, newstring, writer)

with sns.plotting_context('talk', font_scale=1.5):
    g = sns.lmplot(x="Quasik",
                   y=r'$\Delta \Psi$ Unscaled',
                   col="media",
                   hue='media',
                   data=dfcellcon,
                   col_wrap=2,
                   size=5,
                   scatter_kws={"s": 25, "alpha": 1})
    g.set(ylim=(0, 2000))
    g.set(xlim=(0.))
    for i in g.axes:
        i.get_lines()[0].set_visible(False)
#    g.axes[3].get_lines()[0].set_visible(False)

cols = [u'$\Delta \Psi$ Unscaled',
        'Quasik',
        'mito_beta_geo',
        'mito_beta_top',
        'mito_phi',
        'mito_pk3',
        'mito_avgdeg',
        'mito_edgenum',
        'mito_tubew']

print("\nstat test for global connectivity vs quasidensity\n" + 79 * '=')
data1 = []
for mem in sorted(pd.unique(dfcellcon.media)):
    mask = dfcellcon.media == mem
    dftest = dfcellcon.loc[mask]
    for mea in cols:
        res = sp.pearsonr(dftest.loc[mask, 'Quasik'],
                          dftest.loc[mask, mea])
        data1.append((mem, mea, res[0], res[1]))
#        print '%-4s:%-15s:%6.3f:%15.13f' % (mem, mea, res[0], res[1])
dfg = pd.DataFrame(data1)
dfg.columns = ['media', 'type', 'r', 'p']
dfg = dfg.pivot(index='type', columns='media', values='r')
dfg.to_clipboard(float_format='%6.4f')


colors = ["medium green", "pale red", "dusty purple"]
with sns.plotting_context('talk', font_scale=1.3):
    with sns.color_palette(sns.xkcd_palette(colors)):

        g = sns.FacetGrid(
            dfcellcon.ix[
                dfcellcon.media != 'YPD'],
            row='media',
            hue='media')
        g = g.map(sns.distplot, "Quasik", kde=False)

#==============================================================================
# subset of volume ratios
#==============================================================================
dfsubset = dfcellcon.ix[(dfcellcon['Quasik'] > .5) &
                        (dfcellcon['Quasik'] < .9)]
a = dfsubset.groupby('media').count()
b = dfcellcon.groupby('media').count()
c = 1. * a['Quasik'] / b['Quasik']
print "\npercentage of population in subset [%4.2f - %4.2f] \n %s" %\
    (dfsubset.Quasik.min(), dfsubset.Quasik.max(), c)

with sns.plotting_context('talk', font_scale=1.5):
    with sns.color_palette(sns.xkcd_palette(colors)):
        g = sns.lmplot(x="Quasik",
                       y='mito_avgdeg',
                       col="media",
                       hue='media',
                       data=dfsubset,
                       col_wrap=2,
                       size=5,
                       scatter_kws={"s": 25, "alpha": 1})
        g.set(ylim=(1, 3.5))

        f = sns.lmplot(x="Quasik",
                       y='$\Delta \Psi$ Unscaled',
                       col="media",
                       hue='media',
                       data=dfsubset,
                       col_wrap=2,
                       size=5,
                       scatter_kws={"s": 25, "alpha": 1})
        for ax in f.axes:
            ax.get_lines()[0].set_visible(False)

cols = dfsubset.columns[:8]
print"\nstat test for subset [%4.2f - %4.2f] of volume ratios\n" %\
    (dfsubset.Quasik.min(), dfsubset.Quasik.max())
data2 = []
for mea in cols:
    for mem in sorted(pd.unique(dfsubset.media)):
        mask = dfcellcon.media == mem
        res = sp.pearsonr(dfsubset.loc[mask, 'Quasik'],
                          dfsubset.loc[mask, mea])
        data2.append((mem, mea, res[0], res[1]))
#        print '%-4s: %-15s:%6.3f:%15.13f' % (mem, mea, res[0], res[1])
dfs = pd.DataFrame(data2)
dfs.columns = ['media', 'type', 'r', 'p']
dfs = dfs.pivot(index='type', columns='media', values='r')
dfs.to_clipboard(float_format='%6.4f')


dfsubset2 = dfconn.ix[(dfconn['Quasik'] > .5) &
                      (dfconn['Quasik'] < .9)]
cols = dfsubset2.columns[:7]
data3 = []
for mea in cols:
    for mem in sorted(pd.unique(dfsubset2.media)):
        mask = dfconn.media == mem
        res = sp.pearsonr(dfsubset2.loc[mask, 'Quasik'],
                          dfsubset2.loc[mask, mea])
        data3.append((mem, mea, res[0], res[1]))
#        print '%-4s: %-15s:%6.3f:%15.13f' % (mem, mea, res[0], res[1])
dfs2 = pd.DataFrame(data3)
dfs2.columns = ['media', 'type', 'r', 'p']
dfs2 = dfs2.pivot(index='type', columns='media', values='r')
dfs2.to_clipboard(float_format='%6.4f')

#==============================================================================
#  Vol ratio, RFP and DY vs tube width
#==============================================================================
with sns.plotting_context('talk', font_scale=1.5):
    g = sns.lmplot(x="Quasik",
                   y='mito_avgdeg',
                   col="media",
                   hue='media',
                   data=dfcellcon,
                   col_wrap=2,
                   size=5,
                   scatter_kws={"s": 25, "alpha": 1})
    g.set(ylim=(0, 4.))
    g.set_xlabels("Quasi Density")
#    g.axes[0].get_lines()[0].set_visible(False)


x = pd.DataFrame({'width_cor_rfp': [i[0] for i in df.mito_widcoef],
                  'pval': [i[1] for i in df.mito_widcoef]},
                 index=df.index)

dftubew = pd.concat([x, df.loc[:, 'media']], axis=1)
with sns.plotting_context('talk', font_scale=1.5):
    f, ax0 = plt.subplots(1, 1, figsize=(11, 8.5), sharex=True)
    sns.barplot('media',
                'width_cor_rfp',
                data=dftubew.loc[dftubew.pval < .05],
                hue='media',
                ax=ax0)
    ax0.set_xlabel('')
    ax0.set_ylabel('R value per cell')
    ax0.get_legend().set_visible(False)
    ax0.set_title(
        'Distribution of R values for correlation of cell tubewidth vs cell RFP intensity')

y = pd.DataFrame({'width_cor_dy': [i[0] for i in df.mito_widcoefDY],
                  'pval': [i[1] for i in df.mito_widcoefDY]},
                 index=df.index)

dftubew2 = pd.concat([y, df.loc[:, 'media']], axis=1)

with sns.plotting_context('talk', font_scale=1.5):
    f, ax0 = plt.subplots(1, 1, figsize=(11, 8.5), sharex=True)
    sns.violinplot('media',
                   'width_cor_dy',
                   data=dftubew2.loc[dftubew.pval < .05],
                   hue='media',
                   ax=ax0)
    ax0.set_xlabel('')
    ax0.set_ylabel('R value per cell')
    ax0.get_legend().set_visible(False)
    ax0.set_title(
        'Distribution of R values for correlation of cell tubewidth vs cell DY ')


with sns.plotting_context('talk', font_scale=1.5):
    f, ax0 = plt.subplots(1, 1, figsize=(11, 8.5), sharex=True)
    sns.violinplot('media',
                   'mito_cell_ave_rfp',
                   data=df,
                   hue='media',
                   ax=ax0)
    ax0.set_xlabel('')
    ax0.set_ylim(0, 2000)
    ax0.get_legend().set_visible(False)


# =============================================================================
# O2 data here
# =============================================================================
dfcellcon['Number of Edges'] = dfcellcon.mito_edgenum
dfcellcon['Average Degree'] = dfcellcon.mito_avgdeg
dfcellcon['O2 per mito vol'] = ''
dfcellcon['OCR per cell'] = ''

dfcellcon.loc[dfcellcon.media == 'YPD', ['O2 per mito vol']] = .0056
dfcellcon.loc[dfcellcon.media == 'YPE', ['O2 per mito vol']] = .0113
dfcellcon.loc[dfcellcon.media == 'YPL', ['O2 per mito vol']] = .0176
dfcellcon.loc[dfcellcon.media == 'YPR', ['O2 per mito vol']] = .0180
dfcellcon.loc[dfcellcon.media == 'YPD', ['OCR per cell']] = .066
dfcellcon.loc[dfcellcon.media == 'YPE', ['OCR per cell']] = .238
dfcellcon.loc[dfcellcon.media == 'YPL', ['OCR per cell']] = .321
dfcellcon.loc[dfcellcon.media == 'YPR', ['OCR per cell']] = .540


with sns.plotting_context('talk', font_scale=1.5):
    g = sns.lmplot(x="Quasik",
                   y='Average Degree',
                   col="media",
                   hue='media',
                   data=dfcellcon,
                   col_wrap=2,
                   size=5,
                   scatter_kws={"s": 25, "alpha": 1})
    g.set(ylim=(1, 3.5))


with sns.plotting_context('talk', font_scale=1.5):
    f, (ax1, ax0) = plt.subplots(2, 1,
                                 figsize=(11, 8.5),
                                 sharex=True)
    sns.violinplot('media',
                   'Quasik',
                   data=dfcellcon,
                   hue='media',
                   ax=ax0)
#    ax0.set_title(r'OCR per cell$\ (\times\ 10^{7})$ ')
    ax0.set_xlabel('')
    ax0.set_ylabel('Quasi Density')
    ax0.get_legend().set_visible(False)

    sns.barplot('media',
                'O2 per mito vol',
                data=dfcellcon,
                hue='media',
                ax=ax1)
#    ax1.set_title(r'OCR per mito vol /$\mu m^{3}$')
    ax1.set_ylabel(r'OCR per mito vol /$\mu m^{3}$')
    ax1.set_xlabel('')
    ax1.get_legend().set_visible(False)


with sns.plotting_context('talk', font_scale=1.5):
    f, (ax1, ax3, ax2) = plt.subplots(3, 1,
                                      figsize=(11, 8.5),
                                      sharex=True)
    sns.barplot('media',
                'O2 per mito vol',
                data=dfcellcon,
                hue='media',
                ax=ax1)
    ax1.set_ylabel(r'OCR per mito vol /$\mu m^{3}$')
    ax1.set_xlabel('')
    ax1.get_legend().set_visible(False)

    sns.violinplot('media',
                   'mito_cell_avedyr',
                   data=dfcellcon,
                   hue='media',
                   ax=ax2)
    ax2.set_ylabel(r'$\Delta \Psi$ Unscaled')
    ax2.set_ylim(0, 2000)
    ax2.set_xlabel('')
    ax2.get_legend().set_visible(False)

    sns.violinplot('media',
                   'Quasik',
                   data=dfcellcon,
                   hue='media',
                   estimator=np.median,
                   ax=ax3)
    ax3.set_ylabel('Quasi Density')
    ax3.set_ylim(0)
    ax3.set_xlabel('')
    ax3.get_legend().set_visible(False)

    sns.despine(bottom=True)
#    plt.setp(f.axes, yticks=[])
    plt.tight_layout(h_pad=3)
#==============================================================================
#       Misc
#==============================================================================

# need loc method to avoid chain indexing error message
#DF = df.loc[:, ['mito_edge_avedy', 'media']]
# DF[r'mean edge $\Delta \Psi$'] = \
#    DF['mito_edge_avedy'].apply(np.mean)  # easy agg. of edge to cell values
#md.boxviol(DF, r'mean edge $\Delta \Psi$', 'media')
#
#DF = df.loc[:, ['mito_totlen', 'media']]
# DF[r'Mito Cell Vol $\Mu m^3$'] = \
#    DF['mito_totlen'].apply(lambda x: x*np.pi*.15**2)
#md.boxviol(DF, r'Mito Cell Vol $\Mu m^3$', 'media', pltLims=(0, 10))


# =============================================================================
#   DY vs edgelens
# =============================================================================
#df5 = df.loc[:, ['mito_edge_avedyr', 'mito_edgelen']]
#dat = pd.DataFrame()
# for cell in df5.index:
#    dtemp = df5.ix[cell]  # Series of edge data per cell
#    dfedges = pd.DataFrame({i: pd.Series(dtemp[i]) for i
#                            in dtemp.index[:2]})
#    dfedges['media'] = cell[:3]
#    dat = dat.append(dfedges, ignore_index=True)
#
#dat_long = dat[(dat['mito_edgelen'] > 0.) & (dat['mito_edgelen'] < 5.5)]
# gr2 = sns.lmplot('mito_edgelen',
#                 'mito_edge_avedyr',
#                 hue='media',
#                 data=dat_long)
# plt.ylim(0)
# plt.xlim(0)

#dft = convert_longform(df, 'mito_edge_avedy')
#dfl = convert_longform(df, 'mito_edgelen')
#dft['treshold'] = None
#dft['edgelen'] = dfl.loc[:, 'mito_edgelen']
#dft.loc[dft.edgelen < 1.5, ['treshold']] = 'Short'
#dft.loc[dft.edgelen > 3.5, ['treshold']] = 'Long'
# dft.loc[(dft.edgelen > 1.5) & (dft.edgelen < 3.5),
#        ['treshold']] = 'Medium'
# sns.violinplot(x='treshold',
#               y='mito_edge_avedy',
#               hue='cat',
#               data=dft[(dft.cat == 'YPE') |
#                        (dft.cat == 'YPL') |
#                        (dft.cat == 'YPR')],
#               order=['Short', 'Medium', 'Long'],
#               saturation=.75,
#               scale_hue=True)
