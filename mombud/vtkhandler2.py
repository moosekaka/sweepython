# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 16:10:16 2016
Module for plots of analysis of mother bud function in budding yeast
"""
import os.path as op
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import mombud.mungedf as munge
import mombud.functions.vtk_mbplots as mbfuncs
# pylint: disable=C0103
sns.set_style('whitegrid')
plt.rcParams['font.family'] = 'DejaVu Sans'
labelhandler, plviol, plbox, plfacet = (mbfuncs.labelhandler,
                                        mbfuncs.plviol,
                                        mbfuncs.plbox,
                                        mbfuncs.plfacet)

COL_ODR = [u'ΔMFB1', u'ΔNUM1', u'ΔYPT11', u'WT_YPE', u'WT_YPL', u'WT_YPR', ]
COL_ODR_F = ['mfb1', 'num1', 'ypt11', u'WT_YPE', u'WT_YPL', u'WT_YPR', ]
savefolder = r"C:\Users\sweel_Rafelski\Dropbox\SusanneSweeShared\aftermeet"
mombud_dy_vars = ['DY_median_mom', 'DY_median_bud']

plt.close('all')
labelFacet = labelhandler('facet')
labelNormal = labelhandler()
labelRowFacet = labelhandler('rowfacet')

with open('df_mombud_filtered.pkl', 'rb') as inpt:
    df = pickle.load(inpt)

newlabels = {'YPR': u'WT_YPR', 'YPL': u'WT_YPL', 'YPE': u'WT_YPE',
             'WT': u'WT_YPE', 'NUM1': u'ΔNUM1', 'MFB1': u'ΔMFB1',
             'YPT11': u'ΔYPT11'}
df['media_new'] = df.media.map(newlabels)

# =============================================================================
# Run postprocess_df(), which calls inputdata(), mungedata()
# =============================================================================
# #MOM segment into three parts (mbax)
inp_args = {'inpdatpath': 'celldata.pkl',
            # pos. mom/bud cell
            'mbax': [round(i, 2) for i in np.linspace(0, 1, 4)],
            'cellax': np.linspace(0, 1., 11),  # whole cell pos.
            'binsvolbud': np.linspace(0, 40, 5),  # vol bins for bud
            'binsvolmom': np.array([0, 30, 40, 80.]),
            }
# BUD segment into two parts (mbax)
inp_args2 = {'inpdatpath': 'celldata.pkl',
             'mbax': np.linspace(0., 1., 3),  # pos. mom/bud cell
             'cellax': np.linspace(0, 1., 11),  # whole cell pos.
             'binsvolbud': np.linspace(0, 40, 5),  # vol bins for bud
             'binsvolmom': np.array([0, 30, 40, 80.]),
             }
# Run postprocess using two different binnings for mom and buds
outputargs = munge.postprocess_df(**inp_args)
outputargs2 = munge.postprocess_df(**inp_args2)

# =============================================================================
# dfmom and dfbud Dataframes, filtered using keys from df
# =============================================================================
dfmom = outputargs.get('dfmom')
dfmom = dfmom[dfmom.index.isin(df.index)].reset_index()
dfbud = outputargs2.get('dfbud')
dfbud = dfbud[dfbud.index.isin(df.index)].reset_index()

dfmom.replace({'media': newlabels}, inplace=True)
dfbud.replace({'media': newlabels}, inplace=True)

# ============================================================================
# merge dfmom and dfbud, pick medium buds
# =============================================================================
dfmerge = dfmom.merge(dfbud, on='name', suffixes=['_mom', '_bud'])
dfmerge.rename(columns={0.5: 1.5,
                        u'1.0_bud': 2.0,
                        u'1.0_mom': 1.0},
               inplace=True)
medbud = dfmerge.groupby('media_bud').budvol.quantile([0.35, 0.5, 0.65])
grp = dfmerge.groupby('media_bud')
meds = pd.concat([grp
                  .get_group(key)
                  [(grp.get_group(key).budvol > medbud[key, 0.35]) &
                      (grp.get_group(key).budvol < medbud[key, 0.65])]
                  for key in grp.groups.keys()])
meds.reset_index(drop=True, inplace=True)
meds = meds.loc[:,
                meds.columns.isin([0.33, 0.67, 1.0, 1.5, 2.0,
                                   u'media_bud', u'binvol_bud'])]
# =============================================================================
# Data long form
# =============================================================================
frac = pd.melt(df, id_vars=['media_new'],
               value_vars=['frac'])

mb_dy = pd.melt(df, id_vars=['media_new'],
                value_vars=mombud_dy_vars)

size = pd.melt(df, id_vars=['media_new'],
               value_vars='budvol')

outkws1 = dict(default_ylims=[0.05, 0.95],
               labeller=labelNormal, col_order=COL_ODR)

outkws2 = dict(default_ylims=[0.15, 0.9],
               labeller=labelFacet, col_order=COL_ODR_F)
# =============================================================================
# mombudDY plot
# =============================================================================
with sns.plotting_context('talk', font_scale=1.2):
    _, ax0 = plt.subplots(figsize=(20, 16))
    set0 = dict(
        x='media_new',
        y='value',
        hue='variable',
        data=mb_dy,
        width=.75,
        order=COL_ODR,
        notch=True,
        bootstrap=10000,
        medianprops={
            'c': '#ffb16d',
            'markeredgewidth': 2})
    g = sns.boxplot(ax=ax0, **set0)
    plt.savefig(op.join(savefolder, 'boxplot_mombud.png'))

set1 = dict(x='media_new', y='value',
            hue='variable', group_key='media_new',
            title='mombud', ylim=(0, 1.),
            )

plv1 = plviol(**outkws1)
plv1.plt(data=mb_dy, **set1)
plv1.turn_off_legend()
plv1.save_figure(op.join(savefolder, 'violin_mombud.png'))

# =============================================================================
# fracDY plot
# =============================================================================
with sns.plotting_context('talk', font_scale=1.2):
    _, ax1 = plt.subplots(figsize=(20, 16))
    set1 = dict(x='media_new', y='value', data=frac, width=.75,
                order=COL_ODR, join=False, n_boot=1000)
    g = sns.pointplot(ax=ax1, estimator=np.median, **set1)
    g.set(ylim=(0, 2.5))
    plt.savefig(op.join(savefolder, 'CI.png'))

    _, ax2 = plt.subplots(figsize=(20, 16))
    set2 = dict(x='media_new', y='value', data=frac, width=.75,
                order=COL_ODR, notch=True, bootstrap=20000,
                medianprops={'c': '#ed0dd9', 'markeredgewidth': 2})
    g = sns.boxplot(ax=ax2, **set2)
    g.set(ylim=(0, 2.5))
    g.axhline(1.0)
    plt.savefig(op.join(savefolder, 'fracDY.png'))

with sns.plotting_context('talk', font_scale=1.15):
    plt.rcParams['figure.figsize'] = (20, 16)
    set2v = dict(x='media_new', y='value', size=(20, 16),
                 group_key='media_new', inner=None,
                 title='fracDY', ylim='auto')
    plv2 = plviol(**outkws1)
    plv2.plt(data=frac, **set2v)
    plv2.ax.axhline(1.0)
    plv2.turn_off_legend()
    plv2.save_figure(op.join(savefolder, 'test.png'))

# =============================================================================
# mom, bud,  cell axis DY plots
# seaborn map unable to handle unicode in values, so map back greeks
# =============================================================================
labels = {u'ΔNUM1': u'num1', u'ΔMFB1': u'mfb1', u'ΔYPT11': u'ypt11'}
momdy = pd.melt(dfmom,
                id_vars=['media', 'binvol'],
                var_name='mom axis position',
                value_name=u'ΔΨ scaled',
                value_vars=outputargs['mbax'][1:])
momdy.replace({'media': labels}, inplace=True)

buddy = pd.melt(dfbud,
                id_vars=['media', 'binvol'],
                var_name='bud axis position',
                value_name=u'ΔΨ scaled',
                value_vars=outputargs2['mbax'].tolist()[1:])
buddy.replace({'media': labels}, inplace=True)

buddy_meds = pd.melt(meds,
                     id_vars=['media_bud', 'binvol_bud'],
                     var_name='cell axis position',
                     value_name=u'ΔΨ scaled',
                     value_vars=meds.columns[:5].values.tolist())
buddy_meds.replace({'media_bud': labels}, inplace=True)

# =============================================================================
# PLOTS FACETTED
# =============================================================================
set5 = dict(col_wrap=3, col='media', hue='media',
            sharex=True, sharey=True, col_order=COL_ODR_F,
            ylim=(0.0, 1.))

plv5 = plfacet(plt_type='pointplot', **outkws2)
plv5.plt(data=momdy,
         mapargs=['mom axis position', u'ΔΨ scaled'],
         **set5)
plv5.save_figure(op.join(savefolder, 'momDY.png'))

set6 = dict(col_wrap=3, col='media', hue='media',
            sharex=True, sharey=True, col_order=COL_ODR_F,
            ylim=(0.0, 1.0),
            )

plv6 = plfacet(plt_type='pointplot', **outkws2)
plv6.plt(data=buddy,
         mapargs=['bud axis position', u'ΔΨ scaled'],
         **set6)
plv6.save_figure(op.join(savefolder, 'budDY.png'))

set7 = dict(col_wrap=3, col='media_bud', hue='media_bud',
            sharex=True, sharey=True, col_order=COL_ODR_F,
            ylim=(0.0, 1.0),
            )

plv7 = plfacet(plt_type='pointplot', **outkws2)
plv7.plt(data=buddy_meds,
         mapargs=['cell axis position', u'ΔΨ scaled'],
         **set7)
plv7.save_figure(op.join(savefolder, 'medbuds.png'))


#BOX PLOTS VERSION
g0 = sns.factorplot('mom axis position',
                   u'ΔΨ scaled',
                   data=momdy,
                   kind='box', col='media',
                   col_order=COL_ODR_F, notch=True, col_wrap=3)
g0.set(ylim=tuple([0, 1.0]))
plt.savefig(op.join(savefolder, 'box mom dy.png'))

g1 = sns.factorplot('bud axis position',
                   u'ΔΨ scaled',
                   data=buddy,
                   kind='box', col='media',
                   col_order=COL_ODR_F, notch=True, col_wrap=3)
g1.set(ylim=tuple([0, 1.0]))
plt.savefig(op.join(savefolder, 'box bud dy.png'))

g2 = sns.factorplot('cell axis position',
                   u'ΔΨ scaled',
                   data=buddy_meds,
                   kind='box', col='media_bud',
                   col_order=COL_ODR_F, notch=True, col_wrap=3)
g2.set(ylim=tuple([0, 1.0]))
plt.savefig(op.join(savefolder, 'box medbuds.png'))
