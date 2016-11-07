# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 16:10:16 2016
Module for plots of analysis of mother bud function in budding yeast
"""
import os
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
HUE_ODR = [u'WT_YPE', u'ΔMFB1', u'ΔNUM1', u'ΔYPT11']
savefolder = op.expanduser(os.sep.join(['~', 'Dropbox', 'SusanneSweeShared',
                                        'figures_for_mutants']))
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
lowerthresh = 0.25

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
df2 = outputargs['data']
df = pd.concat([df, df2.loc[:, 'bud_dy_var']], axis=1, join='inner')
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
medbud = dfmerge.groupby('media_bud').budvol.quantile([lowerthresh,
                                                       0.5, 0.999])
grp = dfmerge.groupby('media_bud')
meds = pd.concat([grp
                  .get_group(key)
                  [(grp.get_group(key).budvol >= medbud[key, lowerthresh]) &
                      (grp.get_group(key).budvol < medbud[key, 0.999])]
                  for key in grp.groups.keys()])
meds.reset_index(drop=True, inplace=True)
meds = meds.loc[:,
                meds.columns.isin([0.33, 0.67, 1.0, 1.5, 2.0,
                                   u'media_bud', u'binvol_bud'])]
alls = dfmerge.loc[:,
                   dfmerge.columns.isin([0.33, 0.67, 1.0, 1.5, 2.0,
                                        u'media_bud', u'binvol_bud'])]
# =============================================================================
# Data long form
# =============================================================================
frac = pd.melt(df, id_vars=['media_new'],
               value_name=u'ΔΨ scaled',
               value_vars=['frac'])

relabels = {'DY_median_mom': u'ΔΨ Mom', 'DY_median_bud': u'ΔΨ Bud'}
mb_dy = pd.melt(df, id_vars=['media_new'],
                value_name=u'ΔΨ scaled',
                value_vars=mombud_dy_vars)
mb_dy.replace({'variable': relabels}, inplace=True)

size = pd.melt(df, id_vars=['media_new', 'bud_dy_var'],
               value_vars='bud_diameter')
size.rename(columns={'value': 'bud_diameter',
                     'bud_dy_var': u'ΔΨ Bud (var)'},
            inplace=True)
size = size[size['media_new'].isin([u'WT_YPE',
                                    u'ΔMFB1',
                                    u'ΔNUM1',
                                    u'ΔYPT11'])]

outkws1 = dict(default_ylims=[0.05, 0.95], plt_type='boxplot',
               labeller=labelNormal, col_order=COL_ODR)

outkws2 = dict(default_ylims=[0.1, 0.9],
               labeller=labelFacet, col_order=COL_ODR)

outkws3 = dict(default_ylims=[0.05, 0.95], plt_type='boxplot',
               labeller=labelNormal, hue_order=HUE_ODR)

# test of sizes


def vertical_mean_line(x, **kwargs):
    plt.axvline(x.quantile(0.1), color='r', linewidth=1),
    plt.axvline(x.quantile(0.15), linewidth=1),
    plt.axvline(x.quantile(0.2), color='g', linewidth=1)
    plt.axvline(x.quantile(0.25), color='m', linewidth=1)
g = sns.FacetGrid(size, col='media_new', col_wrap=2)
g.map(vertical_mean_line, 'bud_diameter')
g.map(plt.scatter, 'bud_diameter', u'ΔΨ Bud (var)')

# =============================================================================
# mombudDY plot
# =============================================================================
with sns.plotting_context('talk', font_scale=1.5):
    plt.rcParams['figure.figsize'] = (16, 11)
    set1 = dict(x='media_new', y=u'ΔΨ scaled',
                hue='variable', group_key='media_new',
                title=u'ΔΨ mother vs bud', ylim=(0, 1.), notch=True,
                )

    plv1 = plviol(**outkws1)
    plv1.plt(data=mb_dy, **set1)
    plv1.ax.xaxis.label.set_visible(False)
    plv1.ax.get_legend().set_title('')
    plv1.save_figure(op.join(savefolder, 'boxplot mombud dy.png'))

# =============================================================================
# fracDY plot
# =============================================================================

with sns.plotting_context('talk', font_scale=1.5):
    plt.rcParams['figure.figsize'] = (16, 11)
    set2v = dict(x='media_new', y=u'ΔΨ scaled',
                 group_key='media_new', inner=None, notch=True,
                 title=u'',
                 ylim=[0, 2.5])

    plv2 = plviol(**outkws1)
    plv2.plt(data=frac, **set2v)
    plv2.ax.axhline(1.0)
    plv2.ax.xaxis.label.set_visible(False)
    plv2.turn_off_legend()
    plv2.save_figure(op.join(savefolder, 'boxplot frac dy.png'))

# =============================================================================
# mom, bud,  cell axis long form DF
# =============================================================================
momdy = pd.melt(dfmom,
                id_vars=['media', 'binvol'],
                var_name='mom axis position',
                value_name=u'ΔΨ scaled',
                value_vars=outputargs['mbax'][1:])

buddy = pd.melt(dfbud,
                id_vars=['media', 'binvol'],
                var_name='bud axis position',
                value_name=u'ΔΨ scaled',
                value_vars=outputargs2['mbax'].tolist()[1:])

buddy_meds = pd.melt(meds,
                     id_vars=['media_bud', 'binvol_bud'],
                     var_name='cell axis position',
                     value_name=u'ΔΨ scaled',
                     value_vars=meds.columns[:5].values.tolist())

allsizes = pd.melt(alls,
                   id_vars=['media_bud', 'binvol_bud'],
                   var_name='cell axis position',
                   value_name=u'ΔΨ scaled',
                   value_vars=meds.columns[:5].values.tolist())
# =============================================================================
# PLOTS FACETTED
# =============================================================================
with sns.plotting_context('talk', font_scale=.95):
    with sns.color_palette('colorblind'):
        plt.rcParams['figure.figsize'] = (16, 11)

        set7 = dict(col_wrap=3, col='media_bud', hue='media_bud',
                    hue_order=HUE_ODR,
                    sharex=True, sharey=True, col_order=COL_ODR[:3],
                    ylim=(0.0, 1.0))

        plv7 = plfacet(plt_type='pointplot', **outkws2)
        plv7.plt(data=buddy_meds,
                 mapargs=['cell axis position', u'ΔΨ scaled'],
                 **set7)
        for i in plv7.facet_obj.axes:
            wt = sns.pointplot('cell axis position', u'ΔΨ scaled',
                               data=buddy_meds[buddy_meds.media_bud == u'WT_YPE'],
                               ax=i, markers='x',)
            [j.set_alpha(.75) for j in wt.axes.collections]
            [j.set_alpha(.75) for j in wt.axes.lines]

        plv7.save_figure(op.join(savefolder, 'medbuds.png'))

        set8 = dict(col_wrap=3, col='media_bud', hue='media_bud',
                    hue_order=HUE_ODR,
                    sharex=True, sharey=True, col_order=COL_ODR[:3],
                    ylim=(0.0, 1.0),
                    )

        plv8 = plfacet(plt_type='pointplot', **outkws2)
        plv8.plt(data=allsizes,
                 mapargs=['cell axis position', u'ΔΨ scaled'],
                 **set8)
        for i in plv8.facet_obj.axes:
            wt = sns.pointplot('cell axis position', u'ΔΨ scaled',
                               data=allsizes[allsizes.media_bud == u'WT_YPE'],
                               ax=i, markers='x')
            [j.set_alpha(.75) for j in wt.axes.collections]
            [j.set_alpha(.75) for j in wt.axes.lines]
        plv8.save_figure(op.join(savefolder, 'allsizes ptplt.png'))

# BOX PLOTS VERSION
with sns.plotting_context('talk', font_scale=1.1):

    g2 = sns.factorplot('cell axis position',
                        u'ΔΨ scaled',
                        data=buddy_meds,
                        kind='box', col='media_bud',
                        col_order=COL_ODR, notch=True, col_wrap=3)
    labelFacet(g2, mbfuncs.get_group_counts(g2.data))
    g2.set(ylim=tuple([0, 1.25]))
    plt.savefig(op.join(savefolder, 'box med sized buds.png'))

    g3 = sns.factorplot('cell axis position',
                        u'ΔΨ scaled',
                        data=allsizes,
                        kind='box', col='media_bud',
                        col_order=COL_ODR, notch=True, col_wrap=3)
    labelFacet(g3, mbfuncs.get_group_counts(g3.data))
    g3.set(ylim=tuple([0, 1.25]))
    plt.savefig(op.join(savefolder, 'box all bud sizes.png'))

with sns.plotting_context('talk', font_scale=1.25):
    with sns.color_palette('colorblind'):
        plt.rcParams['figure.figsize'] = (16, 11)
        set10 = dict(x='cell axis position', y=u'ΔΨ scaled',
                     hue='media_bud', hue_order=HUE_ODR,
                     title=(u'Average population ΔΨ along cell axis '
                            '(cells with all sizes of buds)'),
                     group_key='media_bud', notch=True, ylim=(0.0, 1.25))

        plv10 = plviol(**outkws3)
        plv10.plt(data=allsizes, **set10)
        plv10.ax.legend(loc=2, title='')
        plt.savefig(op.join(savefolder, 'box all onerow.png'))

        set11 = dict(x='cell axis position', y=u'ΔΨ scaled',
                     hue='media_bud', hue_order=HUE_ODR,
                     title=(u'Average population ΔΨ along cell axis '
                            '(cells with {:2.0f}th percentile bud diameters '
                            'and up)'
                            ).format(lowerthresh*100),
                     notch=True, ylim=(0.0, 1.25))
        plv11 = plviol(**outkws3)
        plv11.plt(data=buddy_meds, **set11)
        plv11.ax.legend(loc=2, title='')
        plt.savefig(op.join(savefolder,
                            ('box {:2.0f}th and up buds onerow.png')
                            .format(lowerthresh*100)))
