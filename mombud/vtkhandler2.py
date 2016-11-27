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
from scipy import stats as ss
import mombud.mungedf as munge
import mombud.functions.vtk_mbplots as mbfuncs
# pylint: disable=C0103
sns.set_style('whitegrid')
plt.rcParams['font.family'] = 'DejaVu Sans'
labelhandler, plviol, plbox, plfacet = (mbfuncs.labelhandler,
                                        mbfuncs.plviol,
                                        mbfuncs.plbox,
                                        mbfuncs.plfacet)

COL_ODRfra = [u'ΔMFB1', u'ΔNUM1', u'ΔYPT11', u'WT_YPE', u'WT_YPL', u'WT_YPR']
COL_ODR = [u'WT_YPE', u'WT_YPL', u'WT_YPR', u'ΔMFB1', u'ΔNUM1', u'ΔYPT11']
HUE_ODR = [u'ΔMFB1', u'ΔNUM1', u'ΔYPT11', u'WT_YPE', u'WT_YPL', u'WT_YPR']
savefolder = op.expanduser(os.sep.join(['~', 'Dropbox', 'SusanneSweeShared',
                                        'figures_for_mutants2']))
mombud_dy_vars = ['DY_median_mom', 'DY_median_bud']

plt.close('all')
labelFacet = labelhandler('facet')
labelNormal = labelhandler()
labelRowFacet = labelhandler('rowfacet')


with open(op.join('cellVolume_complete.pkl'), 'rb') as INPT:
    dfsize = pickle.load(INPT)


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

dfmom['media_new'] = dfmom.media.map(newlabels)
dfbud['media_new'] = dfbud.media.map(newlabels)
# ============================================================================
# merge dfmom and dfbud, pick medium buds
# =============================================================================
cutoff = 1.5
volr = df.loc[:, ['media_new', 'bud_diameter', 'budvolratio', 'bud_dy_var']]
unbud = ((volr.budvolratio < 0.025) &
         ((volr.media_new == 'WT_YPE') | (volr.media_new == 'WT_YPL') |
          (volr.media_new == 'WT_YPR')))  # lo budvolratio WT
volr2 = volr.copy()
volr2.loc[unbud, 'budvolratio'] = 0.0

novar = (volr2.bud_dy_var.isnull()) & (volr2.bud_diameter < cutoff)
ypt = (volr2.media_new == u'ΔYPT11') & (volr2.bud_diameter < cutoff)
volr2.loc[novar, 'budvolratio'] = 0.0  # var == NAN
volr2.loc[ypt, 'budvolratio'] = 0.0  # ypt11 small
volr2 = volr2[volr2.budvolratio < 1.1]  # restrict max budvolratio to q90
volr3 = volr2[volr2.bud_diameter > cutoff]  # df based on the cutoff
num1 = (((volr3.bud_diameter < 2) | (volr3.bud_diameter > 4.3)) &
        (volr3.media_new == u'ΔNUM1'))
volr3 = volr3.loc[~num1]

dfmerge = dfmom.merge(dfbud, on='name', suffixes=['_mom', '_bud'])
dfmerge.rename(columns={0.5: 1.5,
                        u'1.0_bud': 2.0,
                        u'1.0_mom': 1.0},
               inplace=True)
dfmerge = dfmerge.loc[dfmerge.name.isin(volr3.index)]

meds = dfmerge.loc[:,
                   dfmerge.columns.isin([0.33, 0.67, 1.0, 1.5, 2.0,
                                         u'media_new_bud', u'binvol_bud'])]

alls = dfmerge.loc[:,
                   dfmerge.columns.isin([0.33, 0.67, 1.0, 1.5, 2.0,
                                         u'media__new_bud', u'binvol_bud'])]
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
mb_dy.rename({'variable': relabels}, inplace=True)

size = df.loc[:, ['media_new', 'bud_diameter', 'bud_dy_var']]

outkws1 = dict(default_ylims=[0.05, 0.95], plt_type='boxplot',
               labeller=labelNormal, col_order=COL_ODRfra)

outkws2 = dict(default_ylims=[0.1, 0.9],
               labeller=labelFacet, col_order=COL_ODR)

outkws3 = dict(default_ylims=[0.05, 0.95], plt_type='boxplot',
               labeller=labelNormal, hue_order=HUE_ODR)

# =============================================================================
#  test of sizes
# =============================================================================
callout = ['YPL_052315_025_RFPstack_036',
           'YPL_042515_021_RFPstack_029',
           'YPE_042715_013_RFPstack_048',
           'YPE_042715_007_RFPstack_033',
           'YPT11_071016_025_RFPstack_147',
           'YPT11_071016_005_RFPstack_127',
           'YPT11_071016_029_RFPstack_150']
volr2['mark'] = volr2.loc[volr2.index.intersection(callout)].budvolratio
with sns.color_palette('colorblind'):
    h = sns.FacetGrid(volr2, col='media_new',
                      col_wrap=3, col_order=HUE_ODR[:], hue='media_new',
                      size=3, aspect=1.5)
    for i in h.axes.flat:
        i.axvline(x=cutoff, linewidth=1, color='m')
    h.map(plt.scatter, 'bud_diameter', 'budvolratio')
    h.map(plt.scatter, 'bud_diameter', 'mark', s=40, color='r', marker='v')
    h.set(ylim=[0, 1.1])
    plt.savefig(op.join(savefolder, 'budratio_bud_diam.png'))

# =============================================================================
# mombud ΔΨ plot
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
    plv1.save_figure(op.join(savefolder, 'boxplot mombud dy.svg'))

# =============================================================================
# frac ΔΨ plot
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
    plv2.save_figure(op.join(savefolder, 'boxplot frac dy.svg'))

# =============================================================================
# mom, bud,  cell axis long form DF
# =============================================================================
momdy = pd.melt(dfmom,
                id_vars=['media_new', 'binvol'],
                var_name='mom axis position',
                value_name=u'ΔΨ scaled',
                value_vars=outputargs['mbax'][1:])

buddy = pd.melt(dfbud,
                id_vars=['media_new', 'binvol'],
                var_name='bud axis position',
                value_name=u'ΔΨ scaled',
                value_vars=outputargs2['mbax'].tolist()[1:])

buddy_meds = pd.melt(meds,
                     id_vars=['media_new_bud', 'binvol_bud'],
                     var_name='cell axis position',
                     value_name=u'ΔΨ scaled',
                     value_vars=meds.columns[:5].values.tolist())

allsizes = pd.melt(alls,
                   id_vars=['media_new_bud', 'binvol_bud'],
                   var_name='cell axis position',
                   value_name=u'ΔΨ scaled',
                   value_vars=meds.columns[:5].values.tolist())

# =============================================================================
# ranksum test for frac ΔΨ and ΔΨ along cell axis
# =============================================================================
# frac ΔΨ
t1 = (frac.groupby('media_new').get_group(u'WT_YPE')[u'ΔΨ scaled'].dropna())
print u'frac ΔΨ YPE ranksum tests'
print '{}'.format('*'*79)
for t in sorted(COL_ODR[1:]):
    t2 = (frac.groupby('media_new').get_group(t)[u'ΔΨ scaled'].dropna())
    _, p = ss.ranksums(t1, t2)
    print u'p-val vs {:6s} {:6.4f}'.format(t, p)

# ΔΨ along cell axis
for t in sorted(COL_ODR[1:]):
    print u'\nΔΨ WT_YPE against {} along cell axis'.format(t)
    print '{}'.format('*'*79)
    for i in [0.33, 0.67, 1.0, 1.5, 2.0]:
        t1 = (buddy_meds.groupby(['media_new_bud', 'cell axis position'])
              .get_group((u'WT_YPE', i))[u'ΔΨ scaled'].dropna())
        t2 = (buddy_meds.groupby(['media_new_bud', 'cell axis position'])
              .get_group((t, i))[u'ΔΨ scaled'].dropna())
        _, p = ss.ranksums(t1, t2)
        print 'p-val at pos {:4.2f} : {:6.4f}'.format(i, p)

# =============================================================================
# PLOTS FACETTED
# =============================================================================
buds_ype = buddy_meds.media_new_bud == u'WT_YPE'

with sns.plotting_context('talk', font_scale=1.25):
    with sns.color_palette('colorblind'):
        plt.rcParams['figure.figsize'] = (16, 11)

        set7 = dict(col_wrap=3, col='media_new_bud', hue='media_new_bud',
                    hue_order=HUE_ODR, size=4, aspect=1.2,
                    sharex=True, sharey=True, col_order=COL_ODR[:],
                    ylim=(0.0, 1.0))

        plv7 = plfacet(plt_type='pointplot', **outkws2)
        plv7.plt(data=buddy_meds,
                 mapargs=['cell axis position', u'ΔΨ scaled'],
                 **set7)
        for ax in plv7.facet_obj.axes.flat:
            wt = sns.pointplot('cell axis position', u'ΔΨ scaled',
                               data=buddy_meds.loc[buds_ype],
                               ax=ax, markers='x')
            [j.set_alpha(.5) for j in ax.collections]
            [j.set_alpha(.35) for j in ax.lines]

        plv7.facet_obj.set(ylabel='', xlabel='')
        plv7.save_figure(op.join(savefolder, 'medbuds.svg'))
