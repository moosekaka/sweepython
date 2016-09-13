# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 16:10:16 2016
Module for plots of analysis of mother bud function in budding yeast
"""
import sys
import os.path as op
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from wrappers import UsageError
import mombud.mungedf as munge
import mombud.functions.vtk_mbplots as mbfuncs

labelhandler, plviol, plbox, plfacet = (mbfuncs.labelhandler,
                                        mbfuncs.plviol,
                                        mbfuncs.plbox,
                                        mbfuncs.plfacet)
COL_ODR = munge.COL_ODR
COL_ODR = ['MFB1', 'NUM1', 'NORM. NUM1', 'YPT11', 'WT_COMBINED', 'YPL', 'YPR']
HUE_ODR = munge.HUE_ODR
savefolder = munge.datadir
mombud_dy_vars = ['DY_median_mom', 'DY_median_bud']

date_order = ['042515', '042715', '052315','032016', '071016', '071116']
normal_num1 = ['NUM1_032016_003_RFPstack_021',
               'NUM1_032016_007_RFPstack_027',
               'NUM1_032016_009_RFPstack_029',
               'NUM1_032016_013_RFPstack_032',
               'NUM1_032016_013_RFPstack_034',
               'NUM1_032016_015_RFPstack_035',
               'NUM1_071016_011_RFPstack_098',
               'NUM1_071016_013_RFPstack_099',
               'NUM1_071016_016_RFPstack_102',
               'NUM1_071016_018_RFPstack_105',
               'NUM1_071016_018_RFPstack_106',
               'NUM1_071016_020_RFPstack_107',
               'NUM1_071016_032_RFPstack_115',
               'NUM1_071016_032_RFPstack_116']

def ranger(lst, mult):
    """
    return array for y-ticks or x-tic+ks with multiples of `mult`
    from limits in `lst`
    """
    yt = [i - i % mult for i in lst]
    return np.arange(yt[0], yt[1], mult)


#def main():
#    """
#    Main
#    """
plt.close('all')
#try:
labelFacet = labelhandler('facet')
labelNormal = labelhandler()
labelRowFacet = labelhandler('rowfacet')

# args for postprocess_df(), which calls inputdata(), mungedata()
inp_args = {'inpdatpath': 'celldata.pkl',
            'mbax': np.linspace(0., 1., 6),  # pos. mom/bud cell
            'cellax': np.linspace(0, 1., 11),  # whole cell pos.
            'binsvolbud': np.linspace(0, 40, 5),  # vol bins for bud
            'binsvolmom': np.array([0, 30, 40, 80.]),
            }
outputargs = munge.postprocess_df(**inp_args)
#        df = outputargs.get('data')

with open('df_mombud_filtered.pkl', 'rb') as inpt:
    df = pickle.load(inpt)

df['media_combined'] = df.media
df.loc[df.media_combined.isin(['YPE', 'WT']), 'media_combined'] = 'WT_COMBINED'
df.loc[df.index.isin(normal_num1), 'media_combined'] = 'NORM. NUM1'

dfmom = outputargs.get('dfmom')
dfbud = outputargs.get('dfbud')
dfmom = dfmom[dfmom.index.isin(df.index)]
dfbud = dfbud[dfbud.index.isin(df.index)]
# =============================================================================
# Data long form
# =============================================================================
frac = pd.melt(df, id_vars=['media_combined'],
               value_vars=['frac'])

mb_dy = pd.melt(df, id_vars=['media', 'date'],
                value_vars=mombud_dy_vars)


momdy = pd.melt(dfmom,
                id_vars=['media', 'binvol'],
                var_name='mom axis position',
                value_name=r'$\Delta\Psi$ scaled',
                value_vars=outputargs['mbax'].tolist())

buddy = pd.melt(dfbud,
                id_vars=['media', 'binvol'],
                var_name='bud axis position',
                value_name=r'$\Delta\Psi$ scaled',
                value_vars=outputargs['mbax'].tolist())

outkws_wt = dict(default_ylims=[0.05, 0.95],
                 labeller=labelNormal, col_order=date_order)

outkws1 = dict(default_ylims=[0.05, 0.95],
               labeller=labelNormal, col_order=COL_ODR)

outkws2 = dict(default_ylims=[0.15, 0.9],
               labeller=labelFacet, col_order=COL_ODR)

set_wt = dict(x='date', y='value',
            hue='variable', group_key='date',
            title='DY by date', ylim=(0, 1.),
            )
plv_wt = plviol(**outkws_wt)
plv_wt.plt(data=mb_dy, **set_wt)
plv_wt.save_figure(op.join(savefolder, 'wt_violin_mombud.png'))


set1 = dict(x='media', y='value',
            hue='variable', group_key='media',
            title='mombud', ylim=(0, 1.),
            )
plv1 = plviol(**outkws1)
plv1.plt(data=mb_dy, **set1)
plv1.turn_off_legend()
plv1.save_figure(op.join(savefolder, 'violin_mombud.png'))

set2 = dict(x='media_combined', y='value', ylim='auto',
            hue='media_combined', group_key='media_combined', inner=None,
            title='fracDY')
plv2 = plbox(**outkws1)
plv2.plt(data=frac, **set2)
plv2.turn_off_legend()
plv2.save_figure(op.join(savefolder, 'fracDY.png'))

set5 = dict(col_wrap=4, col='media', hue='media',
            sharex=True, sharey=True, col_order=COL_ODR,
            ylim=(0.0, 1.),
            )
plv5 = plfacet(plt_type='pointplot', **outkws2)
plv5.plt(data=momdy,
         mapargs=['mom axis position', r'$\Delta\Psi$ scaled'],
         **set5)
plv5.save_figure(op.join(savefolder, 'momDY.png'))

set6 = dict(col_wrap=4, col='media', hue='media',
            sharex=True, sharey=True, col_order=COL_ODR,
            ylim=(0.0, 1.0),
            )
plv6 = plfacet(plt_type='pointplot', **outkws2)
plv6.plt(data=buddy,
         mapargs=['bud axis position', r'$\Delta\Psi$ scaled'],
         **set6)
plv6.save_figure(op.join(savefolder, 'budDY.png'))

set7 = dict(x='bin_budprog', y='frac',
            hue='media',
            title=u"Δψ vs bud progression\n ", ylim=(0., 3.),
            setargs=dict(xlabel="bud progression",
                         ylabel=u"Δψ bud/Δψ mom",)
            )

# facetted buds plot
bigbinsbud = pd.melt(dfbud,
                     id_vars=['media', 'binvol'],
                     var_name='bud axis position',
                     value_vars=dfbud.columns[:5].tolist())
bigbinsbud = bigbinsbud.dropna()

with sns.plotting_context('talk', font_scale=.85):
    ylims2 = (bigbinsbud['value']
              .quantile([0.025, 0.975]).tolist())
    yt = ranger(ylims2, 0.25)
    m0 = sns.FacetGrid(bigbinsbud,
                       row='media',
                       col="binvol",
                       hue='media',
                       row_order=COL_ODR,
                       )
    m0 = (m0.map(sns.pointplot, 'bud axis position', 'value')
          .set(yticks=yt, ylim=tuple(ylims2)))

