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
sns.set_style('whitegrid')
plt.rcParams['font.family'] = 'DejaVu Sans'
labelhandler, plviol, plbox, plfacet = (mbfuncs.labelhandler,
                                        mbfuncs.plviol,
                                        mbfuncs.plbox,
                                        mbfuncs.plfacet)
COL_ODR = ['MFB1', 'DEFECT. NUM1', 'NORM. NUM1', 'NUM1', 'YPT11',
           'WT', 'YPE', 'WT_COMBINED',  'YPL', 'YPR', ]

HUE_ODR = munge.HUE_ODR
savefolder = r"C:\Users\sweel_Rafelski\Dropbox\SusanneSweeShared\091316"
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


plt.close('all')
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

with open('df_mombud_filtered.pkl', 'rb') as inpt:
    df = pickle.load(inpt)

df['media2'] = np.nan
df.loc[df.media.isin(['YPE', 'WT']), 'media2'] = 'WT_COMBINED'
df.loc[df.index.isin(normal_num1), 'media2'] = 'NORM. NUM1'
num1_all = df.loc[df.media == 'NUM1'].index
num1_mutant = num1_all[~num1_all.isin(normal_num1)]
df.loc[num1_mutant, 'media2'] = 'DEFECT. NUM1'

dfmom = outputargs.get('dfmom')
dfmom = dfmom[dfmom.index.isin(df.index)].reset_index()
dfbud = outputargs.get('dfbud')
dfbud = dfbud[dfbud.index.isin(df.index)].reset_index()

dfmom['media2'] = np.nan
dfmom.loc[dfmom.media.isin(['YPE', 'WT']), 'media2'] = 'WT_COMBINED'
dfmom.loc[dfmom.name.isin(normal_num1), 'media2'] = 'NORM. NUM1'
num1_all = dfmom.loc[dfmom.media == 'NUM1'].name
num1_mutant = num1_all[~num1_all.isin(normal_num1)]
dfmom.loc[dfmom.name.isin(num1_mutant), 'media2'] = 'DEFECT. NUM1'

dfbud['media2'] = np.nan
dfbud.loc[dfbud.media.isin(['YPE', 'WT']), 'media2'] = 'WT_COMBINED'
dfbud.loc[dfbud.name.isin(normal_num1), 'media2'] = 'NORM. NUM1'
num1_all = dfbud.loc[dfbud.media == 'NUM1'].name
num1_mutant = num1_all[~num1_all.isin(normal_num1)]
dfbud.loc[dfbud.name.isin(num1_mutant), 'media2'] = 'DEFECT. NUM1'
# =============================================================================
# Data long form
# =============================================================================
frac = pd.melt(df, id_vars=['media'],
               value_vars=['frac'])
frac1 = pd.melt(df, id_vars=['media2'],
               value_vars=['frac']).dropna()
frac1.rename(columns = {'media2':'media'}, inplace=True)
frac = frac.append(frac1)

mb_dy = pd.melt(df, id_vars=['media'],
                value_vars=mombud_dy_vars)
mb_dy1 = pd.melt(df, id_vars=['media2'],
                value_vars=mombud_dy_vars).dropna()
mb_dy1.rename(columns = {'media2':'media'}, inplace=True)
mb_dy = mb_dy.append(mb_dy1)

outkws1 = dict(default_ylims=[0.05, 0.95],
               labeller=labelNormal, col_order=COL_ODR)

outkws2 = dict(default_ylims=[0.15, 0.9],
               labeller=labelFacet, col_order=COL_ODR)

 # FracDY plot
set1 = dict(x='media', y='value',
            hue='variable', group_key='media',
            title='mombud', ylim=(0, 1.),
            )

plv1 = plviol(**outkws1)
plv1.plt(data=mb_dy, **set1)
plv1.turn_off_legend()
plv1.save_figure(op.join(savefolder, 'violin_mombud.png'))

 # mombudDY plot
with sns.plotting_context('talk', font_scale=1.2):
    _, ax2 = plt.subplots(figsize=(20,16))
    set2 = dict(x='media', y='value', data=frac, width=.75,
                 order=COL_ODR, notch=True, bootstrap=100000)
    g = sns.boxplot(ax=ax2, **set2)
    g.set(ylim=(0, 3.5))
    plt.savefig(op.join(savefolder, 'fracDY.png'))

with sns.plotting_context('talk', font_scale=1.15):
    plt.rcParams['figure.figsize']=(20,16)
    set2v = dict(x='media', y='value', size=(20,16),
               group_key='media', inner=None,
                title='fracDY', ylim='auto')
    plv2 = plviol(**outkws1)
    plv2.plt(data=frac, **set2v)
    plv2.turn_off_legend()
    plv2.save_figure(op.join(savefolder, 'test.png'))

 # mom and bud cell axis DY plots
momdy = pd.melt(dfmom,
                id_vars=['media', 'binvol'],
                var_name='mom axis position',
                value_name=r'$\Delta\Psi$ scaled',
                value_vars=outputargs['mbax'].tolist())

momdy1 = pd.melt(dfmom,
                 id_vars=['media2', 'binvol'],
                 var_name='mom axis position',
                 value_name=r'$\Delta\Psi$ scaled',
                 value_vars=outputargs['mbax'].tolist()).dropna()
momdy1.rename(columns = {'media2':'media'}, inplace=True)
momdy = momdy.append(momdy1)

buddy = pd.melt(dfbud,
                id_vars=['media', 'binvol'],
                var_name='bud axis position',
                value_name=r'$\Delta\Psi$ scaled',
                value_vars=outputargs['mbax'].tolist())

buddy1 = pd.melt(dfbud,
                 id_vars=['media2', 'binvol'],
                 var_name='bud axis position',
                 value_name=r'$\Delta\Psi$ scaled',
                 value_vars=outputargs['mbax'].tolist()).dropna()
buddy1.rename(columns = {'media2':'media'}, inplace=True)
buddy = buddy.append(buddy1)


set5 = dict(col_wrap=5, col='media', hue='media',
            sharex=True, sharey=True, col_order=COL_ODR,
            ylim=(0.0, 1.),
            )
plv5 = plfacet(plt_type='pointplot', **outkws2)
plv5.plt(data=momdy,
         mapargs=['mom axis position', r'$\Delta\Psi$ scaled'],
         **set5)
plv5.save_figure(op.join(savefolder, 'momDY.png'))

set6 = dict(col_wrap=5, col='media', hue='media',
            sharex=True, sharey=True, col_order=COL_ODR,
            ylim=(0.0, 1.0),
            )
plv6 = plfacet(plt_type='pointplot', **outkws2)
plv6.plt(data=buddy,
         mapargs=['bud axis position', r'$\Delta\Psi$ scaled'],
         **set6)
plv6.save_figure(op.join(savefolder, 'budDY.png'))
