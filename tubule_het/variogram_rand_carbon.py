# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 01:26:31 2015
Plot the rms of 'gradient' of two intensity edges/vectors shifted by a lag k
@author: sweel
"""
import os
import os.path as op
import cPickle as pickle
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import wrappers as wr
from tubule_het.autoCor.AutoPopFunc import iterlagspd
# import numpy as np
# pylint: disable=C0103
# pylint: disable=R0204
sns.set_context("talk")
sns.set(style="whitegrid")
sns.set(rc={"legend.markerscale": 3})
savefolder = op.expanduser(os.sep.join(['~', 'Dropbox', 'SusanneSweeShared',
                                        'figures_for_hetero']))
normal_num1 = ['NUM1_032016_003_RFPstack_021',
               'NUM1_032016_007_RFPstack_027',
               'NUM1_032016_009_RFPstack_029',
               'NUM1_032016_013_RFPstack_032',
               'NUM1_032016_013_RFPstack_034',
               'NUM1_032016_015_RFPstack_035',
               'NUM1_032716_022_RFPstack_048',
               'NUM1_032716_020_RFPstack_046',
               'NUM1_032716_018_RFPstack_042',
               'NUM1_032716_014_RFPstack_038',
               'NUM1_032716_009_RFPstack_035',
               'NUM1_071016_011_RFPstack_098',
               'NUM1_071016_013_RFPstack_099',
               'NUM1_071016_016_RFPstack_102',
               'NUM1_071016_018_RFPstack_105',
               'NUM1_071016_018_RFPstack_106',
               'NUM1_071016_020_RFPstack_107',
               'NUM1_071016_032_RFPstack_115',
               'NUM1_071016_032_RFPstack_116',
               'NUM1_032016_003_RFPstack_022',
               'NUM1_071016_009_RFPstack_096',
               'NUM1_032716_003_RFPstack_030',
               'NUM1_032016_007_RFPstack_026',
               'NUM1_032016_001_RFPstack_018',
               'NUM1_032016_003_RFPstack_023',
               'NUM1_032716_018_RFPstack_044',
               'NUM1_032016_001_RFPstack_017',
               'NUM1_032716_016_RFPstack_040',
               'NUM1_032716_009_RFPstack_034']

highvar = ['MFB1_032016_005_RFPstack_002',
           'MFB1_032016_013_RFPstack_009',
           'MFB1_032016_017_RFPstack_014',
           'MFB1_032716_001_RFPstack_001',
           'MFB1_032716_010_RFPstack_011',
           'MFB1_032716_010_RFPstack_012',
           'MFB1_032716_012_RFPstack_014',
           'MFB1_032716_016_RFPstack_018',
           'MFB1_032716_018_RFPstack_024',
           'MFB1_032716_022_RFPstack_026',
           'MFB1_071016_004_RFPstack_045',
           'MFB1_071016_006_RFPstack_046',
           'MFB1_071016_008_RFPstack_047',
           'MFB1_071016_008_RFPstack_048',
           'MFB1_071016_008_RFPstack_049',
           'MFB1_071016_010_RFPstack_051',
           'MFB1_071016_014_RFPstack_054',
           'MFB1_071016_014_RFPstack_055',
           'MFB1_071016_019_RFPstack_059',
           'MFB1_071016_021_RFPstack_063',
           'MFB1_071016_023_RFPstack_064',
           'MFB1_071016_023_RFPstack_065',
           'MFB1_071016_039_RFPstack_081',
           'MFB1_071016_039_RFPstack_082',
           'MFB1_071016_041_RFPstack_084']

wt_norm = ['WT_032016_004_RFPstack_039',
           'WT_032016_004_RFPstack_040',
           'WT_032016_004_RFPstack_041',
           'WT_032016_006_RFPstack_042',
           'WT_032016_006_RFPstack_043',
           'WT_032016_006_RFPstack_044',
           'WT_032016_008_RFPstack_045',
           'WT_032016_008_RFPstack_046',
           'WT_032016_010_RFPstack_047',
           'WT_032016_012_RFPstack_048',
           'WT_032016_012_RFPstack_049',
           'WT_032016_017_RFPstack_050',
           'WT_032016_017_RFPstack_052',
           'WT_032016_017_RFPstack_053',
           'WT_032016_019_RFPstack_054',
           'WT_032016_021_RFPstack_057',
           'WT_032016_021_RFPstack_058',
           'WT_032016_023_RFPstack_059',
           'WT_032716_001_RFPstack_051',
           'WT_032716_001_RFPstack_052',
           'WT_032716_001_RFPstack_053',
           'WT_032716_003_RFPstack_054',
           'WT_032716_003_RFPstack_055',
           'WT_032716_005_RFPstack_056',
           'WT_032716_005_RFPstack_057',
           'WT_032716_007_RFPstack_059',
           'WT_032716_007_RFPstack_060',
           'WT_032716_009_RFPstack_061',
           'WT_032716_009_RFPstack_062',
           'WT_032716_009_RFPstack_063',
           'WT_071016_004_RFPstack_001',
           'WT_071016_008_RFPstack_003',
           'WT_071016_008_RFPstack_004',
           'WT_071016_010_RFPstack_005',
           'WT_071016_010_RFPstack_006',
           'WT_071016_013_RFPstack_007',
           'WT_071016_013_RFPstack_008',
           'WT_071016_013_RFPstack_009',
           'WT_071016_013_RFPstack_010',
           'WT_071016_013_RFPstack_011',
           'WT_071016_013_RFPstack_012',
           'WT_071016_015_RFPstack_013',
           'WT_071016_017_RFPstack_014',
           'WT_071016_017_RFPstack_015',
           'WT_071016_019_RFPstack_017',
           'WT_071016_019_RFPstack_018',
           'WT_071016_021_RFPstack_019',
           'WT_071016_021_RFPstack_020',
           'WT_071016_021_RFPstack_021',
           'WT_071016_021_RFPstack_022',
           'WT_071016_023_RFPstack_023',
           'WT_071016_023_RFPstack_024',
           'WT_071016_023_RFPstack_025',
           'WT_071016_023_RFPstack_027',
           'WT_071016_025_RFPstack_028',
           'WT_071016_025_RFPstack_029',
           'WT_071016_028_RFPstack_030',
           'WT_071016_028_RFPstack_032',
           'WT_071016_030_RFPstack_035']
# =============================================================================
#    Data initialization
# =============================================================================
plt.close('all')
rawdir = op.join(os.getcwd(), 'old_w_new')
vtkF = wr.ddwalk(op.join(rawdir, 'normalizedVTK'),
                 '*skeleton.vtk', start=5, stop=-13)

DYL = pd.DataFrame()
SHL = pd.DataFrame()
DNL = pd.DataFrame()
DUL = pd.DataFrame()
# =============================================================================
# Calculate the lags /variogram
# ==============================================================================
#L2 = []
for mtype in sorted(vtkF.keys())[:]:
    for cell in vtkF[mtype].keys():
        if ((mtype == 'NUM1' and cell in normal_num1) or
            (mtype != 'NUM1' and mtype != 'MFB1' and mtype != 'WT') or
            (mtype == 'MFB1' and cell not in highvar) or
                (mtype == 'WT' and cell in wt_norm)):
            try:
                with open(op.join(rawdir,
                                  'fitted_data_scaled',
                                  '%s.pkl' % cell), 'rb') as inpt:
                    (lNorm, lNormP, randNDY,
                     randUDY, llineId) = pickle.load(inpt)
#                    temp =iterlagspd(lNorm, mtype)
#                    temp['cellname'] = cell
#                    DYL = DYL.append(temp, ignore_index=True)
                    DYL = DYL.append(iterlagspd(lNorm, mtype),
                                     ignore_index=True)
                    SHL = SHL.append(iterlagspd(lNormP, mtype),
                                     ignore_index=True)
                    DNL = DNL.append(iterlagspd(randNDY, mtype),
                                     ignore_index=True)
                    DUL = DUL.append(iterlagspd(randUDY, mtype),
                                     ignore_index=True)
                    print "done {}".format(cell)
#                    L2.append(cell)
            except IOError:
                pass

DYL['type'] = u'ΔΨ scaled'
SHL['type'] = 'Shuffled'
DNL['type'] = 'Normal Dist.'
DUL['type'] = 'Uniform Dist.'

newlabels = {'YPD': 'WT_YPD',
             'YPR': u'WT_YPR', 'YPL': u'WT_YPL', 'YPE': u'WT_YPE',
             'WT': u'WT_YPE', 'NUM1': u'ΔNUM1', 'MFB1': u'ΔMFB1',
             'YPT11': u'ΔYPT11'}

BIG = pd.concat([DYL, SHL, DUL, DNL], ignore_index=True)
A = pd.melt(BIG,
            id_vars=['cat', 'type'],
            var_name='lags/k',
            value_name='F(k)')
WT_old_new = A[A.cat.isin([u'YPR', u'WT', 'YPE', 'YPD', 'YPL'])]
A = A.replace({'cat': newlabels})
A_WT = A[A.cat.isin([u'WT_YPR', u'WT_YPL', 'WT_YPE', 'WT_YPD'])]
A_MUTANTS = A[A.cat.isin([u'ΔNUM1', u'ΔMFB1', 'WT_YPE', u'ΔYPT11'])]
MASK = A['type'] == u'ΔΨ scaled'
B = A[MASK]
C = B.replace({'cat': newlabels})
WT = C[C.cat.isin([u'WT_YPR', u'WT_YPL', 'WT_YPE', 'WT_YPD'])]
MUTANTS = C[C.cat.isin([u'ΔNUM1', u'ΔMFB1', 'WT_YPE', u'ΔYPT11'])]
# =============================================================================
# Plots
# =============================================================================
sns.set(style='white')
# vs random
with sns.plotting_context('talk', font_scale=1.):
    sns.set_style({"legend.markerscale": "1."})
    FIG1a = sns.factorplot(x='lags/k',
                           y='F(k)',
                           col='type',
                           hue='cat',
                           data=A_WT,
                           ci=99, size=5, aspect=1.25,
                           legend_out=True, col_wrap=2)
    for ax in FIG1a.axes.flat:
        [j.set_alpha(.7) for j in ax.axes.collections]
        [j.set_alpha(.7) for j in ax.axes.lines]
    FIG1a._legend.set(title='\n')
    plt.show()
    plt.savefig(op.join(savefolder, 'lags_random_wt.png'))

# vs random MUTANTS
with sns.plotting_context('talk', font_scale=1.):
    with sns.color_palette('colorblind'):
        sns.set_style({"legend.markerscale": "1."})
        FIG1b = sns.factorplot(x='lags/k',
                               y='F(k)',
                               col='type',
                               hue='cat',
                               data=A_MUTANTS,
                               ci=99, size=5, aspect=1.25,
                               legend_out=True, col_wrap=2)
        for ax in FIG1b.axes.flat:
            [j.set_alpha(.7) for j in ax.axes.collections]
            [j.set_alpha(.7) for j in ax.axes.lines]
        FIG1b._legend.set(title='\n')
        plt.show()
        plt.savefig(op.join(savefolder, 'lags_random_mut.png'))

# vs carbon type WT
with sns.plotting_context('talk', font_scale=1.25):
    _, ax1 = plt.subplots(1, 1)
    sns.set_style({"legend.markerscale": "1."})
    FIG2 = sns.pointplot(x='lags/k',
                         y='F(k)',
                         hue='cat',
                         data=WT,
                         ci=99,
                         ax=ax1)
    [j.set_alpha(.7) for j in FIG2.axes.collections]
    [j.set_alpha(.7) for j in FIG2.axes.lines]
    sns.despine(top=True, right=True)
    FIG2.set_ylabel('F(k)')
    FIG2.legend(title='', loc=4)
    plt.show()
    plt.savefig(op.join(savefolder, 'lags_WT.png'))

# vs carbon type MUTANTS
with sns.plotting_context('talk', font_scale=1.25):
    with sns.color_palette('colorblind'):
        _, ax1 = plt.subplots(1, 1)
        sns.set_style({"legend.markerscale": "1."})
        FIG3 = sns.pointplot(x='lags/k',
                             y='F(k)',
                             hue='cat',
                             data=MUTANTS,
                             ci=99,
                             ax=ax1)
        [j.set_alpha(.7) for j in FIG3.axes.collections]
        [j.set_alpha(.7) for j in FIG3.axes.lines]
        sns.despine(top=True, right=True)
        FIG3.set_ylabel('F(k)')
        FIG3.legend(title='', loc=4)
        plt.show()
        plt.savefig(op.join(savefolder, 'lags_MUTANTS.png'))
