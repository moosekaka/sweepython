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
               'NUM1_071016_032_RFPstack_116']
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
# with open(op.join(os.getcwd(), 'df_mombud_filtered.pkl'), 'rb') as inpt:
#    df = pickle.load(inpt)
# =============================================================================
# Calculate the lags /variogram
# ==============================================================================
for mtype in sorted(vtkF.keys())[:]:
    for cell in vtkF[mtype].keys():
        if (mtype == 'NUM1' and cell in normal_num1) or mtype != 'NUM1':
            try:
                with open(op.join(rawdir,
                                  'fitted_data_scaled',
                                  '%s.pkl' % cell), 'rb') as inpt:
                    (lNorm, lNormP, randNDY,
                     randUDY, llineId) = pickle.load(inpt)
                    DYL = DYL.append(iterlagspd(lNorm, mtype),
                                     ignore_index=True)
                    SHL = SHL.append(iterlagspd(lNormP, mtype),
                                     ignore_index=True)
                    DNL = DNL.append(iterlagspd(randNDY, mtype),
                                     ignore_index=True)
                    DUL = DUL.append(iterlagspd(randUDY, mtype),
                                     ignore_index=True)
                    print "done {}".format(cell)
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
    FIG1a._legend.set(title='\n')
    plt.show()
    plt.savefig(op.join(savefolder, 'lags_random_wt.png'))

# vs random
with sns.plotting_context('talk', font_scale=1.):
    sns.set_style({"legend.markerscale": "1."})
    FIG1b = sns.factorplot(x='lags/k',
                           y='F(k)',
                           col='type',
                           hue='cat',
                           data=A_MUTANTS,
                           ci=99, size=5, aspect=1.25,
                           legend_out=True, col_wrap=2)
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
    sns.despine(top=True, right=True)
    FIG2.set_ylabel('F(k)')
    FIG2.legend(title='', loc=4)
    plt.show()
    plt.savefig(op.join(savefolder, 'lags_WT.png'))

# vs carbon type MUTANTS
with sns.plotting_context('talk', font_scale=1.25):
    _, ax1 = plt.subplots(1, 1)
    sns.set_style({"legend.markerscale": "1."})
    FIG3 = sns.pointplot(x='lags/k',
                         y='F(k)',
                         hue='cat',
                         data=MUTANTS,
                         ci=99,
                         ax=ax1)
    sns.despine(top=True, right=True)
    FIG3.set_ylabel('F(k)')
    FIG3.legend(title='', loc=4)
    plt.show()
    plt.savefig(op.join(savefolder, 'lags_MUTANTS.png'))

#vs carbon type WT + yPE
#wt_ype = WT_old_new[WT_old_new.type==u'ΔΨ scaled']
#wt_ype['F(k)'] = np.where(wt_ype.cat=='WT',
#                          wt_ype['F(k)']+.0125, wt_ype['F(k)'])
#with sns.plotting_context('talk', font_scale=1.25):
#    _, ax1 = plt.subplots(1, 1)
#    sns.set_style({"legend.markerscale": "1."})
#    FIG2a = sns.pointplot(x='lags/k',
#                         y='F(k)',
#                         hue='cat',
#                         data=wt_ype,
#                         ci=99,
#                         ax=ax1)
#    sns.despine(top=True, right=True)
#    FIG2a.set_ylabel('F(k)')
#    plt.show()
#    plt.savefig(op.join(savefolder, 'lags_WT_oldnew.png'))
