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
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wrappers as wr
from tubule_het.autoCor.AutoPopFunc import iterlagspd
# import numpy as np
# pylint: disable=C0103
# pylint: disable=R0204
HUE_ODR_WT = ['WT_YPD', u'WT_YPE', u'WT_YPL', u'WT_YPR']
HUE_ODR = [u'ΔMFB1', u'WT_YPE', u'ΔNUM1', u'ΔYPT11']
sns.set_context("talk")
sns.set(style="whitegrid")
sns.set(rc={"legend.markerscale": 3})
savefolder = op.expanduser(os.sep.join(['~', 'Dropbox', 'SusanneSweeShared',
                                        'figures_for_hetero']))
wtypes = [u'WT_YPR', u'WT_YPL', 'WT_YPE', 'WT_YPD']
mtypes = [u'ΔNUM1', u'ΔMFB1', 'WT_YPE', u'ΔYPT11']
newlabels = {'YPD': 'WT_YPD',
             'YPR': u'WT_YPR', 'YPL': u'WT_YPL', 'YPE': u'WT_YPE',
             'WT': u'WT_YPE', 'NUM1': u'ΔNUM1', 'MFB1': u'ΔMFB1',
             'YPT11': u'ΔYPT11'}
# =============================================================================
#    Data initialization (DataFrame of paths to cell files)
# =============================================================================
plt.close('all')
rawdir = op.join(os.getcwd(), 'old_w_new')
vtkF = wr.swalk(op.join(rawdir, 'normalizedVTK'),
                '*skeleton.vtk', start=5, stop=-13)

df = pd.DataFrame(pd.Series(vtkF))
df['media'] = df.index.map(lambda x: x.partition('_')[0])
df.rename(columns={0: 'path'}, inplace=True)
#==============================================================================
# filtering out dataset of interest (mutants cells)
#==============================================================================
with open('list_hetero.pkl', 'rb') as inpt:
    normal_num1, highvar, wt_norm = pickle.load(inpt)

#  get rid of mfb1 highvar
mfb1 = df.index.symmetric_difference(highvar)
df = df.loc[mfb1]
#  filter for wt
wt = df.loc[df.media == 'WT'].index.difference(wt_norm)
df = df.loc[df.index.symmetric_difference(wt)]
#  filter for num1
num1 = df.loc[df.media == 'NUM1'].index.difference(normal_num1)
df = df.loc[df.index.symmetric_difference(num1)]
# =============================================================================
# Calculate the lags /variogram
# ==============================================================================
DYL = pd.DataFrame()
SHL = pd.DataFrame()
DNL = pd.DataFrame()
DUL = pd.DataFrame()
#L2 = []
try:
    with open('tubule_hetero.pkl', 'rb') as inpt:
        BIG = pickle.load(inpt)
        print "pickle found, data loaded!"

except IOError:  # if the pkl file does not exist regen data
    print "pickle not found, regen data\n"
    for mtype, grp in df.groupby('media'):
        print 'Done {}\n{}'.format(mtype, '*' * 79)
        for path in grp.index:
            try:
                with open(op.join(rawdir,
                                  'fitted_data_scaled',
                                  '%s.pkl' % path), 'rb') as inpt:
                    (lNorm, lNormP, randNDY,
                     randUDY, llineId) = pickle.load(inpt)
    #                    temp =iterlagspd(lNorm, mtype)
    #                    temp['cellname'] = path
    #                    DYL = DYL.append(temp, ignore_index=True)
                    DYL = DYL.append(iterlagspd(lNorm, mtype),
                                     ignore_index=True)
                    SHL = SHL.append(iterlagspd(lNormP, mtype),
                                     ignore_index=True)
                    DNL = DNL.append(iterlagspd(randNDY, mtype),
                                     ignore_index=True)
                    DUL = DUL.append(iterlagspd(randUDY, mtype),
                                     ignore_index=True)
#                    print "done {}".format(path)
    #                    L2.append(path)
            except IOError:
                pass

    DYL['type'] = u'ΔΨ scaled'
    SHL['type'] = 'Shuffled'
    DNL['type'] = 'Normal Dist.'
    DUL['type'] = 'Uniform Dist.'
    BIG = pd.concat([DYL, SHL, DUL, DNL], ignore_index=True)

alldist = pd.melt(BIG,
                  id_vars=['cat', 'type'],
                  var_name='lags/k',
                  value_name='F(k)')
alldist.replace({'cat': newlabels}, inplace=True)
alldist_WT = alldist[alldist.cat.isin(wtypes)]
alldist_MUTANTS = alldist[alldist.cat.isin(mtypes)]
real = alldist['type'] == u'ΔΨ scaled'
realdist = alldist[real]
WT = realdist[realdist.cat.isin(wtypes)]
MUTANTS = realdist[realdist.cat.isin(mtypes)]

# =============================================================================
# Plots
# =============================================================================
sns.set(style='white')
# vs random
typefilter = alldist_WT.type.isin(alldist_WT.type.unique()[1:])
with sns.plotting_context('talk', font_scale=1.):
    sns.set_style({"legend.markerscale": "1."})
    FIG1a = sns.factorplot(x='lags/k',
                           y='F(k)',
                           row='type',
                           hue='cat',
                           hue_order=HUE_ODR_WT,
                           data=alldist_WT[typefilter],
                           ci=99, size=5, aspect=1.25,
                           legend=False)
    for ax in FIG1a.axes.flat:
        [j.set_alpha(.7) for j in ax.axes.collections]
        [j.set_alpha(.7) for j in ax.axes.lines]
    FIG1a.set(yticks=np.arange(0, 0.22, 0.05))
    plt.legend(title='', loc=8)
    plt.show()
    plt.savefig(op.join(savefolder, 'lags_random_wt.svg'))

# vs random MUTANTS
typefilter = alldist_MUTANTS.type.isin(alldist_MUTANTS.type.unique()[1:])
with sns.plotting_context('talk', font_scale=1.):
    with sns.color_palette('colorblind'):
        sns.set_style({"legend.markerscale": "1."})
        FIG1b = sns.factorplot(x='lags/k',
                               y='F(k)',
                               row='type',
                               hue='cat',
                               data=alldist_MUTANTS[typefilter],
                               hue_order=HUE_ODR,
                               ci=99, size=5, aspect=1.25,
                               legend=False)
        for ax in FIG1b.axes.flat:
            [j.set_alpha(.7) for j in ax.axes.collections]
            [j.set_alpha(.7) for j in ax.axes.lines]
        FIG1b.set(yticks=np.arange(0, 0.22, 0.05))
        plt.legend(title='', loc=8)
        plt.show()
        plt.savefig(op.join(savefolder, 'lags_random_mut.svg'))

# vs carbon type WT
with sns.plotting_context('talk', font_scale=1.25):
    _, ax1 = plt.subplots(1, 1)
    sns.set_style({"legend.markerscale": "1."})
    FIG2 = sns.pointplot(x='lags/k',
                         y='F(k)',
                         hue='cat',
                         hue_order=HUE_ODR_WT,
                         data=WT,
                         ci=99,
                         ax=ax1)
    [j.set_alpha(.7) for j in FIG2.axes.collections]
    [j.set_alpha(.7) for j in FIG2.axes.lines]
    sns.despine(top=True, right=True)
    FIG2.set_ylabel('F(k)')
    FIG2.legend(title='', loc=8)
    FIG2.set(yticks=np.arange(0, 0.22, 0.05))
    plt.show()
    plt.savefig(op.join(savefolder, 'lags_WT.svg'))

# vs carbon type MUTANTS
with sns.plotting_context('talk', font_scale=1.25):
    with sns.color_palette('colorblind'):
        _, ax1 = plt.subplots(1, 1)
        sns.set_style({"legend.markerscale": "1."})
        FIG3 = sns.pointplot(x='lags/k',
                             y='F(k)',
                             hue='cat',
                             data=MUTANTS,
                             hue_order=HUE_ODR,
                             ci=99,
                             ax=ax1)
        [j.set_alpha(.7) for j in FIG3.axes.collections]
        [j.set_alpha(.7) for j in FIG3.axes.lines]
        sns.despine(top=True, right=True)
        FIG3.set_ylabel('F(k)')
        FIG3.legend(title='', loc=8)
        FIG3.set(yticks=np.arange(0, 0.22, 0.05))
        plt.show()
        plt.savefig(op.join(savefolder, 'lags_MUTANTS.svg'))
