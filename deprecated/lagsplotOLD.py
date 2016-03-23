# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 01:26:31 2015
Plot the rms of 'gradient' of two intensity edges/vectors shifted by a lag k
@author: sweel
"""
import matplotlib.pyplot as plt
import os
import cPickle as pickle
import seaborn as sns
import pandas as pd
from lags.lagsFunctions import iterLagsPD
sns.set_context("talk")
plt.close('all')

# =============================================================================
#    Data input
# =============================================================================
dirList = []
for root, dirs, files in os.walk(os.getcwd()):
    for f in dirs:
        dirList.append(
            os.path.join(root, f))

pdDYL = pd.DataFrame()
pdShL = pd.DataFrame()
pdNL = pd.DataFrame()
pdUL = pd.DataFrame()

for media in dirList[:]:
    labs = media[-3:]
    print('\nNow on %s' % labs + "\n" + "="*79)
    #   make sure the pkl file below exists, run MakeInputForLags.py otherwise
    with open('%s_lags.pkl' % labs, 'rb') as inpt:
        (randNDY, randUDY, Norm, NormPermute) = pickle.load(inpt)

    framesDY = [pd.Series([edge for edge in Norm[cells]],
                          name=cells)
                for cells in Norm.keys()]

    framesShf = [pd.Series([edge for edge in NormPermute[cells]],
                           name=cells)
                 for cells in NormPermute.keys()]

    framesU = [pd.Series([edge for edge in randUDY[cells]],
                         name=cells)
               for cells in randUDY.keys()]

    framesN = [pd.Series([edge for edge in randNDY[cells]],
                         name=cells)
               for cells in randNDY.keys()]

    dfDY = pd.concat(framesDY, axis=1)
    cols1 = dfDY.columns
    dfShf = pd.concat(framesShf, axis=1)
    cols2 = dfShf.columns
    dfU = pd.concat(framesU, axis=1)
    cols3 = dfU.columns
    dfN = pd.concat(framesN, axis=1)
    cols4 = dfN.columns

#    tempDY = iterLags(dfDY, cols1)
#    tempShf = iterLags(dfShf, cols2)
#    tempU = iterLags(dfU, cols3)
#    tempN = iterLags(dfN, cols4)
#
#    tempDY = iterLagsPD(dfDY, cols1)
#    tempShf = iterLagsPD(dfShf, cols2)
#    tempU = iterLagsPD(dfU, cols3)
#    tempN = iterLagsPD(dfN, cols4)

#    print('Stacking DY')
#    d1 = stackMeans(tempDY, cols1, labs)
#    print('Stacking Shf')
#    d2 = stackMeans(tempShf, cols2, labs)
#    print('Stacking Uni')
#    d3 = stackMeans(tempU, cols3, labs)
#    print('Stacking Norm.')
#    d4 = stackMeans(tempN, cols4, labs)
#
#    pdDY = pdDY.append(d1, ignore_index=True)
#    pdSh = pdSh.append(d2, ignore_index=True)
#    pdU = pdU.append(d3, ignore_index=True)
#    pdN = pdN.append(d4, ignore_index=True)
    pdDY = iterLagsPD(dfDY, cols1, labs)
    print('DY complete')
    pdSh = iterLagsPD(dfShf, cols2, labs)
    print('Shuffle complete')
    pdU = iterLagsPD(dfU, cols3, labs)
    print('Uniform complete')
    pdN = iterLagsPD(dfN, cols4, labs)
    print('Normal complete')
    pdDYL = pdDYL.append(pdDY, ignore_index=True)
    pdShL = pdShL.append(pdSh, ignore_index=True)
    pdUL = pdUL.append(pdU, ignore_index=True)
    pdNL = pdNL.append(pdN, ignore_index=True)

pdDYL['type'] = r'$\Delta \Psi$ expt.'
pdShL['type'] = 'Shuffled'
pdNL['type'] = 'Normal Dist.'
pdUL['type'] = 'Uniform Dist.'

big = pd.concat([pdDYL, pdShL, pdUL, pdNL], ignore_index=True)
A = pd.melt(big,
            id_vars=['cat', 'type'],
            var_name='lags/k',
            value_name='F(k)')

mask = A['type'] == '$\Delta \Psi$ expt.'
B = A[mask]

#with sns.axes_style("dark"):
with sns.plotting_context('talk', font_scale=1.25):
    g = sns.factorplot(x='lags/k',
                       y='F(k)',
                       col='cat',
                       hue='type',
                       data=A,
                       col_wrap=2,
                       ci=95,
                       scale=.65)

#    g.despine(right=True, top=True)
    plt.show()
    plt.savefig('lags2.png')

sns.set(style='white')
plt.figure()
with sns.plotting_context('talk', font_scale=1.25):
    g1 = sns.pointplot(x='lags/k',
                       y='F(k)',
                       hue='cat',
                       data=B,
                       palette='GnBu_d')
    sns.despine(top=True, right=True)
    g1.set_ylabel('F(k)')
    plt.show()
    plt.savefig('lags3.png')
