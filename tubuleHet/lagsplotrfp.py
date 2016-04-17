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
from collections import defaultdict
from tubuleHet.autoCor.lagsfun import iterlagspd
sns.set_context("talk")
plt.close('all')


def sclminmax(data):
    """return a scaled min max collection of vtk cellarray data
    """
    vtkscaled = defaultdict(dict)
    for key in data.keys():
        flat = [el for lis in data[key] for el in lis]
        fmax = max(flat)
        fmin = min(flat)
        vtkscaled[key] = []
        for line in data[key]:
            vtkscaled[key].append([(el - fmin) /
                                   (fmax - fmin) for el in line])
    return vtkscaled


# =============================================================================
#    Data input
# =============================================================================
# pylint: disable=C0103
dirlist = []
# pylint: disable=C0103
for root, dirs, files in os.walk(os.getcwd()):
    for f in dirs:
        if f.startswith('YP'):
            dirlist.append(
                os.path.join(root, f))

DYL = pd.DataFrame()
SHL = pd.DataFrame()
DNL = pd.DataFrame()
DUL = pd.DataFrame()

for media in dirlist[:]:
    labs = media[-3:]
    print'\nNow on %s' % labs + "\n" + "="*79
    #   make sure the pkl file below exists, run MakeInputForLags.py otherwise
    with open('%s_lagsRFP.pkl' % labs, 'rb') as inpt:
        (randNDY, randUDY, Norm, NormPermute, data) = pickle.load(inpt)

    randNDY = sclminmax(randNDY)
    randUDY = sclminmax(randUDY)
    Norm = sclminmax(Norm)
    NormPermute = sclminmax(NormPermute)
    dfDY = pd.DataFrame({cells: pd.Series([edge for edge
                                           in Norm[cells]]) for cells
                         in Norm.keys()})
    dfShf = pd.DataFrame({cells: pd.Series([edge for edge
                                            in NormPermute[cells]]) for cells
                          in Norm.keys()})
    dfU = pd.DataFrame({cells: pd.Series([edge for edge
                                          in randUDY[cells]]) for cells
                        in Norm.keys()})
    dfN = pd.DataFrame({cells: pd.Series([edge for edge
                                          in randNDY[cells]]) for cells
                        in Norm.keys()})
    cols1 = dfDY.columns
    cols2 = dfShf.columns
    cols3 = dfU.columns
    cols4 = dfN.columns

    pdDY = iterlagspd(dfDY, cols1, labs)
    print'DY complete'
    pdSh = iterlagspd(dfShf, cols2, labs)
    print'Shuffle complete'
    pdU = iterlagspd(dfU, cols3, labs)
    print'Uniform complete'
    pdN = iterlagspd(dfN, cols4, labs)
    print'Normal complete'
    DYL = DYL.append(pdDY, ignore_index=True)
    SHL = SHL.append(pdSh, ignore_index=True)
    DUL = DUL.append(pdU, ignore_index=True)
    DNL = DNL.append(pdN, ignore_index=True)

DYL['type'] = r'$\Delta \Psi$ expt.'
SHL['type'] = 'Shuffled'
DNL['type'] = 'Normal Dist.'
DUL['type'] = 'Uniform Dist.'

BIG = pd.concat([DYL, SHL, DUL, DNL], ignore_index=True)
A = pd.melt(BIG,
            id_vars=['cat', 'type'],
            var_name='lags/k',
            value_name='F(k)')

MASK = A['type'] == r'$\Delta \Psi$ expt.'
B = A[MASK]

with sns.plotting_context('talk', font_scale=1.25):
    FIG1 = sns.factorplot(x='lags/k',
                          y='F(k)',
                          col='type',
                          hue='cat',
                          data=A,
                          col_wrap=2,
                          ci=99,
                          scale=.65)
    plt.show()

sns.set(style='white')
plt.figure()
with sns.plotting_context('talk', font_scale=1.25):
    FIG2 = sns.pointplot(x='lags/k',
                         y='F(k)',
                         hue='cat',
                         data=B,
                         ci=99)
    sns.despine(top=True, right=True)
    FIG2.set_ylabel('F(k)')
    plt.show()
