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
from tubuleHet.autoCor.AutoPopFunc import iterlagspd
# pylint: disable=C0103
# pylint: disable=R0204
sns.set_context("talk")
sns.set(style="whitegrid")
sns.set(rc={"legend.markerscale": 3})

# =============================================================================
#    Data initialization
# =============================================================================
plt.close('all')
datadir = op.join(os.getcwd(), 'data')
rawdir = op.join(os.getcwd(), 'output')
vtkF = wr.ddwalk(op.join(rawdir, 'normalizedVTK'),
                 '*skeleton.vtk', start=5, stop=-13)

DYL = pd.DataFrame()
SHL = pd.DataFrame()
DNL = pd.DataFrame()
DUL = pd.DataFrame()

# =============================================================================
# Calculate the lags /variogram
# ==============================================================================
for mtype in sorted(vtkF.keys())[:]:
    for cell in vtkF[mtype].keys():
        with open(op.join(rawdir,
                          'fitted_data_scaled',
                          '%s.pkl' % cell), 'rb') as inpt:
            (lNorm, lNormP, randNDY, randUDY, llineId) = pickle.load(inpt)

            DYL = DYL.append(iterlagspd(lNorm, mtype), ignore_index=True)
            SHL = SHL.append(iterlagspd(lNormP, mtype), ignore_index=True)
            DNL = DNL.append(iterlagspd(randNDY, mtype), ignore_index=True)
            DUL = DUL.append(iterlagspd(randUDY, mtype), ignore_index=True)
            print "done {}".format(cell)

DYL['type'] = r'$\Delta \Psi$ actual'
SHL['type'] = 'Shuffled'
DNL['type'] = 'Normal Dist.'
DUL['type'] = 'Uniform Dist.'

BIG = pd.concat([DYL, SHL, DUL, DNL], ignore_index=True)
A = pd.melt(BIG,
            id_vars=['cat', 'type'],
            var_name='lags/k',
            value_name='F(k)')

MASK = A['type'] == r'$\Delta \Psi$ actual'
B = A[MASK]

# =============================================================================
# Plots
# =============================================================================
sns.set(style='white')
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
    plt.savefig('lags2.png')

with sns.plotting_context('talk', font_scale=1.25):
    _, ax1 = plt.subplots(1, 1)
    FIG2 = sns.pointplot(x='lags/k',
                         y='F(k)',
                         hue='cat',
                         data=B,
                         ci=99,
                         ax=ax1)
    sns.despine(top=True, right=True)
    FIG2.set_ylabel('F(k)')
    plt.show()
    plt.savefig('lags3.png')

g = sns.FacetGrid(B, col="cat", col_wrap=2)
g = g.map(sns.violinplot, 'lags/k', 'F(k)')
