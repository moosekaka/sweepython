# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 20:39:03 2015
PLOT DY dist for one cell from four different type of populations and their
respective random fitted distributions
@author: sweel
"""
import matplotlib.pyplot as plt
import glob
import os
import cPickle as pickle
import seaborn as sns
from tubule_het.autoCor.fitDistr import fitDist as fitd
import numpy as np
import pandas as pd
sns.set_context("talk")
sns.set(style="darkgrid")
# pylint: disable=C0103

randNDY = {}
randUDY = {}

dirlist = {}
for root, dirs, files in os.walk(os.getcwd()):
    for f in dirs:
        if f.startswith('YP'):
            dirlist[f[:3]] = (os.path.join(root, f))

# =============================================================================
#   make Dataframe df of the lines for the cells from four populations
# =============================================================================
dic = {'1YPE_042515_001_RFPstack_000': 17,
       '0YPD_042515_003_RFPstack_001': 28,
       '2YPL_042515_001_RFPstack_002': 62,
       '3YPR_042715_003_RFPstack_004': 27}
df = pd.DataFrame()

for k in sorted(dic.keys()):
    labs = k[1:4]
    print'\nNow on %s' % labs+"\n"+"="*79
    files = glob.glob(dirlist[labs]+r'\Norm*vtk')
    with open(
        os.path.join(
            dirlist[labs], '%s_grph.pkl' % labs), 'rb') as inpt:
        G = pickle.load(inpt)[2]
        Graphs = {i.graph['cell']: i for i in G}
        output = fitd(files, Graphs)
        data = output[0]
        sampN, sampU, Norm, NormPermute = output[1:5]

    cell = k[1:]

    for line in range(data[cell].number_of_lines):
        M = data[cell].get_cell(line).number_of_points
        randNDY\
            .setdefault(cell, [])\
            .append(sampN[cell].rvs(size=M))
        randUDY\
            .setdefault(cell, [])\
            .append(sampU[cell].rvs(size=M))

    dftemp = pd.DataFrame({'Normal': pd.Series(randNDY[cell][dic[k]]),
                           'Uniform': pd.Series(randUDY[cell][dic[k]]),
                           'Actual': pd.Series(Norm[cell][dic[k]]),
                           'Shuffled': pd.Series(NormPermute[cell][dic[k]])})
    dftemp['media'] = labs
    df = df.append(dftemp)
df['pos'] = df.index

with sns.plotting_context('talk', font_scale=1.3):
    #  for random and actual distributions, four diff populations
    MELT = pd.melt(df,
                   id_vars=['media', 'pos'],
                   var_name='type',
                   value_name='DYscaled')

#    for random type distributions
    MELTr = MELT.loc[(MELT.type == 'Shuffled') |
                     (MELT.type == 'Normal') |
                     (MELT.type == 'Uniform')]
    FIG1 = sns.factorplot(x='pos',
                          y='DYscaled',
                          hue='media',
                          col='type',
                          data=MELTr[MELTr.pos <= 30],
                          scale=.35)
    FIG1.axes[0][0].set_ylabel(r'$\Delta \Psi$ scaled')
    for subp in FIG1.axes[0]:
        subp.set_xticks(np.arange(0, 60, 10))
        subp.set_xticklabels(np.arange(0, 60, 10))
    plt.show()

#   for actual distribution
    MELTa = MELT.loc[(MELT.type == 'Actual')]
    FIG2 = sns.factorplot(x='pos',
                          y='DYscaled',
                          hue='media',
                          col='type',
                          data=MELTa,
                          scale=.35)
    FIG2.axes[0][0].set_ylabel(r'$\Delta \Psi$ scaled')
    for subp in FIG2.axes[0]:
        subp.set_xticks(np.arange(0, 60, 10))
        subp.set_xticklabels(np.arange(0, 60, 10))
    plt.show()
