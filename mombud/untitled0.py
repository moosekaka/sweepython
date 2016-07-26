# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 21:43:21 2016

@author: sweel_Rafelski
"""
import os
import os.path as op
import cPickle as pickle
import seaborn as sns
import numpy as np
import pandas as pd
import mombud.mungedf as mdf
import mombud.functions.vtk_mbplots as mbfuncs
labelhandler, plviol, plbox, plfacet = (mbfuncs.labelhandler,
                                        mbfuncs.plviol,
                                        mbfuncs.plbox,
                                        mbfuncs.plfacet)
labelFacet = labelhandler('facet')
labelNormal = labelhandler()
# pylint: disable=C0103
COL_ODR = ['MFB1', 'NUM1', 'YPT11', 'WT', 'YPE', 'YPL', 'YPR']
HUE_ODR = ['DY_abs_mean_mom', 'DY_abs_mean_bud', 'whole_cell_abs']
datadir = op.join(os.getcwd(), 'mutants', 'transformedData', 'filtered')
datadir_old = op.join(os.getcwd(), 'data', 'transformedData')
mutants=['MFB1', 'NUM1', 'WT', 'YPT11']
mombud_dy_vars = ['DY_median_mom', 'DY_median_bud']

kwargs = {}
def_args = {'regen':True,
            'save': False,
            'inpdatpath': 'celldata.pkl',
            'mbax': np.linspace(0., 1., 6),
            'cellax': np.linspace(0, 1., 11),
            'binsvolbud': np.linspace(0, 40, 5),
            'binsvolmom': np.array([0, 30, 40, 80.])}

def_args.update(kwargs)
outputargs = mdf.postprocess_df(**def_args)
df = outputargs['data']
df1 = df[df.media.isin(mutants)]

subset = {}
for med in ['MFB1', 'NUM1', 'WT', 'YPT11']:
    subset[med] = df1[(df1['media'] == med) & (df1['date'] == '071016')]

rand_subset = {}
for i in subset.keys():
    rand_subset[i] = subset[i].sample(frac=0.5)

plt.close('all')
changedates = [vals for key in D.keys() for vals in D[key].index.values]
rowind = df1.index[df1.index.isin(changedates)]
df1.loc[rowind, 'date'] = '071116'
#df1 = df1[~(df1.date=='032716')]
dat = pd.melt(df1, ['media', 'date'], mombud_dy_vars)


seta = dict(x='media', y='value', hue='variable')
plva = sns.FacetGrid(data=dat, row='date')
plva.map(sns.violinplot, 'media', 'value', 'variable',
         order=COL_ODR[:4]).set(ylim=(0, 1.))
