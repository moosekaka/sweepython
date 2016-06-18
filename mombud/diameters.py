# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 00:50:45 2016

@author: sweel_Rafelski
"""
from collections import defaultdict
import numpy as np
import pandas as pd
import seaborn as sns
from mombud.functions import vtk_mbfuncs as vf
# pylint: disable=C0103
COL_ODR = ['MFB1', 'NUM1', 'YPT11', 'WT', 'YPE', 'YPL', 'YPR']

params = {'inpdatpath': 'celldata.pkl',
          'outdatapath': 'bootneck.pkl'}
vtkdf = vf.wrapper(**params)

diameters = defaultdict(dict)
for k in sorted(vtkdf.keys()):
    cell = vtkdf[k]['df']
    celldims = cell.groupby(['type']).x.agg([np.max, np.min])
    diameters[k] = (celldims.amax-celldims.amin).to_dict()

alldims = pd.DataFrame.from_dict(diameters, orient='index')
alldims = alldims.reset_index()
alldims['type'] = alldims['index'].apply(lambda x: x.split('_')[0])
melt = pd.melt(alldims,
               id_vars='type',
               value_vars=['bud', 'mom'])

with sns.plotting_context('talk'):

    g = sns.FacetGrid(melt,
                      col='type',
                      col_wrap=4,
                      hue="variable",
                      col_order=COL_ODR,
                      size=3,
                      aspect=1.5)
    g = (g.map(sns.stripplot,
               "value", jitter=0.1)).set(xlim=(0.),
                                         xlabel='diameter/microns')
## _____________________________________________________________________________
#if __name__ == '__main__':
#    main()
