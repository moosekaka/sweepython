# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 06:04:40 2015

@author: sweel
"""
import matplotlib.pyplot as plt
import cPickle as pickle
import seaborn as sns
import pandas as pd
import numpy as np
sns.set_context("talk")
sns.set(style="darkgrid")
sns.set(rc={"legend.markerscale": 3})
plt.close('all')
colors = ["denim blue",
          "medium green"]

with open('autocorNorm.pkl', 'rb') as inpt:
    datanorm = pickle.load(inpt)

with open('autocorRFP.pkl', 'rb') as inpt:
    datarfp = pickle.load(inpt)

datarfp['valtype'] = 'background subt. RFP'
datanorm['valtype'] = r'$\Delta \Psi$ norm'
datacp = pd.concat([datarfp.loc[datarfp.thresh == 40],
                    datanorm.loc[datanorm.thresh == 40]])

with sns.plotting_context('talk', font_scale=1.4):
    FIGM = sns.factorplot(x='lag',
                          y='auto_cor',
                          col='type',
                          hue='valtype',
                          palette=sns.xkcd_palette(colors),
                          col_wrap=2,
                          scale=.5,
                          data=datacp)
    plt.show()
    FIGM.despine(left=True)
    FIGM.set_ylabels('Autocorr. Coeff.')
    for subp in FIGM.axes:
        subp.set_xticks(np.arange(0, 15, 2))
        subp.set_xticklabels(np.arange(0, 15, 2))
