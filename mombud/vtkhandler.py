# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 16:10:16 2016
Module for plots of analysis of mother bud function in budding yeast
"""
import sys
import os
import os.path as op
import inspect
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from wrappers import UsageError
import mombud.vtk_mom_bud_analyse_refactored as mb
from mombud.vtk_mom_bud_analyse_refactored import postprocess_df, _dySet
from mombud.functions.vtk_mbplots import (labelhandler,
                                          plviol, plbox, plfacet)

COL_ODR = ['MFB1', 'NUM1', 'YPT11', 'WT', 'YPE', 'YPL', 'YPR']
HUE_ODR = ['DY_abs_mean_mom', 'DY_abs_mean_bud', 'whole_cell_abs']
savefolder = mb.datadir
gfp_plot_vars =  ['DY_abs_mean_mom',
                  'DY_abs_mean_bud',
                  'DY_abs_cell_mean']
mombud_dy_vars = ['DY_median_mom', 'DY_median_bud']

def main():
    plt.close('all')
    try:
        labelNormal = labelhandler()
        labelBudVol = labelhandler('rowfacet')
        labelFacet = labelhandler('facet')
        labelRsqr = labelhandler('rsqr')

        inp_args = {'regen': False,
                    'save': False,  # toggle to save plots
                    'inpdatpath': 'celldata.pkl',
                    'mbax': np.linspace(0., 1., 6),  # pos. along mom/bud cell
                    'cellax': np.linspace(0, 1., 11),  # position along whole cell
                    'binsvolbud': np.linspace(0, 40, 5),  # vol binning for bud
                    'binsvolmom': np.array([0, 30, 40, 80.]),

                    }

        outputargs = postprocess_df(**inp_args)  # calls inputdata(), mungedata()
        df = outputargs.get('data')
#==============================================================================
# Data long form
#==============================================================================
        frac = pd.melt(df, id_vars=['media'],
                       value_vars=['frac'])
        mb_dy = pd.melt(df, id_vars=['media', 'date'],
                        value_vars=mombud_dy_vars)
        diameters = pd.melt(df, id_vars='media',
                            value_vars=['bud_diameter', 'mom_diameter'])

#==============================================================================
# mutant only subset
#==============================================================================
        mutants = (mb_dy[mb_dy['media']
                   .isin(['NUM1', 'MFB1', 'YPT11', 'WT'])]
                   .reset_index(drop=True))
        group = mutants.groupby('date')
        outkws2 = dict(default_ylims=[0.05, 0.95],
                       labeller=labelNormal)
        set1 = dict(data=mutants, x='media', y='value',
                    hue='variable', no_legend=False, #save=True,
                    group_key='media')

        plv = plviol(**outkws2)
        plv.plt(**set1)

        set2 = dict(data=frac, x='media', y='value',
                    hue='media',# save=True,
                    group_key='media', inner=None, no_legend=True)
        plv2 = plbox(**outkws2)
        plv2.plt(**set2)

        set3 = dict(data=diameters, col_wrap=4, col='media',
                    size=3, aspect=1.5, hue='variable',  # save=True,
                    mapargs=['value',],
                    setargs=dict(xlim=(0.), xlabel='diameter'))
        plv3 = plfacet()
        plv3.plt(**set3)

        print "Finished plotting!"
        return 0

    except UsageError as e:
        print e
        return 1


if __name__=='__main__':
    sys.exit(main())
