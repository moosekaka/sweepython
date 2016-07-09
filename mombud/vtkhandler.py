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
from mombud.vtk_mom_bud_analyse_refactored import postprocess_df
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
    """
    Main

    kwargs
    ------
    regen, save : Bool
        toggle for vf.gen_data()
    inpdatpath: Str
        path for vf.gen_data individual celldf pickled data
    cellax, mbax : np array
        range for mombud and cell-axis axes
    binsvolbud, binsvolmom: np array
        cell size bins
    dy_dict : Str
        type of DY for mombud plots

    """
    plt.close('all')
    try:
        labelFacet = labelhandler('facet')
        labelNormal = labelhandler()
        labelRowFacet = labelhandler('rowfacet')
        labelRsqr = labelhandler('rsqr')

        inp_args = {'inpdatpath': 'celldata.pkl',
                    'mbax': np.linspace(0., 1., 6),  # pos. mom/bud cell
                    'cellax': np.linspace(0, 1., 11),  # whole cell pos.
                    'binsvolbud': np.linspace(0, 40, 5),  # vol bins for bud
                    'binsvolmom': np.array([0, 30, 40, 80.]),
                    }
        # calls inputdata(), mungedata()
        outputargs = postprocess_df(**inp_args)
        df = outputargs.get('data')
        dfmom = outputargs.get('dfmom')
        dfbud= outputargs.get('dfbud')
#==============================================================================
# Data long form
#==============================================================================
        frac = pd.melt(df, id_vars=['media'],
                       value_vars=['frac'])
        mb_dy = pd.melt(df, id_vars=['media', 'date'],
                        value_vars=mombud_dy_vars)
        diameters = pd.melt(df, id_vars='media',
                            value_vars=['bud_diameter', 'mom_diameter'])
        size = pd.melt(df.loc[:, ['media','budvol', 'momvol']],'media')
        momdy = pd.melt(dfmom,
                        id_vars=['media', 'binvol'],
                        var_name='mom axis position',
                        value_name=r'$\Delta\Psi$ scaled gradient',
                        value_vars=outputargs['mbax'].tolist())
        momdy = momdy.dropna()
        buddy= pd.melt(dfbud,
                       id_vars=['media', 'binvol'],
                       var_name='bud axis position',
                       value_name=r'$\Delta\Psi$ scaled gradient',
                       value_vars=outputargs['mbax'].tolist())
        buddy = buddy.dropna()
    # mutant only subset
        mutants = (mb_dy[mb_dy['media']
                   .isin(['NUM1', 'MFB1', 'YPT11', 'WT'])]
                   .reset_index(drop=True))
        group = mutants.groupby('date')
#==============================================================================
        outkws1 = dict(default_ylims=[0.025, 0.975],
                       labeller=labelNormal, col_order=COL_ODR)

        outkws2 = dict(default_ylims=[0.01, 0.99],
                       labeller=labelFacet, col_order=COL_ODR)

        set1 = dict(data=mutants, x='media', y='value',
                    hue='variable', group_key='media',
                    title='mombud', ylims=(0, 1.)
                    )
        plv = plviol(**outkws1)
        plv.plt(**set1)
        plv.turn_off_legend()
        plv.save_figure(op.join(savefolder, 'mutantsDY_mombud.png'))

        set2 = dict(data=frac, x='media', y='value',
                    hue='media', group_key='media', inner=None
                    )
        plv2 = plbox(**outkws1)
        plv2.plt(**set2)
        plv2.turn_off_legend()
        plv2.save_figure(op.join(savefolder, 'fracDY.png'))

        set3 = dict(data=diameters, col_wrap=4, col='media',
                    size=3, aspect=1.5, hue='variable',
                    mapargs=['value', ], col_order=COL_ODR,
                    setargs=dict(xlim=(0.), xlabel='diameter')
                    )
        plv3 = plfacet(plt_type='distplot', **outkws2)
        plv3.plt(**set3)
        plv3.save_figure(op.join(savefolder, 'diameters.png'))

        set4 = dict(data=size,
                    col_wrap=4, col='media', hue='variable',
                    mapargs=['value', ], col_order=COL_ODR,
                    setargs=dict(xlim=(0.),  xlabel='volume')
                    )
        plv4 = plfacet(plt_type='distplot', **outkws2)
        plv4.plt(**set4)
        plv4.save_figure(op.join(savefolder, 'size.png'))

        set5 = dict(data=momdy,
                    col_wrap=4, col='media', hue='media',
                    sharex=True, sharey=True, col_order=COL_ODR,
                    ylim=tuple(buddy[r'$\Delta\Psi$ scaled gradient']
                               .quantile([0.025, 0.975])),
                    mapargs=['mom axis position',
                             r'$\Delta\Psi$ scaled gradient'],
                    )
        plv5 = plfacet(plt_type='pointplot')
        plv5.plt(**set5)
        plv5.save_figure(op.join(savefolder, 'momDY.png'))

        set6 = dict(data=buddy,
                    col_wrap=4, col='media', hue='media',
                    sharex=True, sharey=True, col_order=COL_ODR,
                    ylim=tuple(buddy[r'$\Delta\Psi$ scaled gradient']
                               .quantile([0.025, 0.975])),
                    mapargs=['bud axis position',
                             r'$\Delta\Psi$ scaled gradient'],
                    )
        plv6 = plfacet(plt_type='pointplot')
        plv6.plt(**set6)
        plv6.save_figure(op.join(savefolder, 'budDY.png'))


        print "Finished plotting!"
        return 0

    except UsageError as e:
        print e
        return 1

if __name__=='__main__':
    sys.exit(main())
