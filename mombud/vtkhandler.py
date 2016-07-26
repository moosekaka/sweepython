# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 16:10:16 2016
Module for plots of analysis of mother bud function in budding yeast
"""
import sys
import os.path as op
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from wrappers import UsageError
import mombud.mungedf as munge
import mombud.functions.vtk_mbplots as mbfuncs

labelhandler, plviol, plbox, plfacet = (mbfuncs.labelhandler,
                                        mbfuncs.plviol,
                                        mbfuncs.plbox,
                                        mbfuncs.plfacet)
COL_ODR = munge.COL_ODR
HUE_ODR = munge.HUE_ODR
savefolder = munge.datadir
gfp_plot_vars =  ['DY_abs_mean_mom',
                  'DY_abs_mean_bud',
                  'DY_abs_cell_mean']
mombud_dy_vars = ['DY_median_mom', 'DY_median_bud']

def main():
    """
    Main
    """
    plt.close('all')
    try:
        labelFacet = labelhandler('facet')
        labelNormal = labelhandler()
        labelRowFacet = labelhandler('rowfacet')

        # args for postprocess_df(), which calls inputdata(), mungedata()
        inp_args = {'inpdatpath': 'celldata.pkl',
                    'mbax': np.linspace(0., 1., 6),  # pos. mom/bud cell
                    'cellax': np.linspace(0, 1., 11),  # whole cell pos.
                    'binsvolbud': np.linspace(0, 40, 5),  # vol bins for bud
                    'binsvolmom': np.array([0, 30, 40, 80.]),
                    }
        outputargs = munge.postprocess_df(**inp_args)
        df = outputargs.get('data')
        dfmom = outputargs.get('dfmom')
        dfbud= outputargs.get('dfbud')
#==============================================================================
# Data long form
#==============================================================================
        frac = pd.melt(df, id_vars=['media'],
                       value_vars=['frac'])
        frac_cellcycle =  df[['media','bin_budprog', 'frac']].reset_index(drop=True)
        mb_dy = pd.melt(df, id_vars=['media', 'date'],
                        value_vars=mombud_dy_vars)
        mb_gfp = pd.melt(df, id_vars=['media', 'date'],
                        value_vars=gfp_plot_vars)
        diameters = pd.melt(df, id_vars='media',
                            value_vars=['bud_diameter', 'mom_diameter'])
        size = pd.melt(df.loc[:, ['media','budvol', 'momvol']],'media')
        momdy = pd.melt(dfmom,
                        id_vars=['media', 'binvol'],
                        var_name='mom axis position',
                        value_name=r'$\Delta\Psi$ scaled',
                        value_vars=outputargs['mbax'].tolist())
#        momdy = momdy.dropna()
        buddy= pd.melt(dfbud,
                       id_vars=['media', 'binvol'],
                       var_name='bud axis position',
                       value_name=r'$\Delta\Psi$ scaled',
                       value_vars=outputargs['mbax'].tolist())
#        buddy = buddy.dropna()
    # mutant only subset
        mutants = (mb_dy[mb_dy['media']
                   .isin(['NUM1', 'MFB1', 'YPT11', 'WT'])]
                   .reset_index(drop=True))
        group = mutants.groupby('date')
#==============================================================================
        outkws1 = dict(default_ylims=[0.05, 0.95],
                       labeller=labelNormal, col_order=COL_ODR)

        outkws2 = dict(default_ylims=[0.15, 0.9],
                       labeller=labelFacet, col_order=COL_ODR)

        seta = dict(x='media', y='value',
                    hue='variable', group_key='media',
                    title='gfp uptake', ylim='auto',
                    )
        plva = plviol(**outkws1)
        plva.plt(data=mb_gfp, **seta)
        plva.turn_off_legend()
        plva.save_figure(op.join(savefolder, 'gfp_mombud.png'))

        set0 = dict(x='media', y='value',
                    hue='variable', group_key='media',
                    title='mombud', ylim=(0, 1.),
                    )
        plv = plviol(**outkws1)
        plv.plt(data=group.get_group('032016'), **set0)
        plv.turn_off_legend()
        plv.save_figure(op.join(savefolder, 'mutantsDY_mombud.png'))

        set1 = dict(x='media', y='value',
                    hue='variable', group_key='media',
                    title='mombud', ylim=(0, 1.),
                    )
        plv1 = plviol(**outkws1)
        plv1.plt(data=mb_dy, **set1)
        plv1.turn_off_legend()
        plv1.save_figure(op.join(savefolder, 'violin_mombud.png'))

        set2 = dict(x='media', y='value', ylim='auto',
                    hue='media', group_key='media', inner=None,
                    title='fracDY')
        plv2 = plbox(**outkws1)
        plv2.plt(data=frac, **set2)
        plv2.turn_off_legend()
        plv2.save_figure(op.join(savefolder, 'fracDY.png'))

        set3 = dict(col_wrap=4, col='media',
                    size=3, aspect=1.5, hue='variable',
                    col_order=COL_ODR,
                    setargs=dict(xlim=(0.), xlabel='diameter')
                    )
        plv3 = plfacet(plt_type='distplot', **outkws2)
        plv3.plt(data=diameters, mapargs=['value', ], **set3)
        plv3.save_figure(op.join(savefolder, 'diameters.png'))

        set4 = dict(col_wrap=4, col='media', hue='variable',
                    col_order=COL_ODR,
                    setargs=dict(xlim=(0.),  xlabel='volume')
                    )
        plv4 = plfacet(plt_type='distplot', **outkws2)
        plv4.plt(data=size, mapargs=['value', ], **set4)
        plv4.save_figure(op.join(savefolder, 'size.png'))

        set5 = dict(col_wrap=4, col='media', hue='media',
                    sharex=True, sharey=True, col_order=COL_ODR,
                    ylim=(0.0, 1.),
                    )
        plv5 = plfacet(plt_type='pointplot', **outkws2)
        plv5.plt(data=momdy,
                 mapargs=['mom axis position', r'$\Delta\Psi$ scaled'],
                 **set5)
        plv5.save_figure(op.join(savefolder, 'momDY.png'))

        set6 = dict(col_wrap=4, col='media', hue='media',
                    sharex=True, sharey=True, col_order=COL_ODR,
                    ylim=(0.0, 1.0),
                    )
        plv6 = plfacet(plt_type='pointplot', **outkws2)
        plv6.plt(data=buddy,
                 mapargs=['bud axis position', r'$\Delta\Psi$ scaled'],
                 **set6)
        plv6.save_figure(op.join(savefolder, 'budDY.png'))

        set7 = dict(x='bin_budprog', y='frac',
                    hue='media',
                    title=u"Δψ vs bud progression\n ", ylim=(0., 3.),
                    setargs=dict(xlabel="bud progression",
                                 ylabel=u"Δψ bud/Δψ mom",)
                    )
        plv7 = plviol(plt_type='pointplot', **outkws1)
        plv7.plt(data=frac_cellcycle, **set7)
        plv7.save_figure(op.join(savefolder, 'dy_cellcycle.png'))

        set8 = dict(col_wrap=4, col='media', hue='media',
                    sharex=True, sharey=True, col_order=COL_ODR,
                    title=u"Δψ vs bud progression\n ", ylim=(0., 3.),
                    setargs=dict(xlabel="bud progression",
                                 )
                    )
        plv8 = plfacet(plt_type='pointplot', **outkws2)
        plv8.plt(data=frac_cellcycle,
                 mapargs=['bin_budprog', 'frac'], **set8)

        print "Finished plotting!"
        return 0

    except UsageError as e:
        print e
        return 1

if __name__=='__main__':
    sys.exit(main())
