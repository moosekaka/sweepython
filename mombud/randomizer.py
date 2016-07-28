# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 21:43:21 2016

@author: sweel_Rafelski
"""
import os.path as op
import cPickle as pickle
import matplotlib.pyplot as plt
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
mutants=['MFB1', 'NUM1', 'WT', 'YPT11']
gfp_plot_vars =  ['DY_abs_mean_mom',
                  'DY_abs_mean_bud',
                ]
mombud_dy_vars = ['DY_median_mom', 'DY_median_bud']

kwargs = {}
def_args = {'regen':True,
            'save': False,
            'inpdatpath': 'celldata.pkl',
            'mbax': np.linspace(0., 1., 6),
            'cellax': np.linspace(0, 1., 11),
            'binsvolbud': np.linspace(0, 40, 5),
            'binsvolmom': np.array([0, 30, 40, 80.])}


def main(randomize=False, **kwargs):
    plt.close('all')
    def_args.update(kwargs)
    outputargs = mdf.postprocess_df(**def_args)
    df = outputargs['data']
    df_mutants = df[df.media.isin(mutants)]

    subset = {}
    for med in ['MFB1', 'NUM1', 'WT', 'YPT11']:
        subset[med] = df_mutants[(df_mutants['media'] == med) &
                                 (df_mutants['date'] == '071016')]

    try:
        with open('rowindex.pkl', 'rb') as inpt:
            rowind = pickle.load(inpt)
    except IOError:
        raise
#        if randomize:
#            rand_subset = {}
#        for i in subset.keys():
#            rand_subset[i] = subset[i].sample(frac=0.5)
#            changedates = [vals
#                           for key in rand_subset.keys()
#                           for vals in rand_subset[key].index.values]
#            rowind = df_mutants.index[df_mutants.index.isin(changedates)]
#            with open('rowindex.pkl', 'wb') as outp:
#                pickle.dump(rowind, outp)


    df_mutants.loc[rowind, 'date'] = '071116'
    df_mutants = df_mutants[~(df_mutants.date == '032716')]
    hiind = df_mutants.DY_abs_cell_mean.sort_values()
    df_mutants = df_mutants.loc[~(hiind>7000)]

    mb_dy_date = pd.melt(df_mutants, ['media', 'date'], mombud_dy_vars)
    mb_gfp = pd.melt(df_mutants, ['media', 'date'], gfp_plot_vars)
    mb_dy = pd.melt(df_mutants, ['media'], mombud_dy_vars)

    outkws = dict(default_ylims=[0.025, 0.975],
                  labeller=labelNormal, col_order=COL_ODR)
    set1 = dict(x='media', y='value',
                hue='variable', group_key='media',
                 ylim=(0,5000),
                )

    # mom bud gfp uptake
    for i in df_mutants.date.unique():
        plv = plviol(**outkws)
        plv.plt(data=mb_gfp[mb_gfp.date==i], title=i, **set1)
        plv.save_figure(op.join(mdf.datadir,
                                'gfpuptake_{}.png'.format(i)))

    #  mom bud dy
    set2 = dict(x='media', y='value',
                hue='variable', group_key='media',
                 ylim=(0, 1.),
    )
    for i in df_mutants.date.unique():
        plv = plviol(**outkws)
        plv.plt(data=mb_dy_date[mb_dy_date.date==i], title=i, **set2)
        plv.save_figure(op.join(mdf.datadir,
                                'mombud_dy_{}.png'.format(i)))

    plv = plviol(**outkws)
    plv.plt(data=mb_dy, title='mom-bud-DY', **set2)
    plv.save_figure(op.join(mdf.datadir,
                            'mombud_dy.png'.format(i)))

    df_mutants_and_old = pd.concat([df_mutants, df[~df.media.isin(mutants)]])
    with open('df_mombud_filtered.pkl', 'wb') as outp:
        pickle.dump(df_mutants_and_old, outp)

if __name__=='__main__':
    main(randomize=False)
