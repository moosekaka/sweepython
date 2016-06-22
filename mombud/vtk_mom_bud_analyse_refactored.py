# -*- coding: utf-8 -*-
"""
Main module to analyze mom bud asymmetry
"""
import os
import os.path as op
import cPickle as pickle
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mombud.functions import vtk_mbfuncs as vf
from mombud.functions import vtk_mbplots as vp
import wrappers as wr
# pylint: disable=C0103


class UsageError(Exception):
    """
    Class for user-facing (non-programming) errors
    """
    pass


def getdata():
    """
    Get input data from specified work dirs

    Returns
    -------
    filekeys_f, datadir : dict(Str), Str
        filepaths of inputs and folder for data files
    dfmb : DataFrame
        cell volume data
    """

    datadir = op.join(os.getcwd(), 'mutants', 'transformedData', 'filtered')

    # old data
    datadir_old = op.join(os.getcwd(), 'data', 'transformedData')

    # DataFrames for new and old cell picked point
    with open(op.join(datadir, 'mombudtrans_new.pkl'), 'rb') as inpt:
        dfmb = pickle.load(inpt)  # columns base, neck, tip, media, bud, mom

    with open(op.join(datadir_old, 'mombudtrans.pkl'), 'rb') as inpt:
        dfmb_o = pickle.load(inpt)  # columns base, neck, tip, media, bud, mom

    # reject candidates
    rejectfold = op.join(datadir, os.pardir, 'reject')
    reject = wr.swalk(rejectfold, '*png', stop=-4)

    # VTK files for new and old data
    try:
        filext = "*vtk"
        vtkF = wr.ddwalk(datadir, filext, stop=-4)
    except:
        raise UsageError(
            "filetypes {} not found in {}".format(filext, datadir))

    try:
        filext = "*vtk"
        vtkF_old = wr.swalk(datadir_old, filext, stop=-4)
    except:
        raise UsageError(
            "filetypes {} not found in {}".format(filext, datadir))

    # file paths for VTKs
    filekeys_old = {item: vtkF_old[item] for item
                    in sorted(vtkF_old.keys()) if
                    item.split('_')[0] != 'YPD' and
                    item not in reject}

    filekeys = {item: vtkF[media][item] for media
                in sorted(vtkF.keys()) for item
                in sorted(vtkF[media].keys())}

    filekeys_f = {f: filekeys[f] for f in filekeys if f not in reject}
    filekeys_f.update(filekeys_old)  # add YPE to dict
    # dataframe of neck, mom and bud tip positions, bud and mom volumes
    dfmb = dfmb.append(dfmb_o[dfmb_o.media != 'YPD'])
    return filekeys_f, dfmb, datadir


def process_ind_df(vtkdf, mbax=None, cellax=None, **kwargs):
    """
    bin Δψ distrbution along cellaxis for each ind. cell DataFrame and append
    to the respective DataFrames

    Parameters
    ----------
    vtkdf : dict
        dict of DataFrames of ind. cell data
    mbax : np.array
        bin cell position for ind. mom/bud cell
    cellax : np. array
        bin cell position for whole cell
    Returns
    -------
    dicout : dict
        dictionary of DataFrames for mom bud analyses
    """

    if mbax is None:
        raise UsageError('please specify bin range for mom bud axis')
    if cellax is None:
        raise UsageError('please specify bin range for whole cell axis')

    # Dicts for budding progression and budratio DataFrames, etc.
    dicout = defaultdict(dict)
    dicint = defaultdict(dict)
    keys = sorted(vtkdf.keys())

    for k in keys:
        cell = vtkdf[k]['df']  # cell is ref/view, modification is inplace
        # update with whole cell stat. data
        dicint['cell'][k] = vtkdf[k]['celldata']
        # set index to cellname so that concatenate is possible
        cell['name'] = k
        cell.set_index('name', inplace=True)

    # Concat of ALL cell dfs into one giant DataFrame
    df_concat = pd.concat([vtkdf[k]['df'] for k in keys])
    df_concat.reset_index(inplace=True)
    groups = df_concat.groupby('name')
    # bin the dataframe according to individual (mom/bud) axis
    df_concat['ind_cell_binpos'] = (groups['ind_cell_axis']
                                    .apply(pd.cut,
                                           bins=mbax, labels=mbax[1:]))
    # bin the dataframe according to individual entire cell axis
    df_concat['whole_cell_binpos'] = (groups['whole_cell_axis']
                                      .apply(pd.cut,
                                             bins=cellax, labels=cellax[1:]))

    # groupby mom/buds , get agg. stats for Δψ
    df_agg = (df_concat.groupby(['name', 'type'])[['DY', 'DY_abs', 'DY_unscl']]
              .agg([np.mean, np.median]).unstack())
    df_agg.columns = (['index'] +
                      ['_'.join(c) for c in df_agg.columns.values[1:]])

    # DataFrame for all ind. cells
    dicout['dfcell'] = pd.DataFrame.from_dict(dicint['cell'], orient='index')
    dicout['dfcell'] = dicout['dfcell'].merge(df_agg,
                                              left_index=True,
                                              right_index=True)

    # Scaling factor for raw GFP daily variations
#    df2 = dicout['dfcell'].groupby(['date']).whole_cell_abs.mean()
#    df3 = df2 / df2.min()
#    dicout['dfcell'] = dicout['dfcell'].join(df3, on='date',
#                                             lsuffix='_org',
#                                             rsuffix='_scaling_fact')

    # bin by ind cell position and scale by whole cell mean
    dfbinned = (df_concat.groupby(['name', 'type', 'ind_cell_binpos']).
                DY.mean().unstack(level='ind_cell_binpos'))
    dfbinned.columns = dfbinned.columns.astype('float')
    df = dicout['dfcell']['whole_cell_mean']  # whole cell mean Δψ
    for i in ['dfbud', 'dfmom']:
        dicout[i] = dfbinned.xs(i[2:], level='type')
        dicout[i] = dicout[i].div(df, axis=0)   # scaling by whole cell mean Δψ
    return dicout


def postprocess_df(**kwargs):
    """
    Set population level data ,update parameters dict for plotting and filter
    unwanted data

    Returns
    -------
    outputdic : dict
        dictionary of DataFrames for mom bud analyses
    """

    kwargs['filekeys_f'], kwargs['dfmb'], kwargs['savefolder'] = getdata()
    Dout = process_ind_df(vf.gen_data(dfvoldata=kwargs['dfmb'],
                                      fkeys=kwargs['filekeys_f'],
                                      **kwargs), **kwargs)
    cellall = Dout['dfcell']
    cellposmom = Dout['dfmom']
    cellposbud = Dout['dfbud']

    # add cell volume data from cell tracing data
    cellall['budvol'] = kwargs['dfmb'].bud
    cellall['momvol'] = kwargs['dfmb'].mom

    # strip 'c' from some dates
    cellall['date'] = cellall.date.apply(
        lambda x: x.replace(x[0], '0') if x.startswith('c') else x)

    # YPE subdataset
    YPE = cellall[(cellall.type == 'YPE') | (cellall.type == 'WT')]
    YPE = YPE.reset_index(drop=True)

    # ratio of mom/bud Δψ
    cellall['frac'] = (cellall.loc[:, 'DY_median_bud'] /
                       cellall.loc[:, 'DY_median_mom'])

    #  90th percentile bud volume of each media type
    q90 = cellall.groupby('type').budvol.quantile(.90)
    cellall = cellall.join(q90, on='type', rsuffix='_q90')
    #  budvolratio is based on the largest 10% cells
    cellall['budvolratio'] = cellall.budvol / cellall.budvol_q90

    # concatenate type and volume data from cellall DataFrame to mom and bud df
    dic = dict(zip(['bud', 'mom'], [cellposbud, cellposmom]))
    for key, val in dic.iteritems():
        dic[key] = pd.concat([val, cellall.loc[:, ['type']]], axis=1)
        dic[key]['%svol' % key] = cellall['%svol' % key]
        dic[key]['binvol'] = vf.bincell(
            dic[key], '%svol' % key, kwargs['binsvol%s' % key])

    # Add bins used for plotting budding progression
    # add the 2. cat. for cells that are larger than the 90th percentile
    binsaxisbig = kwargs['cellax']
    cellall['bin_budprog'] = vf.bincell(cellall,
                                        'budvolratio',
                                        np.r_[binsaxisbig, [2.]])
    cellall['binbudvol'] = dic['bud']['binvol']

    # Filter conditions, reject large cells
    filt_size = (cellall.momvol > 100) & (cellall.type != 'YPD')
#    filt_type = (
#                  ((cellall.type == 'YPE') & (cellall.date == '052315')) |
#                  ((cellall.type == 'WT') & (cellall.date == '032716'))
#                 )
    mask = ~(filt_size)
#    mask = ~(filt_size) & ~(filt_type)  # (opt. reject YPE)
    cellall = cellall[mask]

    # Output dict labels
    outputdic = {'data': cellall, 'data_ype': YPE}
    for i in ['dfmom', 'dfbud']:
        outputdic[i] = dic[i[2:]][mask]
    outputdic['counts'] = cellall.groupby('type').size().to_dict()  # type
    outputdic['counts_ype'] = YPE.groupby('date').size().to_dict()
    outputdic['counts_buds'] = cellall.groupby(['type', 'binbudvol']).size()
    outputdic['counts_date'] = cellall.groupby('date').size().to_dict()
    outputdic.update(kwargs)
    return outputdic


def main(**kwargs):
    """
    Main

    kwargs
    ------
    plotlist : List
        plotting function names
    """

    wd = os.path.expanduser(os.sep.join(
        ('~', 'Documents', 'Github', 'sweepython', 'WorkingData')))
    os.chdir(wd)
    def_args = {'regen': False,
                'save': False,  # toggle to save plots
                'inpdatpath': 'celldata.pkl',
                'mbax': np.linspace(0., 1., 6),  # pos. along mom/bud cell
                'cellax': np.linspace(0, 1., 11),  # position along whole cell
                'binsvolbud': np.linspace(0, 40, 5),  # vol binning for bud
                'binsvolmom': np.array([0, 30, 40, 80.]),
                'COL_ODR': ['MFB1', 'NUM1', 'YPT11',
                            'WT', 'YPE', 'YPL', 'YPR'],
                'HUE_ODR': ['DY_abs_mean_mom',
                            'DY_abs_mean_bud',
                            'whole_cell_abs']}
    def_args.update(kwargs)  # override default args with user kwargs, if any

    plot_list = kwargs.get('plotlist',
                           ['plotDyAxisDist',  # 0
                            'plotSizeDist',  # 1
                            'plotBudProgr',  # 2
                            'plotGFP',     # 3
                            'plotViolins',  # 4
                            'plotRegr',
                            'plotDims'])

    outputargs = postprocess_df(**def_args)  # call getdata(), process_ind_df()
    # plots
    for f in plot_list:
        getattr(vp, f)(**outputargs)
        print 'finished {}'.format(f)
# _____________________________________________________________________________
if __name__ == '__main__':
    plt.close('all')
    L = ('plotDyAxisDist',  # 0
         'plotSizeDist',  # 1
         'plotBudProgr',  # 2
         'plotGFP',     # 3
         'plotViolins',  # 4
         'plotRegr',
         'plotDims')
    # labs == first two letters after plotXXX
    labs = (l.lower().partition('plot')[2][:2] for l in L)
    D = dict(zip(labs, L))
    main(regen=False, plotlist=[D['di']], save=True)
#    main(plotlist=D.values()[1:-1], save=False)
#    main()
