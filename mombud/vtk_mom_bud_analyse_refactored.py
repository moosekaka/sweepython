# -*- coding: utf-8 -*-
"""
Main module to analyze mom bud asymmetry
"""
import sys
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


def mungedata(vtkdf, **kwargs):
    """
    compute Δψ distrbution along cellaxis for each ind. cell and append
    to the resepective DataFrames
    """
    # bins for binning the bud progression ratio
    mbax = kwargs.get('binsaxis', 0)  # ind mom/bud cell bins
    cellax = kwargs.get('binsaxisbig', 0)  # whole cell bins

    # Dicts for budding progression and budratio DataFrames, etc.
    dicout = defaultdict(dict)
    dicout['dfmom_fp'] = pd.DataFrame()
    dicdf = defaultdict(dict)
    keys = sorted(vtkdf.keys())

    for k in keys:
        # get Dataframe of pos along x-axis for inidivual mom/bud cells
        cell = vtkdf[k]['df']

        # bin the dataframe according to individual (mom/bud) axis
        cell['ind_cell_binpos'] = vf.bincell(cell, 'ind_cell_axis', mbax)

        # bin the dataframe according to individual entire cell axis
        cell['whole_cell_binpos'] = vf.bincell(cell,
                                               'whole_cell_axis',
                                               cellax)

        #  DataFrame for mom bin points + first point on bud
#        fp = cell[cell['type'] == 'bud'][:1].reset_index(drop=True)
#        Xcell = cell.groupby('whole_cell_binpos').DY.mean()
#        Xcell[Xcell > cell.neckpos_cellaxis.max()] = np.nan
#        Xcell = Xcell.reset_index()
#        dfXcell = vf.xcell(Xcell, fp)
#        dfXcell['type'] = celltype
#        dfmom_fp = dfmom_fp.append(dfXcell)  # DataFrame output

        # Δψ Series data to DataFrames
        dy_wholecell = cell.mean()[['DY', 'DY_abs']]
        dy_wholecell.rename({'DY': 'whole_cell_mean',
                             'DY_abs': 'whole_cell_abs'}, inplace=True)
        dicdf['cell'][k] = dy_wholecell.to_dict()
        dicdf['cell'][k].update({'type': k.split('_')[0],
                                 'date': k.split('_')[1],
                                 'neckpos': vtkdf[k]['neckpos'],
                                 'neckpos_cellaxis': vtkdf[k]['neckpos_s']})

        # Δψ mom-bud axis data
        X = vf.mombudscale(cell, k, dy_wholecell.whole_cell_mean)
        dicdf['bud'][k] = X['bud']
        dicdf['mom'][k] = X['mom']

        # neckregion analy.
#        vf.neckDY(k, cell, vtkdf[k]['neckpos'], outdic=dicdf['neck'])

        # set common index for ind cell so that concatenate is possible
        vtkdf[k]['df']['name'] = k
        vtkdf[k]['df'].set_index('name', inplace=True)

    # get agg. stat. measures for concatenate DataFrame of all cells
    dftemp = pd.concat([vtkdf[k]['df'] for k in keys])
    dftemp = dftemp.reset_index()
    df_agg = (dftemp.groupby(['name', 'type'])[['DY', 'DY_abs']]
              .agg([np.mean, np.median]).unstack().reset_index())
    df_agg.columns = (['index'] +
                      ['_'.join(c) for c in df_agg.columns.values[1:]])

    dicout['dfcell'] = pd.DataFrame.from_dict(dicdf['cell'], orient='index')
    dicout['dfcell'] = dicout['dfcell'].merge(
        df_agg, left_index=True, right_on='index')
    dicout['dfbud'] = pd.DataFrame.from_dict(dicdf['bud'], orient='index')
    dicout['dfmom'] = pd.DataFrame.from_dict(dicdf['mom'], orient='index')

    return dicout


def inputdata():
    """
    Get input data from specified work dirs
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
    except LookupError:
        sys.exit("error filetypes {} not found in {}".format(filext, datadir))

    try:
        filext = "*vtk"
        vtkF_old = wr.swalk(datadir_old, filext, stop=-4)
    except LookupError:
        sys.exit("error filetypes {} not found in {}".format(filext, datadir))

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

# =============================================================================
#     Data Input
# =============================================================================
def setupDataFrames(**kwargs):
    """
    Set population level data and update parameters dict for plotting
    """
    kwargs['filekeys_f'], kwargs['dfmb'], kwargs['savefolder'] = inputdata()
    Dout = mungedata(vf.wrapper(**kwargs), **kwargs)
    cellall = Dout['dfcell']
    cellposmom = Dout['dfmom']
    cellposbud = Dout['dfbud']

    # cleanup and add. labels for dataframes, calculate aggr measures etc.
    cellall = cellall.set_index('index')
    cellall['budvol'] = kwargs['dfmb'].bud
    cellall['momvol'] = kwargs['dfmb'].mom
    for i in [cellall, cellposbud, cellposmom]:
        i.reset_index(inplace=True)

    # strip 'c' from some dates
    stripc = lambda x: x.replace(x[0], '0') if x.startswith('c') else x
    cellall['date'] = cellall.date.apply(lambda x: stripc(x))

    # YPE subdataset
    YPE = cellall[(cellall.type == 'YPE') | (cellall.type == 'WT')]
    YPE = YPE.reset_index(drop=True)

    # remove high and low Δψ day
#    cellall = cellall.ix[~(((cellall.type == 'YPE') &
#                            (cellall.date == '052315')) |
#                           ((cellall.type == 'WT') &
#                            (cellall.date == '032716')))]

    cellall['frac'] = (cellall.ix[:, 'DY_median_bud'] /
                       cellall.ix[:, 'DY_median_mom'])

    #  90th percentile bud volume of each media type
    Q = cellall.groupby('type').quantile(.90)
    cellall['q90'] = cellall.type.apply(lambda x: Q.ix[x].budvol)

    #  budvolratio is based on the largest 10% cells
    cellall['budvolratio'] = cellall.budvol / cellall.q90

    # concatenate type and volume data from cellall DataFrame to mom and bud df
    dic = dict(zip(['bud', 'mom'], [cellposbud, cellposmom]))
    for key, val in dic.iteritems():
        dic[key] = pd.concat([val, cellall.ix[:, ['type']]], axis=1)
        dic[key]['%svol' % key] = cellall['%svol' % key]
        dic[key]['binvol'] = vf.bincell(
            dic[key], '%svol' % key, kwargs['binsvol%s' % key])

    # Add bins used for plotting budding progression
    # add the 2. cat. for cells that are larger than the 90th percentile
    binsaxisbig = kwargs['binsaxisbig']
    cellall['bin_budprog'] = vf.bincell(cellall,
                                        'budvolratio',
                                        np.r_[binsaxisbig, [2.]])
    cellall['binbudvol'] = dic['bud']['binvol']

    # reject super large cells
    rejectlist = dic['mom'].ix[(np.asarray(dic['mom'].momvol) > 100) &
                               (dic['mom'].type != 'YPD'), 'index']
    cellall = cellall.ix[~ cellall.ix[:, 'index'].isin(rejectlist)]

    # Output dict labels
    outputdic = {'data': cellall,
                 'data_ype': YPE}
    outputdic['dfmom'] = (dic['mom'].
                          ix[~dic['mom'].ix[:, 'index'].isin(rejectlist)])
    outputdic['dfbud'] = (dic['bud'].
                          ix[~dic['bud'].ix[:, 'index'].isin(rejectlist)])
    outputdic['counts'] = cellall.groupby('type').size().to_dict()  # type
    outputdic['counts_ype'] = YPE.groupby('date').size().to_dict()
    outputdic['counts_buds'] = cellall.groupby(['type', 'binbudvol']).size()
    outputdic['counts_date'] = cellall.groupby('date').size().to_dict()
    outputdic.update(kwargs)
    return outputdic


def main(*plot_func_names, **kwargs):
    """
    Main
    """
    wd = os.path.expanduser(os.sep.join(
        ('~', 'Documents', 'Github', 'sweepython', 'WorkingData')))
    os.chdir(wd)
    args = {'regen': False,
            'save': False,  # toggle to save plots
            'inpdatpath': 'celldata.pkl',
            'binsaxis': np.linspace(0., 1., 6),  # pos. along mom/bud cell
            'binsaxisbig': np.linspace(0, 1., 11),  # position along whole cell
            'binsvolbud': np.linspace(0, 40, 5),  # vol binning for bud
            'binsvolmom': np.array([0, 30, 40, 80.])}

    # update defaults to  specified pars. in main()
    for arg in args:
        if arg in kwargs:
            args[arg] = kwargs[arg]

    outputargs = setupDataFrames(**args)

    # plots
    for f in plot_func_names:
        getattr(vp, f)(**outputargs)

# _____________________________________________________________________________
if __name__ == '__main__':
    plt.close('all')
    L = ['plotDyAxisDist',  # 0
         'plotSizeDist',  # 1
         'plotBudProgr',  # 2
         'plotGFP',     # 3
         'plotViolins',  # 4
         'plotRegr']
    D = dict(zip(range(len(L)), L))
    main(D[4], save=True)
