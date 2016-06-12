# -*- coding: utf-8 -*-
"""
Main module to analyze mom bud asymmetry
"""
import sys
import os
import os.path as op
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mombud.functions import vtk_mbfuncs as vf
from mombud.functions import vtk_mbplots as vp
import wrappers as wr
# pylint: disable=C0103


def mungedata(filepaths, df, **kwargs):
    """
    compute Δψ distrbution along cellaxis for each ind. cell and append
    to the resepective DataFrames
    """
    # bins for binning the bud progression ratio
    mbax = kwargs.pop('binsaxis', 0)
    cellax = kwargs.pop('binsaxisbig', 0)

    # DataFrames for budding progression and budratio, size distr., frac Δψ
    dfcell = pd.DataFrame()
    dfmom = pd.DataFrame()  # for Δψ distributino along mom/bud cell axis
    dfbud = pd.DataFrame()
#    cellposmom_fp = pd.DataFrame()  # for Δψ distributino mom + fp of bud
    dfneck = pd.DataFrame()  # for average Δψ around the neck region
    for k in sorted(filepaths):

        # returns Dataframe of pos along x-axis for inidivual mom/bud cells
        cell = vf.cellpos(filekeys_f[k], df)

        # bin the dataframe according to individual (mom/bud) axis
        cell['ind_cell_binpos'] = vf.bincell(cell, 'ind_cell_axis', mbax)


        # bin the dataframe according to individual entire cell axis
        cell['whole_cell_binpos'] = vf.bincell(cell,
                                               'whole_cell_axis',
                                               cellax)
    #  first point on bud
#        fp = cell[cell['type'] == 'bud'][:1].reset_index(drop=True)
        Xcell = cell.groupby('whole_cell_binpos').DY.mean()
#        Xcell = Xcell[Xcell < cell.neckpos_cellaxis.max()].reset_index()
        Xcell[Xcell > cell.neckpos_cellaxis.max()] = np.nan
        Xcell = Xcell.reset_index()
#        xc = vf.xcell(Xcell, fp)
#        cellposmom_fp = cellposmom_fp.append(xc)  # momcell + first bud point

        # DY Series data to DataFrames
        DYseries, dy_wholecell = vf.dyseries(cell, k)
        dfcell = dfcell.append(DYseries)

        # DY mom - bud axis data
        X = vf.mombudscale(cell, k, dy_wholecell.whole_cell_mean)
        dfbud = dfbud.append(X['bud'])
        dfmom = dfmom.append(X['mom'])

        # neckregion analy.
        neckreg = vf.neckDY(k, cell)
        dfneck = dfneck.append(neckreg, ignore_index=False)
    return dfcell, dfmom, dfbud, dfneck

# _____________________________________________________________________________
if __name__ == '__main__':
    plt.close('all')
    bins = {'binsaxis': np.linspace(0., 1., 6),  # pos. along mom/bud cell
            'binsaxisbig': np.linspace(0, 1., 11),  # position along whole cell
            'binsvolbud': np.linspace(0, 40, 5),  # vol binning for bud
            'binsvolmom': np.array([0, 30, 40, 80.])}  # vol binning for mom

# =============================================================================
#     Data Input
# =============================================================================
    # new data
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
        sys.exit("error filetypes %s not found in %s" % (filext, datadir))

    try:
        filext = "*vtk"
        vtkF_old = wr.swalk(datadir_old, filext, stop=-4)
    except LookupError:
        sys.exit("error filetypes %s not found in %s" % (filext, datadir))

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

# =============================================================================
#    Call mungedata()
# =============================================================================
    cellall, cellposmom, cellposbud, neckregion = mungedata(filekeys_f, dfmb,
                                                            **bins)

# =============================================================================
#    cleanup and add. labels for dataframes, calculate aggr measures etc.
# =============================================================================
    cellall['budvol'] = dfmb.bud
    cellall['momvol'] = dfmb.mom
    for i in [cellall, cellposbud, cellposmom]:
        i.reset_index(inplace=True)
    cellall['type'] = cellall['index'].apply(lambda x: x.split('_')[0])
    cellall['date'] = cellall['index'].apply(lambda x: x.split('_')[1])

    # strip 'c' from some dates
    stripc = lambda x: x.replace(x[0], '0') if x.startswith('c') else x
    cellall['date'] = cellall.date.apply(lambda x: stripc(x))
    # YPE subdataset
    YPE = cellall[(cellall.type == 'YPE') | (cellall.type == 'WT')]
    YPE = YPE.reset_index(drop=True)

    # remove high and low DY day
    cellall = cellall.ix[~(((cellall.type == 'YPE') &
                            (cellall.date == '052315')) |
                           ((cellall.type == 'WT') &
                            (cellall.date == '032716')))]

#    cellposmom_fp = cellposmom_fp.reset_index()
#    cellposmom_fp = pd.concat(
#        [cellposmom_fp, cellall.ix[:, ['type']]], axis=1)

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
            dic[key], '%svol' % key, bins['binsvol%s' % key])

    # Add bins used for plotting budding progression
    # add the 2. cat. for cells that are larger than the 90th percentile
    binsaxisbig = bins['binsaxisbig']
    cellall['bin_budprog'] = vf.bincell(cellall,
                                        'budvolratio',
                                        np.r_[binsaxisbig, [2.]])
    cellall['binbudvol'] = dic['bud']['binvol']

    # reject super large cells
    rejectlist = dic['mom'].ix[(np.asarray(dic['mom'].momvol) > 100) &
                               (dic['mom'].type != 'YPD'), 'index']
    cellall = cellall.ix[~ cellall.ix[:, 'index'].isin(rejectlist)]
    cellposmom = dic['mom'].ix[~dic['mom'].ix[:, 'index'].isin(rejectlist)]
    cellposbud = dic['bud'].ix[~dic['bud'].ix[:, 'index'].isin(rejectlist)]

    # labels for counts of each condition, facet
    Nype = YPE.groupby('date').size().to_dict()
    Nbud = cellall.groupby(['type', 'binbudvol']).size()
    N = cellall.groupby('type').size().to_dict()  # counts for each type

# =============================================================================
#     Plot using seaborn
# =============================================================================
    plotfuncs = ['plotSizeDist',
                 'plotBudProgr',
                 'plotDyAxisDist',
                 'plotNeck',
                 'plotYPE',
                 'plotRegr']
    params = {'data': cellall, 'neckdata': neckregion, 'data_ype': YPE,
              'counts': N, 'counts_buds': Nbud, 'counts_ype': Nype,
              'dfmom': cellposmom, 'dfbud': cellposbud,
              'savefolder': datadir, 'save': False}
    params.update(bins)
    plt.close('all')

    for f in plotfuncs[:-1]:
        getattr(vp, f)(**params)
