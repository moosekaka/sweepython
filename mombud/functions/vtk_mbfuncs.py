# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 16:45:59 2016
Functions for mom bud analysis in module vtk_mom_bud_analyse.py
@author: sweel
"""
from collections import defaultdict
import pandas as pd
import numpy as np
from tvtk.api import tvtk
import vtk
# pylint: disable=C0103
vtkF = defaultdict(dict)
mombud = defaultdict(dict)


def vtkopen(fpath):
    """
    wrapper to open polydata files
    """
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(fpath)
    reader.Update()
    data = reader.GetOutput()
    return data


def cellpos(cellname, df, **kwargs):
    """
    Return DataFrame of cell along mom-bud axis coords.

    Parameters
    ----------
    cellname : str
               Name of cell
    df : dataFrame
         Dataframe of mom,bud,neck coords

    kwargs:
    -------
    dyscale : str
        `DY_minmax` (default)

    dyraw : str
        `DY_raw` (default)

    Returns
    -------
    celldf : DataFrame
        Columns `DY`, `x`, `wholecell_xaxis`, `type`, `indcell_xaxis`

    """
    dyscale = kwargs.pop("dyscale", "DY_minmax")
    dyraw = kwargs.pop("dyraw", "DY_raw")

    cellkey = cellname.rsplit('\\', 1)[1][:-4]
    data = vtkopen(cellname)
    data = tvtk.to_tvtk(data)
    # is a column vec of R^3 (coordinates in the skel)
    npx = data.points.to_array()
    # indices of npx that would sort npx according to the x-axis
    xind = npx[:, 0].argsort()
    dy = data.point_data.get_array(dyscale).to_array()
    dy_raw = data.point_data.get_array(dyraw).to_array()

    #  individual skeletons xyz and Δψ
    celldf = pd.DataFrame({'x': npx[:, 0][xind],
                           'DY': dy[xind],
                           'DY_abs': dy_raw[xind]})
    xn, _, _ = df.ix[cellkey, 'neck']
    xb, _, _ = df.ix[cellkey, 'base']
    xt, _, _ = df.ix[cellkey, 'tip']

    xn_scaled = (xn - celldf.ix[0, 'x']) / (xt - xb)

    celldf['whole_cell_axis'] = ((celldf.ix[:, 'x'] -
                                  celldf.ix[0, 'x']) / (xt - xb))
    celldf['type'] = ''
    celldf.loc[celldf.x > xn, ['type']] = 'bud'
    celldf.loc[celldf.x <= xn, ['type']] = 'mom'
    #
    celldf.loc[celldf.type == 'bud',
               'ind_cell_axis'] = (celldf.ix[:, 'x']-xn) / (xt-xn)
    celldf.loc[celldf.type ==
               'mom', 'ind_cell_axis'] = (celldf.ix[:, 'x']-xb) / (xn-xb)
    celldf.index.name = cellkey
    return dict(df=celldf, neckpos=xn, neckpos_s=xn_scaled)


def bincell(cellname, col, bins):
    """
    Return a cell DataFrame  and `col` according to `bins`.

    Parameters
    ----------
    cellname : str
           name of cell
    col  : str
           col to be binned
    bins : list
           sequence of scalars for the bin edges

    Returns
    -------
    column of categorical labels for the bins
    """
    binnedcell = pd.cut(cellname.ix[:, col],
                        bins,
                        labels=bins[1:])

    return binnedcell


def mombudscale(df, fname, cellmean):
    """
    create mom/bud Series of average Δψ, scaled to cell minmax values,
    DY_minmax on MOM-BUD cell AXIS.
    """
    D = dict()
    for t in ['mom', 'bud']:
        D[t] = df.ix[df['type'] == t].groupby('ind_cell_binpos').DY.mean()
        D[t] = D[t] / cellmean
        D[t].name = fname
        # scale by minmax of whole cell
        # D['mom'] = lambda x, y :(x-y.min()) / (x.max()-y.min())
        # Xbud = vf.scaleminmax(Xbud, scaled_dy_wholecell)
    return D


def neckDY(fname, celldf, neck_position, outdic, dist=None):
    """
    Return two Series of points within a range of +-dist from neck
    """
    if dist is None:
        dist = [.15, .3, .5]
    tempdic = defaultdict(dict)
    tempdic[fname]['bud'] = {}
    tempdic[fname]['mom'] = {}

    for d in dist:
        tempdic[fname]['bud'][d] = celldf.loc[(celldf.x >= neck_position) &
                                              (celldf.x <
                                              (neck_position + d))].DY.mean()
        tempdic[fname]['mom'][d] = celldf.loc[(celldf.x < neck_position) &
                                              (celldf.x >=
                                              (neck_position - d))].DY.mean()
    return outdic.update(tempdic)


def xcell(x, f):
    """
    return DataFrame of mom binned Δψ + first point of bud
    """
    x['temp'] = x.index.astype('float')

    if len(f.DY.values):
        x = x.append(pd.Series({'DY': f.get_value(0, 'DY')}),
                     ignore_index=True)
    else:
        x = x.append(pd.Series({'DY': 0}),
                     ignore_index=True)

    x.loc[max(x.index), 'temp'] = 'fp'
    x['cellaxis_mom_budfp'] = x.temp.astype('category', ordered=True)
    return x
