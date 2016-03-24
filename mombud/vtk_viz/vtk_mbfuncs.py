# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 16:45:59 2016
Functions for mom bud analysis in module vtk_mom_bud_analyse.py
@author: sweel
"""
from collections import defaultdict
import pandas as pd
from tvtk.api import tvtk
import vtk
# pylint: disable=C0103
# pylint: disable=maybe-no-member
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

    Returns
    -------
    celldf : DataFrame
        Columns `DY`, `x`, `posxcell`, `type`, `posx`, `pos`, `binposx`,
        `binposxcell`
    """
    cellkey = cellname.rsplit('\\', 1)[1][:-4]
    data = vtkopen(cellname)
    data = tvtk.to_tvtk(data)
    # is a column vec of R^3 (coordinates in the skel)
    npx = data.points.to_array()
    # indices of npx that would sort npx according to the x-axis
    xind = npx[:, 0].argsort()
    dy = data.point_data.get_array(u"DY_minmax").to_array()

    #  individual skeletons xyz and Δψ
    celldf = pd.DataFrame({'x': npx[:, 0][xind],
                           'DY': dy[xind]})
    xn, _, _ = df.ix[cellkey, 'neck']
    xb, _, _ = df.ix[cellkey, 'base']
    xt, _, _ = df.ix[cellkey, 'tip']

#    xn_scaled = (xn - cell.ix[0, 'x']) / (xt - xb)

    celldf['posxcell'] = (celldf.ix[:, 'x'] - celldf.ix[0, 'x']) / (xt - xb)
    celldf['type'] = ''
    celldf.loc[celldf.x > xn, ['type']] = 'bud'
    celldf.loc[celldf.x <= xn, ['type']] = 'mom'
    celldf.ix[celldf.type == 'bud', 'posx'] = (celldf.ix[:, 'x']-xn) / (xt-xn)
    celldf.ix[celldf.type == 'mom', 'posx'] = (celldf.ix[:, 'x']-xb) / (xn-xb)
    celldf.reset_index(drop=True, inplace=True)
    celldf['neckpos'] = xn

    return celldf


def neckDY(celldf, neckloc, dist=0.3):
    """
    Return two Series of points within a range of +-dist from neck
    """
    gtneck = celldf.loc[(celldf.x >= neckloc) & (celldf.x < (neckloc+dist))]
    ltneck = celldf.loc[(celldf.x < neckloc) & (celldf.x >= (neckloc-dist))]
    return gtneck.DY.mean(), ltneck.DY.mean()


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


def scaleminmax(cell1, cell2):
    """
    scaling of cell1 with cell2 min/max
    """
    scaled = (cell1-cell2.min()) / (cell2.max()-cell2.min())
    return scaled
