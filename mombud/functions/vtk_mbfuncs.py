# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 16:45:59 2016
Functions for mom bud analysis in module vtk_mom_bud_analyse.py
@author: sweel
"""
import os.path as op
from collections import defaultdict
import cPickle as pickle
import pandas as pd
from tvtk.api import tvtk
import vtk
import numpy as np
# pylint: disable=C0103


class UsageError(Exception):
    """
    Class for user-facing (non-programming) errors
    """
    pass


def gen_data(regen=False, **kwargs):
    """
    wrapper func to call mungedata, pass default params in kwargs and
    regenerate individual vtk DataFrames via vf.cellpos()

    Kwargs
    ------

    inpdatpath : Str
        filepath to celldata pickle file, if not specified then `dfvoldata` and
        `fkeys` must be specified

    dfvoldata : DataFrame
        cell volume data

    fkeys : dict
        dictionary of filepaths to individual cell VTK data

    Returns
    -------

    F : dict
        dictionary of DataFrames data for individual cells, output of calling
        cellpos()
    """

    fpath = kwargs.get('inpdatpath')

    # regenerates pickle file if not exist
    if regen or not op.isfile(fpath):
        F = {}

        for k in ['dfvoldata', 'fkeys']:
            if k not in kwargs:
                raise UsageError('Missing {}'.format(k))
        dfvol = kwargs.get('dfvoldata')
        filepaths = kwargs.get('fkeys')

        for k in sorted(filepaths):
            F[k] = cellpos(filepaths[k], dfvol)
        with open(fpath, 'wb') as out:
            pickle.dump(F, out)
    else:
        with open(fpath, 'rb') as inp:
            F = pickle.load(inp)
    return F


def vtkopen(fpath):
    """
    wrapper to open polydata files
    """
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(fpath)
    reader.Update()
    data = reader.GetOutput()
    return data


def gen_cell_dict(data, label=None, vtk_label=None):
    """
    Create DataFrame of scalar values according to vtk_label[label]
    """
    npx = data.points.to_array()
    # indices of npx that would sort npx according to the x-axis
    xind = npx[:, 0].argsort()
    labels = zip(label, vtk_label)
    def_dic = defaultdict(dict)
    for h, i in labels:
        def_dic[h][i] = data.point_data.get_array(i).to_array()

    out_dic = ({h: def_dic[h][j][xind]
                for h, i in def_dic.iteritems()
                for j in i.keys()})

    out_dic.update({'x': npx[:, 0][xind]})

    return pd.DataFrame(out_dic)


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
    outdic = {}
    cellkey = cellname.rsplit('\\', 1)[1][:-4]
    data = vtkopen(cellname)
    data = tvtk.to_tvtk(data)

    # labels for types of scalar values to select from VTK file
    lab = {'label': ['DY', 'DY_abs', 'DY_unscl'],
           'vtk_label': ['DY_minmax', 'bkstGFP', 'DY_raw']}

    # Main DataFrame
    celldf = gen_cell_dict(data, **lab)
    # Cell picked points
    xn, xb, xt = [x[0] for x in df.loc[cellkey, ['neck', 'base', 'tip']]]

    celldf['type'] = np.where(celldf['x'] > xn, 'bud', 'mom')
    first_pos = celldf.loc[0, 'x']
    outdic['neckpos_scaled'] = (xn - first_pos) / (xt - xb)
    celldf['whole_cell_axis'] = (celldf.loc[:, 'x'] - first_pos) / (xt - xb)

    # calc. individual cell position grouped by bud or mom
    celldf['ind_cell_axis'] = (celldf.x - xn) / (xt - xn)  # default for buds
    (celldf['ind_cell_axis']
     .where(celldf.type == 'bud', (celldf.x - xb) / (xn - xb), inplace=True))

    celldf.index.name = cellkey
    outdic['bud_diameter'] = xt - xn
    outdic['mom_diameter'] = xn - xb
    outdic['neckpos'] = xn

    return dict(df=celldf, celldata=outdic)


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
