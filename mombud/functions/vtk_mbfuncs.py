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
import traceback
from wrappers import UsageError
# pylint: disable=C0103
class FalseException(Exception):
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
    F = {}
    try:  # open the pickle file if exist
        if regen:
            raise FalseException
        fpath = kwargs['inpdatpath']
        with open(fpath, 'rb') as inp:
            F = pickle.load(inp)

    except (FalseException, AssertionError, KeyError, IOError):

        try:  # try to regen the data
            dfvol = kwargs['dfvoldata']
            filepaths = kwargs['fkeys']
            for k in sorted(filepaths):
                F[k] = cellpos(filepaths[k], dfvol)
            with open('celldata.pkl', 'wb') as out:
                pickle.dump(F, out)

        except KeyError:  # not enough input to regen data
            traceback.print_stack(limit=4)
            raise UsageError("Must have 'dfvoldata' and 'fkeys' kwargs")

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
    maskMom = makeMask('mom')
    maskBud = makeMask('bud')
    for d in dist:
        tempdic[fname]['bud'][d] = celldf.loc[
            maskBud(celldf, neck_position, d)].DY.mean()
        tempdic[fname]['mom'][d] = celldf.loc[
            maskMom(celldf, neck_position, d)].DY.mean()
    return outdic.update(tempdic)


def makeMask(swtch):
    """
    Make mask according to 'mom' or 'bud' region
    """
    if swtch not in ['mom', 'bud']:
        raise ValueError('must be either "mom" or "bud"')

    else:
        def _mask_fn(df, x1, rad):
            if swtch == 'mom':
                return (df.x < x1) & (df.x >= (x1 - rad))
            elif swtch == 'bud':
                return (df.x >= x1) & (df.x < (x1 + rad))
        return _mask_fn
