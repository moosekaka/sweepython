# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 01:24:48 2016
Bootstrap a distribution of mom bud positions
@author: sweel_rafelski
"""
import sys
import os
import os.path as op
import cPickle as pickle
import matplotlib.pyplot as plt
import scipy.stats as sp
import numpy as np
#from mayavi import mlab
import pandas as pd
import seaborn as sns
from tvtk.api import tvtk
from mombud.vtk_viz import vtk_mbfuncs as vf
import wrappers as wr
# pylint: disable=C0103
# pylint: disable=maybe-no-member
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.close('all')

#datadir = op.join(os.getcwd(), 'data', 'transformedData')
datadir = ("C:\\Users\\sweel_rafelski\\Documents\\GitHub\\sweepython\\WorkingData\\data\\transformedData")


def boot(X, p, N):
    """
    Bootstrap the cell values at column p by shuffling the indices

    Parameters
    ----------
    X : DataFrame
        Pandas dataframe to be sampled/bootstrapped
    p : str
        Column name to be shuffled
    N : Int
        Number of bootstrap runs
    Returns
    -------
    result :  list
        List of shuffled vals from column selected by p

    """
    result = []
    for _ in range(N):
        rX = X.ix[:, p].sample(n=len(X), replace=True)
#        idx = rX.index.tolist()
        result.append(rX.values)
    return result


def bstrp(X, N):
    n=0
    while n<N:
        yield np.random.permutation(X)
        np.random.randint


def cellpos(cellname, df, **kwargs):
    """
    Return DataFrame of cell dys along mom-bud axis coords.

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
    data = vf.vtkopen(cellname)
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

    celldf['posxcell'] = (celldf.ix[:, 'x'] - celldf.ix[0, 'x']) / (xt - xb)
    celldf['type'] = ''
    celldf.loc[celldf.x > xn, ['type']] = 'bud'
    celldf.loc[celldf.x <= xn, ['type']] = 'mom'
    celldf.ix[celldf.type == 'bud', 'posx'] = (celldf.ix[:, 'x']-xn) / (xt-xn)
    celldf.ix[celldf.type == 'mom', 'posx'] = (celldf.ix[:, 'x']-xb) / (xn-xb)
    celldf.reset_index(drop=True, inplace=True)
    celldf['neckpos'] = xn
    return celldf


if __name__ == '__main__':
    with open(op.join(datadir, 'mombudtrans.pkl'), 'rb') as inpt:
        dfmb = pickle.load(inpt)  # has cols base, neck, tip, media, bud, mom

    filext = "*vtk"
    try:
        vtkF = wr.ddwalk(datadir, filext, stop=-4)
    except:
        Exception
        sys.exit("error filetypes %s not found in %s" % (filext, datadir))

    filekeys = {item: vtkF[media][item] for media
                in sorted(vtkF.keys()) for item
                in sorted(vtkF[media].keys())}


    for key in sorted(filekeys)[:]:
        cell = vf.cellpos(filekeys[key], dfmb)
        # pos along x-axis for inidivual mom or bud cell
    #    cell['binposx'] = vf.bincell(cell, 'posx', binsaxis)
    #    # scaled Δψ to min-max of the GRADIENT within mom/bud
    #    Xmom = cell.ix[cell['type'] == 'mom'].groupby('binposx').DY.mean()
    #    Xbud = cell.ix[cell['type'] == 'bud'].groupby('binposx').DY.mean()
    #
    #    # pos along x-axis for the whole cell
    #    cell['binposxcell'] = vf.bincell(cell, 'posxcell', binsaxisbig)
    #    # scale Δψ to min-max of the GRADIENT within cell
    #    DYcell = cell.groupby('binposxcell').DY.mean()
    #
    #    Xmom = vf.scaleminmax(Xmom, DYcell)
    #    Xbud = vf.scaleminmax(Xbud, DYcell)
    #    Xbud.name = key
    #    Xmom.name = key
    #    medianDY = cell.groupby('type').median().DY
    #    medianDY.name = key
    #    cellall = cellall.append(medianDY)
    #    cellposbud = cellposbud.append(Xbud)
    #    cellposmom = cellposmom.append(Xmom)
    #    # temp dict of mean Δψ at +- range of dist from budneck
    #    tempdic = {dist: vf.neckDY(cell, cell.neckpos, dist)
    #               for dist in [.15, .3, .5]}
    #    temp = pd.DataFrame({'bud': pd.Series({k: tempdic[k][0] for k in tempdic}),
    #                         'mom': pd.Series({k: tempdic[k][1] for k in tempdic}),
    #                         'cellname': key})
    #    temp['dist'] = temp.index
    #    temp.set_index('cellname', inplace=True)
    #    neckregion = neckregion.append(temp, ignore_index=False)
