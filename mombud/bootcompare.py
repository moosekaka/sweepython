# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 01:24:48 2016
Bootstrap a distribution of mom bud positions
@author: sweel_rafelski
"""
import scipy.stats as sp
import sys
import os
import os.path as op
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
from mayavi import mlab
import pandas as pd
import seaborn as sns
from mombud.vtk_viz import vtk_mbfuncs as vf
import wrappers as wr

# pylint: disable=C0103
# pylint: disable=maybe-no-member
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.close('all')
#datadir = op.join( './data', 'transformedData')
datadir = op.join(os.getcwd(), 'data', 'transformedData')
#datadir = 'C:\\Users\\sweel_rafelski\\Documents\\GitHub\\sweepython\\WorkingData\\data\\transformedData'


def boot(X, p,  N):
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


with open(op.join(datadir, 'mombudtrans.pkl'), 'rb') as inpt:

    dfmb = pickle.load(inpt)  # has columns base, neck, tip, media, bud, mom

try:
    filext = "*vtk"
    vtkF = wr.ddwalk(datadir, filext, stop=-4)
except:
    Exception
    sys.exit("error filetypes %s not found in %s" % (filext, datadir))


filekeys = {item: vtkF[media][item] for media
            in sorted(vtkF.keys()) for item
            in sorted(vtkF[media].keys())}


# compute 'z' positions  dyRaw
cellall = pd.DataFrame(columns=['mom', 'bud'])
cellposmom = pd.DataFrame()
cellposbud = pd.DataFrame()
neckregion = pd.DataFrame()
binsbudprog = np.r_[np.arange(0, 1.1, .1), 2]
binsaxis = np.linspace(0., 1., 6)
binsaxisbig = np.linspace(0, 1., 11)
binsvolbud = np.linspace(0, 40, 5)
binsvolmom = np.array([0, 30, 40, 80.])


for key in sorted(filekeys)[:]:
    try:
        cell = wr.safecall(key, filekeys, dfmb)
    except Exception:
        sys.exit("%s not in filekeys" % key)
    # pos along x-axis for inidivual mom or bud cell
    cell['binposx'] = vf.bincell(cell, 'posx', binsaxis)
    # scaled Δψ to min-max of the GRADIENT within mom/bud
    Xmom = cell.ix[cell['type'] == 'mom'].groupby('binposx').DY.mean()
    Xbud = cell.ix[cell['type'] == 'bud'].groupby('binposx').DY.mean()

    # pos along x-axis for the whole cell
    cell['binposxcell'] = vf.bincell(cell, 'posxcell', binsaxisbig)
    # scale Δψ to min-max of the GRADIENT within cell
    DYcell = cell.groupby('binposxcell').DY.mean()

    Xmom = vf.scaleminmax(Xmom, DYcell)
    Xbud = vf.scaleminmax(Xbud, DYcell)
    Xbud.name = key
    Xmom.name = key
    medianDY = cell.groupby('type').median().DY
    medianDY.name = key
    cellall = cellall.append(medianDY)
    cellposbud = cellposbud.append(Xbud)
    cellposmom = cellposmom.append(Xmom)
    # temp dict of mean Δψ at +- range of dist from budneck
    tempdic = {dist: vf.neckDY(cell, cell.neckpos, dist)
               for dist in [.15, .3, .5]}
    temp = pd.DataFrame({'bud': pd.Series({k: tempdic[k][0] for k in tempdic}),
                         'mom': pd.Series({k: tempdic[k][1] for k in tempdic}),
                         'cellname': key})
    temp['dist'] = temp.index
    temp.set_index('cellname', inplace=True)
    neckregion = neckregion.append(temp, ignore_index=False)
