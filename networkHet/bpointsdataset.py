# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 16:15:25 2015

@author: sweel
"""
"""
    Calculate the statistical and topological measures of interest of cells
    in various carbon sources and munge it into dataframe
"""
import os
import fnmatch
import numpy as np
from tvtk.api import tvtk
import networkx as nx
import scipy.stats as sp
import matplotlib.pyplot as plt
from collections import defaultdict
from mungedata import MungeDataFuncs as md
from numpy.random import choice as samp_no_rep
import pandas as pd
import cPickle as pickle
import seaborn as sns
plt.close('all')


def bpts_inten(vtkdata, bptscoord, radinf=.3):
    """return list of branchpoints ptIDs in nodes with
    averaged intensity values within a radius of influence from data
    """
    big = {}
    polydata = tvtk.PolyData()  # vtk object as input for Locator
    polydata.points = vtkdata.points

    for point in bptscoord:
        result = tvtk.IdList()
        loc = tvtk.PointLocator(data_set=polydata)
        loc.build_locator()
        loc.find_points_within_radius(
            radinf, bptscoord[point], result)
        big[point] = result
    return big
# =============================================================================
#               init vars and get vtk, graph data
# =============================================================================
# pylint: disable=C0103
vtkF = {}
G = {}
backgroundGFP = {}
backgroundRFP = {}
parDir = os.path.dirname(os.getcwd())
for root, dirs, files in os.walk(os.getcwd()):
    for i in files:
        if fnmatch.fnmatch(i, '*vtk'):
            vtkF.setdefault(root.rsplit('\\', 1)[1], []).append(
                os.path.join(root, i))
        if fnmatch.fnmatch(i, '*grph.pkl'):
            G.setdefault(root.rsplit('\\', 1)[1], []).append(
                os.path.join(root, i))

media = sorted(vtkF.keys())

#    networkX graph objects of mitograph
for i in G:
    with open(G[i][0], 'rb') as inpt:
        temp = pickle.load(inpt)[2]
        G[i].append(temp)

#    get metadatas
with open(parDir+'\\'+'fileMetas.pkl', 'rb') as inpt:
    METAS = pickle.load(inpt)
for i in METAS:
    backgroundRFP[i] = METAS[i][0]
    backgroundGFP[i] = METAS[i][1]


# =============================================================================
#               begin pipeline
# =============================================================================
# pylint: disable=C0103
mito_avgdeg = defaultdict(dict)
mito_bpts_dy = defaultdict(dict)
mito_bpts_dyraw = defaultdict(dict)

for mem in media:
    print'\nNow on %s\n' % mem + "=" * 79
    for n, a in enumerate(vtkF[mem]):
        Norm = []
        NormRaw = []
        GFP = []
        RFP = []
        arrPts = []
        W = []
        W2 = []
        rGFP = []
        lineId = {}

        curGrph = G[mem][1][n]
        reader = tvtk.PolyDataReader()
        reader.set(file_name=a)
        reader.update()
        data = reader.output
        scalarsNorm = data.point_data.scalars
        temp = data.point_data
        dyRaw = np.ravel(temp.get_array('DY_raw'))
        rawGFP = np.ravel(temp.get_array('rGFP'))
        rawRFP = np.ravel(temp.get_array('rRFP'))
        WidthEq = np.ravel(temp.get_array('WidthEq'))
        tubeWidth = np.ravel(temp.get_array('tubeWidth'))
        filekey = a.rsplit('\\', 1)[1][5:][:-13]

        if backgroundRFP[filekey] > min(rawRFP):
            minA = backgroundRFP[filekey]-1
        else:
            minA = backgroundRFP[filekey]-1
        minB = min(backgroundGFP[filekey], min(rawGFP))

#    get btps intensity value within radofinfluence
        branchpoints = {j: attr['coord'] for j, attr
                        in curGrph.nodes(data=True)
                        if attr['degree'] > 2}

        bptpid = md.bpts_inten(data, branchpoints)
        bptdy = {key: np.mean([scalarsNorm[el] for el in vals])
                 for key, vals in sorted(bptpid.iteritems())}

        bptdy_raw = {key: np.mean([dyRaw[el] for el in vals])
                     for key, vals in sorted(bptpid.iteritems())}

        bpids = np.unique([el for lis in bptpid.values() for el in lis])
        allids = {i: np.array(data.get_cell(i).point_ids) for i
                  in range(data.number_of_lines)}
        nonbpids = np.unique([el for lis in allids.values() for el
                             in lis if el not in bpids])

        fakebp = {}
        for k, v in bptpid.items():
            fakebp[k] = samp_no_rep(nonbpids, size=len(v), replace=True)

        fakebpdy_raw = {key: np.mean([dyRaw[el] for el in vals])
                        for key, vals in sorted(fakebp.iteritems())}
