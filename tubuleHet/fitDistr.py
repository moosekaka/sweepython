# -*- coding: utf-8 -*-
"""
Created on Thu Jul 02 16:24:22 2015
Fit uniform, normal or shuffled distribution of cell DYs
@author: sweel
"""
from tvtk.api import tvtk
from collections import defaultdict
from vtkfuncs import vtkdata, vtklineids, vtkshuf, vtksamp
# pylint: disable=C0103


def fitDist(files, Graphs):
    """Return fitted distributions for bootstrapping
     branchpoints Delta Psi

    Parameters
    ----------
    files :
        List of Norm vtk files
    Graphs :
        Network x graph objects, generated from
        03createEdgeNodeListMulti.py

    Returns
    -------
    data :
        Dictionary of vtk object indexed by cell name
    sampN :
        Dictionary of normal distribution based on mean and std
        of cell indexed
    sampN :
        Dictionary of uniform distribution based mean-2*std
        (clipped at zero) and mean+4*std of cell indexed
    Norm :
        Dictionary of scaled DY values listed by cell edges of
        cell indexed
    lineId :
        list of lines with the start and end point vtk Index and also
        type of point (branch/end)
    """
    data = {}
    sampN = {}
    sampU = {}
    sampNRaw = {}
    sampURaw = {}
    Norm = defaultdict(dict)
    NormPermute = defaultdict(dict)
    NormRaw = defaultdict(dict)
    NormPerRaw = defaultdict(dict)
    lineId = defaultdict(dict)

    for el in files:
        filekey = el.rsplit('\\', 1)[1][5:][:-13]
        reader = tvtk.PolyDataReader()
        reader.set(file_name=el)
        reader.update()
        data[filekey] = reader.output
        print filekey
#       actual distribution
        Norm[filekey] = vtkdata(data[filekey])
        NormRaw[filekey] = vtkdata(data[filekey], voi='DY_raw')
        lineId[filekey] = vtklineids(data[filekey], Graphs[filekey])
#       shuffle distribution
        NormPermute[filekey] = vtkshuf(data[filekey])
        NormPerRaw[filekey] = vtkshuf(data[filekey], voi='DY_raw')
#       random distributions
        sampN[filekey] = vtksamp(data[filekey])[0]
        sampU[filekey] = vtksamp(data[filekey])[1]
        sampNRaw[filekey] = vtksamp(data[filekey], voi='DY_raw')[0]
        sampURaw[filekey] = vtksamp(data[filekey], voi='DY_raw')[1]

    return(data, sampN, sampU, Norm, NormPermute,
           sampNRaw, sampURaw, NormRaw, NormPerRaw,
           lineId)


