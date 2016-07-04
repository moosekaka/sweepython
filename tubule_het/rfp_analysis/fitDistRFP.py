# -*- coding: utf-8 -*-
"""
Created on Thu Jul 02 16:24:22 2015
Fit uniform, normal or shuffled distribution of cell DYs
@author: sweel
"""
from tvtk.api import tvtk
from collections import defaultdict
from tubule_het.autoCor.fitDistr import vtkdata, vtkshuf, vtksamp
# pylint: disable=C0103


def fitdrfp(files):
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
    sampNRaw = {}
    sampURaw = {}
    NormRaw = defaultdict(dict)
    NormPerRaw = defaultdict(dict)

    for el in files:
        filekey = el.rsplit('\\', 1)[1][5:][:-13]
        reader = tvtk.PolyDataReader()
        reader.set(file_name=el)
        reader.update()
        data[filekey] = reader.output
        print filekey
#       actual distribution
        NormRaw[filekey] = vtkdata(data[filekey], voi='tubeWidth')
#       shuffle distribution
        NormPerRaw[filekey] = vtkshuf(data[filekey], voi='tubeWidth')
#       random distributions
        sampNRaw[filekey] = vtksamp(data[filekey], voi='tubeWidth')[0]
        sampURaw[filekey] = vtksamp(data[filekey], voi='tubeWidth')[1]

    return(data, sampNRaw, sampURaw, NormRaw, NormPerRaw)
