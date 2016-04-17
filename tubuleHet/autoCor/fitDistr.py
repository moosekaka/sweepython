# -*- coding: utf-8 -*-
"""
Created on Thu Jul 02 16:24:22 2015
Fit uniform, normal or shuffled distribution of cell DYs
@author: sweel
"""
from collections import defaultdict
import numpy as np
from numpy.random import choice as samp_no_rep
import scipy.stats as sp
# pylint: disable=C0103


def getBptsEpts(vtkData, curGrph):
    """Return branchpoints and endpoints Indexes

    Parameters
    ----------
    vtkData :
        Vtk polydata cell
    curGrph :
        Network x graph
    """
    bpts = [
        curGrph.node[i]['coord'] for i in curGrph.nodes()
        if curGrph.node[i]['degree'] > 1]
    bptsId = [vtkData.find_point(el) for el in bpts]

    epts = [
        curGrph.node[i]['coord'] for i in curGrph.nodes()
        if curGrph.node[i]['degree'] == 1]
    eptsId = [vtkData.find_point(el) for el in epts]

    return(bptsId, eptsId)


def pointIdsList(data):
    """Return pointIdList

    Parameters
    ----------
    data :
        Vtk polydata cell
    """
    lines = []
    for el in range(data.number_of_lines):
        temp = np.ravel(data.get_cell(el).point_ids)
        lines.append(temp)
    pointIds = np.unique([el for ln in lines for el in ln])
    return pointIds


def vtkdata(data, voi='DY_minmax'):
    """Return flattened list of variable of interest (VOR) from vtkdata
    default vor is DY_minmax if kwarg not specified

    Parameters
    ----------
    data:
        vtk reader output
    voi:
        rRFP|rGFP|TubeWidth|DY_minmax|DY_raw|WidthEq|bkstGFP|bkstRFP

    Returns
    -------
    norm:
        list of edges data values
     """
    norm = []
    pdata = data.point_data
    datavals = np.ravel(pdata.get_array(voi))

    for j in range(data.number_of_lines):
        cids = data.get_cell(j).point_ids
        norm.append([datavals[i] for i in cids])
    return norm


def vtklineids(data, graph):
    """Return flattened list of variable of interest (VOR) from vtkdata
    default vor is DY_minmax if kwarg not specified

    Parameters
    ----------
    data:
        vtk reader output

    Returns
    -------
    lid:
        list of edges data point ids
     """
    lid = defaultdict(dict)
    firstId = 0
    bpts, _ = getBptsEpts(data, graph)
    for j in range(data.number_of_lines):
        lid[j] = []
        cid = data.get_cell(j).point_ids
        temp = [cid[0], cid[-1]]  # first and last point on cell
        if temp[0] in bpts:
            lid[j].append(('b', firstId))
        else:
            lid[j].append(('e', firstId))

        if temp[1] in bpts:
            lid[j].append(('b', firstId+len(cid)))
        else:
            lid[j].append(('e', firstId+len(cid)))

        firstId += len(cid)
    return lid


def vtkshuf(data, voi='DY_minmax'):
    """Return shuffled list of variable of interest (VOR) from vtkdata
    default vor is DY_minmax if kwarg not specified

    Parameters
    ----------
    data:
        vtk reader output
    voi:
        rRFP|rGFP|TubeWidth|DY_minmax|DY_raw|WidthEq|bkstGFP|bkstRFP
    Returns
    -------
    normpermute:
        shuffled distr using np.choice (samp_no_rep)

     """
    normpermute = []
    ptIds = pointIdsList(data)
    pdata = data.point_data
    datavals = np.ravel(pdata.get_array(voi))

    for j in range(data.number_of_lines):
        cids = data.get_cell(j).point_ids

        # sampl cell pointIds w.o replace
        shuffledIds = samp_no_rep(ptIds, len(cids), replace=False)
        normpermute.append(
            [datavals[int(k)] for k in shuffledIds])
    return normpermute


def vtksamp(data, voi='DY_minmax'):
    """Return a fitted uniform and normal list of variable of interest (VOR)
    from vtkdata, default vor is DY_minmax if kwarg not specified

    Parameters
    ----------
    data:
        vtk reader output
    voi:
        rRFP|rGFP|TubeWidth|DY_minmax|DY_raw|WidthEq|bkstGFP|bkstRFP

    Returns
    -------
    sampN:
        scipy norm distr with mean and std of actual data
    sampU:
        scipy unifrm distr bounded by .01 and .99 percentile of act dist
     """
    temp = data.point_data
    sampN = []
    sampU = []
    datavals = np.ravel(temp.get_array(voi))
    cellMeans = np.mean(datavals)
    cellStds = np.std(datavals)

    N = sp.norm(cellMeans, cellStds)  # norm distribution
#    q1 = np.percentile(datavals, .0)
#    q2 = np.percentile(datavals, 1)
#    sampU = sp.uniform(q1, q2)
    a = cellMeans - 1.5 * cellStds
    U = sp.uniform(  # uniform distribution
        a.clip(min=0), cellMeans + 1.5 * cellStds)
    for line in range(data.number_of_lines):
        M = data.get_cell(line).number_of_points
        sampN.append(N.rvs(size=M))
        sampU.append(U.rvs(size=M))
    return (sampN, sampU)


def fitDist(vdata, grph):
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

 #       actual distribution
    Norm = vtkdata(vdata)
    NormRaw = vtkdata(vdata, voi='DY_raw')
    lineId = vtklineids(vdata, grph)

#       shuffle distribution
    NormPermute = vtkshuf(vdata)
    NormPerRaw = vtkshuf(vdata, voi='DY_raw')

#       random distributions
    sampN = vtksamp(vdata)[0]
    sampU = vtksamp(vdata)[1]
    sampNRaw = vtksamp(vdata, voi='DY_raw')[0]
    sampURaw = vtksamp(vdata, voi='DY_raw')[1]

    return(sampN, sampU, Norm, NormPermute,
           sampNRaw, sampURaw, NormRaw, NormPerRaw,
           lineId)
