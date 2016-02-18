# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 00:27:44 2015
@author: sweel
"""
import numpy as np
from collections import defaultdict
from autoCor.graphfuncs import getBptsEpts
from numpy.random import choice as samp_no_rep
import scipy.stats as sp

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


def vtkdata(data, **kwargs):
    """Return flattened list of variable of interest (VOR) from vtkdata
    default vor is DY_minmax if kwarg not specified

    Parameters
    ----------
    data:
        vtk reader output
    kwargs:
        voi='WidthEq'| 'DY_raw'| 'rRFP'| 'rGFP'| 'bkstRFP'| 'bkstGFP'| 'tubeWidth'

    Returns
    -------
    norm:
        list of edges data values
     """
    norm = []
    voi = kwargs.pop("voi", None)
    if voi is None:
        datavals = data.point_data.scalars
    else:
        temp = data.point_data
        datavals = np.ravel(temp.get_array(voi))

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


def vtkshuf(data, **kwargs):
    """Return shuffled list of variable of interest (VOR) from vtkdata
    default vor is DY_minmax if kwarg not specified

    Parameters
    ----------
    data:
        vtk reader output
    kwargs:
        voi='WidthEq'| 'DY_raw'| 'rRFP'| 'rGFP'| 'bkstRFP'| 'bkstGFP'| 'tubeWidth'
    Returns
    -------
    normpermute:
        shuffled distr using np.choice (samp_no_rep)

     """
    normpermute = []
    voi = kwargs.pop("voi", None)
    ptIds = pointIdsList(data)
    if voi is None:
        datavals = data.point_data.scalars
    else:
        temp = data.point_data
        datavals = np.ravel(temp.get_array(voi))

    for j in range(data.number_of_lines):
        cids = data.get_cell(j).point_ids

        # sampl cell pointIds w.o replace
        shuffledIds = samp_no_rep(ptIds, len(cids), replace=False)
        normpermute.append(
            [datavals[int(k)] for k in shuffledIds])
    return normpermute


def vtksamp(data, **kwargs):
    """Return a fitted uniform and normal list of variable of interest (VOR)
    from vtkdata, default vor is DY_minmax if kwarg not specified

    Parameters
    ----------
    data:
        vtk reader output
    kwargs:
        voi='WidthEq'| 'DY_raw'| 'rRFP'| 'rGFP'| 'bkstRFP'| 'bkstGFP'| 'tubeWidth'

    Returns
    -------
    sampN:
        scipy norm distr with mean and std of actual data
    sampU:
        scipy unifrm distr bounded by .01 and .99 percentile of act dist
     """
    voi = kwargs.pop("voi", None)
    if voi is None:
        datavals = data.point_data.scalars
    else:
        temp = data.point_data
        datavals = np.ravel(temp.get_array(voi))

    cellMeans = np.mean(datavals)
    cellStds = np.std(datavals)

    sampN = sp.norm(cellMeans, cellStds)  # norm distribution
#    q1 = np.percentile(datavals, .0)
#    q2 = np.percentile(datavals, 1)
#    sampU = sp.uniform(q1, q2)
    a = cellMeans - 1.5 * cellStds
    sampU = sp.uniform(  # uniform distribution
        a.clip(min=0), cellMeans + 1.5 * cellStds)
    return (sampN, sampU)


def sclminmax(data):
    """return a scaled min max collection of vtk cellarray data
    """
    vtkscaled = defaultdict(dict)
    for key in data.keys():
        flat = [el for lis in data[key] for el in lis]
        fmax = max(flat)
        fmin = min(flat)
        vtkscaled[key] = []
        for line in data[key]:
            vtkscaled[key].append([(el - fmin) /
                                   (fmax - fmin) for el in line])
    return vtkscaled
