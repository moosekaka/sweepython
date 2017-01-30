# -*- coding: utf-8 -*-
"""
Module for functions for mito network normalization
@author: sweel_rafelski
"""
import os
import os.path as op
import vtk
import vtk.util.numpy_support as vnpy
import numpy as np
from tvtk.api import tvtk
# pylint: disable=C0103
datadir = op.join(os.getcwd())


def add_scalars(fskel, fvr, fvg, radius=2.5):
    """
    Add scalar values from voxels data (eg. resampledVTK)

    Parameters
    ----------
    fskel : VTK file
        VTK file from MitoGraph
    fvr, fvg : vtk resampled file
        Obtained from MitoGraph by specifying -export_resampled_image option
    radius : float
        radius of point cloud sphere to average with

    Returns
    -------
    polydata : VTK poly
        polydata object with raw scalar values and width
    """
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(fskel)
    reader.Update()
    reader2 = vtk.vtkStructuredPointsReader()
    reader2.SetFileName(fvg)
    reader2.Update()
    reader3 = vtk.vtkStructuredPointsReader()
    reader3.SetFileName(fvr)
    reader3.Update()

    dataSkel = reader.GetOutput()  # Skel coords value from skel file
    dataGFPV = reader2.GetOutput()  # voxels GFP
    dataRFPV = reader3.GetOutput()
    ptsOld = dataSkel.GetPoints()
    Cell = vtk.vtkCellArray()

    intenGFP = dataGFPV.GetPointData().GetScalars().GetTuple1
    intenRFP = dataRFPV.GetPointData().GetScalars().GetTuple1
    skelWidth = dataSkel.GetPointData().GetScalars('Width')
    loc = vtk.vtkPointLocator()
    loc.SetDataSet(dataGFPV)
    loc.BuildLocator()
    result = vtk.vtkIdList()
    rawRFP = vtk.vtkDoubleArray()
    rawRFP.SetName("rRFP")
    rawGFP = vtk.vtkDoubleArray()
    rawGFP.SetName("rGFP")

#   add the lines/ cells  for connectivity of skel info
    for i in range(dataSkel.GetNumberOfLines()):
        ptID = dataSkel.GetCell(i).GetPointIds()
        oldId = [
            ptID.GetId(pid) for pid in range(ptID.GetNumberOfIds())]
        Cell.InsertNextCell(len(oldId))
        for j in oldId:
            Cell.InsertCellPoint(j)

#   averaging of pts intensity value surrounding each point in skel
    for n in range(dataSkel.GetNumberOfPoints()):

        ptOI = tuple(np.ceil(i/.055) for i in dataSkel.GetPoint(n))

        loc.FindPointsWithinRadius(radius, ptOI, result)

        voxID = [
            result.GetId(i) for i in range(result.GetNumberOfIds())]

        g = np.mean([intenGFP(m) for m in voxID])
        rawGFP.InsertNextValue(g)

        r = np.mean([intenRFP(m) for m in voxID])
        rawRFP.InsertNextValue(r)
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(ptsOld)
    polydata.SetLines(Cell)
    polydata.GetPointData().AddArray(rawRFP)
    polydata.GetPointData().AddArray(rawGFP)
    polydata.GetPointData().AddArray(skelWidth)
    polydata.GetPointData().GetArray(2).SetName("TubeWidth")
    return polydata


def normSkel(polydata, backgrnd):
    """
    Normalize channels to correct for focal plane intensity variations

    Parameters
    ----------
    polydata : vtkPolyData
        Must be a VTK PolyData object (not tVTK)
    backgrnd : tuple
        tuple containing background values for each channel
    """
    polydata = tvtk.to_tvtk(polydata)  # convert to tvtk to have traits
    temp = polydata.point_data
    rawGFP = np.ravel(temp.get_array('rGFP'))
    rawRFP = np.ravel(temp.get_array('rRFP'))
    if backgrnd[0] > min(rawRFP):  # background has higher min than skel rfp
        min_rfp = min(rawRFP)-1.  # minus one to ensure no zero divisions
    else:
        min_rfp = backgrnd[0]-1.  # normally, min should be background RFP val

    min_gfp = min(backgrnd[1], min(rawGFP))
    lines = []

#   method for getting pointIds as unique list
    for el in range(polydata.number_of_lines):
        temp = np.ravel(polydata.get_cell(el).point_ids)
        lines.append(temp)
    pointIds = np.unique([el for ln in lines for el in ln])

#   background Substracted rfp and gfps
    rfp_bk = rawRFP-min_rfp
    gfp_bk = rawGFP-min_gfp
#   width equivalent
    W = rfp_bk / np.min(rfp_bk)
    DY = gfp_bk / W  # raw DY/W normalized values
#   rescale DY to minmax
    minDY = min([DY[i] for i in pointIds])
    maxDY = max([DY[i] for i in pointIds])
    normDY = ((DY-minDY)/(maxDY-minDY))
    return normDY, DY, rfp_bk, gfp_bk, W


def writevtk(data, fname, **kwargs):
    """
    Write out a vtk file using VTK polydata object *dat* and a filename *fname*
    with optional labels dictionary *kwargs* for the outputs

    kwargs
    ------
    Default dictionary keys are:

    * normalized_dy
    * unscaled_dy
    * ch1_bckgrnd
    * ch2_bckgrnd
    * width_eqv'
    """

    for k in sorted(kwargs):
        temp = vnpy.numpy_to_vtk(kwargs[k])
        temp.SetName(k)
        data.GetPointData().AddArray(temp)

    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(fname)
    writer.SetInputData(data)
    writer.Update()
