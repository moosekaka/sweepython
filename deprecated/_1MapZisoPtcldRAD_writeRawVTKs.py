"""
    First step in pipeline, creates 'raw' vtk files before normalizing
"""
import vtk
import glob
import string
import os
import numpy as np
# pylint: disable=C0103
radius = 2.5  # radius of averaging
b = os.getcwd()
fskel = glob.glob(b+'\\'+'*RFP*skel*')
fVR = glob.glob(b+'\\'+'*RFP*resample*')
fVG = glob.glob(b+'\\'+'*GFP*resample*')  # GFP voxels
media = b.rsplit('\\', 1)
# pylint: enable=C0103
for x in range(len(fskel)):
    print'now on cell:%3d' % x
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(fskel[x])
    reader.Update()
    reader2 = vtk.vtkStructuredPointsReader()
    reader2.SetFileName(fVG[x])
    reader2.Update()
    reader3 = vtk.vtkStructuredPointsReader()
    reader3.SetFileName(fVR[x])
    reader3.Update()

    dataSkel = reader.GetOutput()  # Skel coords value from skel file
    dataGFPV = reader2.GetOutput()  # voxels GFP
    dataRFPV = reader3.GetOutput()

    maxR = dataRFPV.GetPointData().GetScalars().GetRange()[1]
    minR = dataRFPV.GetPointData().GetScalars().GetRange()[0]
    maxG = dataGFPV.GetPointData().GetScalars().GetRange()[1]
    minG = dataGFPV.GetPointData().GetScalars().GetRange()[0]

    ptsOld = dataSkel.GetPoints()
    Cell = vtk.vtkCellArray()
    dim = dataRFPV.GetDimensions()
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

# =============================================================================
#   add the lines/ cells  for connectivity of skel info
# =============================================================================
    for i in range(dataSkel.GetNumberOfLines()):
        oldL = dataSkel.GetCell(i).GetPoints()
        ptID = dataSkel.GetCell(i).GetPointIds()

        oldPt = [
            oldL.GetPoint(idx) for idx in range(oldL.GetNumberOfPoints())]

        oldId = [
            ptID.GetId(pid) for pid in range(ptID.GetNumberOfIds())]

        Cell.InsertNextCell(len(oldId))

        for j in oldId:
            Cell.InsertCellPoint(j)

# =============================================================================
#   averaging of pts intensity value surrounding each point in skel
# =============================================================================
    for n in range(dataSkel.GetNumberOfPoints()):

        ptOI = tuple(np.ceil(i/.055) for i in dataSkel.GetPoint(n))

        loc.FindPointsWithinRadius(radius, ptOI, result)

        voxID = [
            result.GetId(i) for i in range(result.GetNumberOfIds())]

        g = np.mean([intenGFP(m) for m in voxID])
        rawGFP.InsertNextValue(g)

        r = np.mean([intenRFP(m) for m in voxID])
        rawRFP.InsertNextValue(r)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(ptsOld)
    polyData.SetLines(Cell)
    polyData.GetPointData().AddArray(rawRFP)
    polyData.GetPointData().AddArray(rawGFP)
    polyData.GetPointData().AddArray(skelWidth)

    list1 = []
    writer = vtk.vtkPolyDataWriter()

    fileName = fskel[x].rsplit('\\', 1)[1]
    fileString = os.path.join(
        os.getcwd(), string.join(
            (str(radius), 'raw', fileName), sep='_'))

    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(fileString)
    writer.SetInput(polyData)
    writer.Update()
