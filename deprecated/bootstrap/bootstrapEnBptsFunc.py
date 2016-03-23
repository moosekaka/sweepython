# -*- coding: utf-8 -*-
"""
Created on Thu Jul 02 13:35:06 2015
@author: sweel
"""
import vtk
import scipy.stats as sp
import numpy as np


def bootStrp(radius, files, Grphs):
    """Bootstrap branchpoints enriched in DY

        Parameters
        ----------
        radius :
            radius of influence
        files :
            list of Norm vtk files
        Graphs :
            Network x graph objects, generated from
            03createEdgeNodeListMulti.py

        Returns
        -------
            data, BrcMeans, cellAPmeans, aveIntAllPoints, aveIntNodes,
            nodesAveVal, enrichBP, diffBP
    """
    enrichBP = []
    diffBP = []
    cellAPmeans = []
    BrcMeans = []
    aveIntNodes = []
    scalars = [[] for el in range(len(files))]
    data = [[] for el in range(len(files))]
    RAD_INF = radius
    aveIntAllPoints = []

#   read in VTK data
    for el in range(len(files)):
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(files[el])
        reader.Update()
        data[el] = reader.GetOutput()
        scalars[el] = data[el].GetPointData().GetScalars("DY_minmax")
        aveIntAllPoints.append(
            [scalars[el].GetTuple1(i) for i
             in range(scalars[el].GetNumberOfTuples())])

#   define branch and endpoints
    for a in range(len(files)):
        nodesAveVal = []
        curGrph = Grphs[a]
        bpoints = [
            j for i, j in enumerate(curGrph)
            if curGrph.node[j]['degree'] > 1]

#   vtk locator object for file 'a'
        loc = vtk.vtkPointLocator()
        loc.SetDataSet(data[a])
        loc.BuildLocator()

#   iterate around the number of branch points
        for i in range(len(bpoints)):
            # randomly select an integer between 0 and maxpointID
            POI = np.random.random_integers(
                0, data[a].GetNumberOfPoints()-1)

            result = vtk.vtkIdList()
            loc.FindPointsWithinRadius(
                RAD_INF, data[a].GetPoint(POI), result)

#   point IDs and the coords
            ptIDOI = [
                result.GetId(el) for el in range(result.GetNumberOfIds())]

#   array of random points averaged around RAD_INF
            nodesAveVal.append(
                np.mean([data[a].GetPointData().GetScalars()
                        .GetTuple1(el) for el in ptIDOI]))

#   append to the giant cell the normalized brchpt inten values
        aveIntNodes.append([x for x in nodesAveVal])

    BrcMeans = [np.mean(i) for i in aveIntNodes]
    cellAPmeans = [np.mean(i) for i in aveIntAllPoints]

#   Test of significance for brcpts > cell mean
    TEST = [
        sp.ranksums(aveIntNodes[i], aveIntAllPoints[i]) for i
        in range(len(files))]

    enrichBP = [
        i for i, j in enumerate(TEST)
        if j[1] < 0.05 and BrcMeans[i]-cellAPmeans[i] > 0]

#   Test of significance for brcpts != cell mean
    TEST3 = [
        sp.ranksums(aveIntNodes[i], aveIntAllPoints[i]) for i
        in range(len(files))]

    diffBP = [i for i, j in enumerate(TEST3) if j[1] < 0.05]

    return(data, BrcMeans, cellAPmeans, aveIntAllPoints, aveIntNodes,
           nodesAveVal, enrichBP, diffBP)
