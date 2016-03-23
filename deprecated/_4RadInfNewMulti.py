"""
    run createEdgeNodeList.py first to get the pickle file, output is for
    *celldata.pkl and *coordinateForMitoskel.pkl
"""
import vtk
import glob
import os
import scipy.stats as sp
import numpy as np
import cPickle as pickle
# pylint: disable=C0103
# =============================================================================
#       input file data and get pickle data input
# =============================================================================
fileLoc = os.getcwd()+r'\Norm*vtk'
files = glob.glob(fileLoc)
media = files[0].rsplit('\\', 2)
C = len(files)
INPUT = open(media[1]+'_grph.pkl', 'rb')  # from 03createEdgeNodeListMulti.py
idx = [i for i, j in enumerate(files)]  # filtered of files in ignorelist
nodesEdges = pickle.load(INPUT)
nodes = [nodesEdges[0][i] for i in idx]
edges = [nodesEdges[0][i] for i in idx]
G = [nodesEdges[2][i] for i in idx]

# =============================================================================
#           Main function block
# =============================================================================


def radInfPipe(radius):
    arrBrc = []
    cellAPmeans = []
    BrcMeans = []
    aveIntAllPoints = []
    brch = []
    ends = []
    aveIntNodes = {}
    NodeBrchMeanVar = {}
# pylint: enable=C0103
    RAD_INF = radius  # to give a sphere of influence equivalent to 0.4 microns

    for a in range(len(files)):
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(files[a])
        reader.Update()
        data = reader.GetOutput()
        scalars = data.GetPointData().GetScalars("DY_minmax")
        Max = scalars.GetRange()[1]
        Min = scalars.GetRange()[0]
        aveIntAllPoints.append(
            [scalars.GetTuple1(i) for i in range(scalars.GetNumberOfTuples())])
#   get  branch and endpoints coordinates
        curGrph = G[a]
        epoints = [
            j for i, j in enumerate(curGrph) if curGrph.node[j]['degree'] == 1]
        bpoints = [
            j for i, j in enumerate(curGrph) if curGrph.node[j]['degree'] > 1]
        brch.append([curGrph.node[j]['coord'] for j in bpoints])
        ends.append([curGrph.node[j]['coord'] for j in epoints])

# =============================================================================
#   add points of interest to the input array for FindPointsWithinNodeRadius
#   locator object , and use the adjacency iterator to get the adjacent edges
#   for a node to average the scalar value
# =============================================================================
#   arrays to hold averaged  and std node intensities in current cell
        nodesAveVal = []
        cloudStdv = []
        NodeBrchMeanVarTmp = []

#   adjacent  and nodes of interest:
        adnoi = [nbrs for n, nbrs in curGrph.adjacency_iter()]
        noi = [n for n, nbrs in curGrph.adjacency_iter()]

#  adjacent edges to the node OI, note different keys for each parallel edge!)
        for h, i in enumerate(adnoi):
            lineOI = [k['cellID'] for j in i.values() for k in j.values()]

#   reset each time we are assigning points from adnoi
            edgePts = vtk.vtkPoints()  # local points in the RAD_INF
            polyData = vtk.vtkPolyData()  # the object as input for Locator
            locInt = vtk.vtkDoubleArray()  # local scalars in the RAD_INF

#   arrays for calculating variance in branchpoints
            tempNBV = []
            temparrScal = [[] for i in range(len(lineOI))]
            temparrPts = [[] for i in range(len(lineOI))]

            for x, j in enumerate(lineOI):
                cellpts = data.GetCell(j).GetPoints()
                cellIds = data.GetCell(j).GetPointIds()

                for k in range(data.GetCell(j).GetNumberOfPoints()):
                    edgePts.InsertNextPoint(cellpts.GetPoint(k))
                    locInt.InsertNextValue(
                        scalars.GetTuple1(cellIds.GetId(k)))

#   arrays for calculating variance in branchpoints
                    temparrScal[x].append(scalars.GetTuple1(cellIds.GetId(k)))
                    temparrPts[x].append(cellpts.GetPoint(k))

#   locator VTK object to find points with radius of RAD_INF
            polyData.SetPoints(edgePts)
            polyData.GetPointData().SetScalars(locInt)
            result = vtk.vtkIdList()
            loc = vtk.vtkPointLocator()
            loc.SetDataSet(polyData)
            loc.BuildLocator()
            loc.FindPointsWithinRadius(
                RAD_INF, curGrph.node[noi[h]]['coord'], result)

#   The line adjacent to the node of interst and their pt IDS
            ptIDOI = [
                result.GetId(el) for el in range(result.GetNumberOfIds())]
            ptsValsRad = [polyData.GetPoint(el) for el in ptIDOI]

#   mean of points in each individual cell/line that is within the rad_inf
            nodesAveVal.append(
                np.mean([polyData.GetPointData().GetScalars().GetTuple1(el)
                        for el in ptIDOI]))

#   std deviation for all points within the region of interest
            cloudStdv.append(
                np.std([polyData.GetPointData().GetScalars().GetTuple1(el)
                       for el in ptIDOI]))

#   mean of points in each individual LINE that is within the radius of inf
            for x in range(len(lineOI)):
                tempNBV.append(
                    [temparrScal[x][i] for i, j
                     in enumerate(temparrPts[x])
                     if j in ptsValsRad])

                tempNBV[x] = [(y-Min)/(Max-Min) for y in tempNBV[x]]

#   the stdv of EACH line in the within the radius of inf
            NodeBrchMeanVarTmp.append(
                np.std([np.mean(i) for i in tempNBV]))

#   normalized inten values and degree of the nodes for each cell 'a'
        aveIntNodes[a] = {}
        NodeBrchMeanVar[a] = {}
        for w, x in enumerate(nodesAveVal):
            aveIntNodes[a][noi[w]] = (
                (x-Min)/(Max-Min), curGrph.node[noi[w]]['degree'])

#   array of stdv of each line within a node of interest
            NodeBrchMeanVar[a][noi[w]] = (
                NodeBrchMeanVarTmp[w], curGrph.node[noi[w]]['degree'])

# =============================================================================
#        Finally, flatten the branch/endpoints and print
# =============================================================================
    meanBP = np.mean(
        [aveIntNodes[a][el][0] for a
         in range(len(aveIntNodes)) for el
         in aveIntNodes[a]
         if aveIntNodes[a][el][1] > 1])

    meanEP = np.mean(
        [aveIntNodes[a][el][0] for a
         in range(len(aveIntNodes)) for el
         in aveIntNodes[a]
         if aveIntNodes[a][el][1] == 1])

# =============================================================================
#       Stat TESTs compare mean of branchpoints and mean of cell
# =============================================================================
#   pop inten mean for cell
    cellAPmeans = [np.mean(i) for i in aveIntAllPoints]

#   array of Branchpts int in that cell
    for a in range(len(G)):
        arrBrc.append(
            [aveIntNodes[a][el][0] for el
             in aveIntNodes[a]
             if aveIntNodes[a][el][1] > 1])

#   branchpoints mean for that cell
    BrcMeans = [np.mean(i) for i in arrBrc]

#   Test of significance for brcpts > cell mean
    TEST = [
        sp.ranksums(arrBrc[i], aveIntAllPoints[i]) for i
        in range(len(files))]

    EnrichBP = [
        i for i, j in enumerate(TEST)
        if j[1] < 0.05 and BrcMeans[i]-cellAPmeans[i] > 0]

#   Test of significance for brcpts != cell mean
    TEST2 = [
        sp.ranksums(arrBrc[i], aveIntAllPoints[i]) for i
        in range(len(files))]

    DiffBP = [i for i, j in enumerate(TEST2)
              if j[1] < 0.05]

    print("Mean of branchpoints =%7.4f\n"
          "Mean of endpoints =%7.4f\n"
          "Mean of allpoints =%7.4f\n"
          "Radius of influence =%4.1f\n"
          "Number of Enriched Branchpts cells =%2s\n"
          "Number of Different Branchpts cells =%2s\n" %
          (meanBP, meanEP,
           np.mean([el for l in aveIntAllPoints for el in l]),
           radius, len(EnrichBP), len(DiffBP)))

    return(EnrichBP, DiffBP, arrBrc, BrcMeans, cellAPmeans,
           NodeBrchMeanVar, brch, ends, curGrph, aveIntAllPoints,
           aveIntNodes, edgePts, locInt)

# =============================================================================
#           Run as Main
# =============================================================================
if __name__ == '__main__':

    (EnrichBP, DiffBP, arrBrc, BrcMeans, cellAPmeans,
     NodeBrchMeanVar, brch, ends, curGrph, aveIntAllPoints,
     aveIntNodes, edgePts, locInt) = radInfPipe(0.3)

    IO = tuple((EnrichBP, DiffBP, arrBrc, BrcMeans, cellAPmeans))
    output = open(media[1]+'_cellData.pkl', 'wb')
    pickle.dump(IO, output)
    output.close()

    IO2 = tuple((brch, ends, NodeBrchMeanVar))
    output2 = open(media[1]+'coordinateForMitoskel.pkl', 'wb')
    pickle.dump(IO2, output2)
    output2.close()
