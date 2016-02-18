"""
     make sure VTK file exist, if not run MapZisoPtcldRAD.py first!
     this module returns edgeslist,nodelist and graph objects fromm VTK files
     outputs the graph objects as a pickle file *_grph.pkl
"""
import os
import os.path as op
import math
import cPickle as pickle
import fnmatch
from collections import defaultdict
import numpy as np
import networkx as nx
import vtk
# pylint: disable=C0103
vtkF = defaultdict(dict)
datadir = op.join(os.getcwd(), 'data')


def disteuc(pt1, pt2):
    """returns euclidian dist btw pt1 and pt2
    """
    return math.sqrt(vtk.vtkMath.Distance2BetweenPoints(pt1, pt2))


def makegraph(vtkdata, graphname, scalartype='DY_raw'):
    """
    Return networkX graph object from vtk skel

    Parameters
    ----------
    vtkdata: vtkPolyData
        must use vtk reader, not mayavi

    graphname : str
        name for graph

    scalartype : str
        point data type to plot on skeleton, can be of type:
        DY_minmax,
        WidthEq,
        DY_raw,
        rRFP,
        rGFP,
        bkstRFP,
        bkstGFP,

    Returns
    -------
    nds, edgs : list
        list of nodes and edges data
    G : networkX
        graph object
    """
    nnodes = 0
    scalars = vtkdata.GetPointData().GetScalars(scalartype)
    num_lines = vtkdata.GetNumberOfLines()
    G = nx.MultiGraph(cell=graphname)
# =============================================================================
#   for every line in each cell, starting from the highest index compare the
#   end points of that line with another line of the next lower index, if
#   pythagorean dist==0, label it a node and add it to the graph object G[a]
# =============================================================================
    for i in range(num_lines-1, -1, -1):
        line1 = vtkdata.GetCell(i).GetPoints()
        np_pointsi = line1.GetNumberOfPoints()
        pointID = vtkdata.GetCell(i).GetPointIds()
        exist1 = exist2 = False
        r1 = line1.GetPoint(0)
        r2 = line1.GetPoint(np_pointsi-1)
        inten1 = scalars.GetTuple1(pointID.GetId(0))
        inten2 = scalars.GetTuple1(pointID.GetId(np_pointsi-1))

        for j in range(i+1):
            line2 = vtkdata.GetCell(j).GetPoints()
            np_pointsj = line2.GetNumberOfPoints()
            s1 = line2.GetPoint(0)
            s2 = line2.GetPoint(np_pointsj-1)

#    compare same line, will identify a line of zero length if true
            if i == j:
                d = disteuc(r1, r2)
                if d == 0:
                    exist2 = True
            else:
                d = disteuc(r1, s1)
                if d == 0:
                    exist1 = True
                d = disteuc(r1, s2)
                if d == 0:
                    exist1 = True
                d = disteuc(r2, s1)
                if d == 0:
                    exist2 = True
                d = disteuc(r2, s2)
                if d == 0:
                    exist2 = True
        if exist1 is False:
            G.add_node(nnodes, coord=r1, inten=inten1)
            nnodes += 1
        if exist2 is False:
            G.add_node(nnodes, coord=r2, inten=inten2)
            nnodes += 1

# =============================================================================
#    for every node identified in the list of nnodes, compare each line in
#    that cell with that node to identify the edges
# =============================================================================

    for i in range(num_lines):
        EDL = 0
        line1 = vtkdata.GetCell(i).GetPoints()
        np_pointsi = line1.GetNumberOfPoints()
        r1 = line1.GetPoint(0)
        r2 = line1.GetPoint(np_pointsi-1)
        pointID = vtkdata.GetCell(i).GetPointIds()
        mdpt = (
            (r1[0]+r2[0]) / 2, (r1[1]+r2[1]) / 2, (r1[2]+r2[2]) / 2)

#   for each line sum thru the line the distance between points to
#   get edge len EDL
        EDL = np.sum(
            [disteuc(line1.GetPoint(pid), line1.GetPoint(pid-1))
             for pid in range(1, np_pointsi)])

        for j in G.nodes_iter():
            nod = G.node[j]
            r = nod['coord']
            d = disteuc(r, r1)
            if d == 0:
                n1 = j
            d = disteuc(r, r2)
            if d == 0:
                n2 = j
        G.add_edge(
            n1, n2, weight=EDL, cellID=i, midpoint=mdpt)

    deg = nx.degree(G)
    for i in G.nodes_iter():
        G.node[i]['degree'] = deg[i]
    nds = G.nodes(data=True)
    edgs = G.edges(data=True)
    return nds, edgs, G

# =============================================================================
#      main
# =============================================================================
# pylint: enable=C0103
if __name__ == '__main__':

    for root, dirs, files in os.walk(op.join(datadir, 'normSkel')):
        for i in files:
            if fnmatch.fnmatch(i, '*skeleton.vtk'):
                media = root.rsplit('\\', 1)[1]
                vtkF[media][i[5:-13]] = os.path.join(root, i)

    for mediatype in sorted(vtkF.keys()):
        nlist = []
        elist = []
        glist = []
        print 'creating edge node lists for %s' % mediatype
        print'number of files = %-3d\n' % len(vtkF[mediatype])
        for files in sorted(vtkF[mediatype]):
            reader = vtk.vtkPolyDataReader()
            reader.SetFileName(vtkF[mediatype][files])
            reader.Update()
            data = reader.GetOutput()
            node_data, edge_data, nxgrph = makegraph(data, files)

            nlist.append(node_data)
            elist.append(edge_data)
            glist.append(nxgrph)

            filename = op.join(datadir, '%s_grph.pkl' % mediatype)
            with open(filename, 'wb') as output:
                pickle.dump((nlist, elist, glist), output)
