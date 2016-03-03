"""
Module to create mitonetwork in graph representation, with nodes, edges and
their respective attributes
"""
import os
import os.path as op
import math
import cPickle as pickle
from collections import defaultdict
import numpy as np
import networkx as nx
import vtk
from tvtk.api import tvtk
import wrappers as wr
# pylint: disable=C0103
vtkF = defaultdict(dict)
datadir = op.join(os.getcwd(), 'data')
rawdir = op.join(os.getcwd(), 'output')


def disteuc(pt1, pt2):
    """
    returns euclidian dist btw pt1 and pt2
    """
    test = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(pt1, pt2))
    return test


def findedges(x, y, grph):
    """
    compares points x, y with every other node in grph to find connected
    edges
    """
    for n, nat in grph.nodes(data=True):
        d = disteuc(nat['coord'], x)  # faster than np.linalg.norm
        if d == 0:
            n1 = n
        d = disteuc(nat['coord'], y)
        if d == 0:
            n2 = n
    return n1, n2


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
    first = {}
    last = {}
    scalars = tvtk.to_tvtk(
        vtkdata.GetPointData().GetScalars(scalartype)).to_array()
    points = tvtk.to_tvtk(vtkdata.GetPoints()).to_array()
    G = nx.MultiGraph(cell=graphname)
    for i in range(vtkdata.GetNumberOfCells()):
        temp = vtkdata.GetCell(i).GetPointIds()
        mx = temp.GetNumberOfIds()
        first[i] = temp.GetId(0)   # first point ids dictionary
        last[i] = temp.GetId(mx-1)  # last point ids dictionary

    fp = points[first.values()]  # first point coordinates
    lp = points[last.values()]  # last point coordinates
    exist = {}

    for i in range(vtkdata.GetNumberOfLines()-1, -1, -1):
        pt_id_1 = first[i]
        pt_id_2 = last[i]
        inten1 = scalars[pt_id_1]
        inten2 = scalars[pt_id_2]
        single_pt = fp[i] - lp[i]
        a = fp[i, np.newaxis]-lp[:i]
        b = fp[i, np.newaxis]-fp[:i]
        c = lp[i, np.newaxis]-fp[:i]
        d = lp[i, np.newaxis]-lp[:i]
        dics = {'a': a, 'b': b, 'c': c, 'd': d}
        exist['single_pt'] = np.all(single_pt == 0)
        for key in sorted(dics):
            exist[key] = np.any(np.all(dics[key] == 0, axis=1))
        if not (exist['a'] or exist['b'] or exist['single_pt']):
            G.add_node(nnodes, coord=tuple(fp[i]), inten=inten1)
            nnodes += 1
        if not (exist['c'] or exist['d']):
            G.add_node(nnodes, coord=tuple(lp[i]), inten=inten2)
            nnodes += 1

    for i in range(vtkdata.GetNumberOfLines()):
        pt_id_1 = first[i]
        pt_id_2 = last[i]
        r1 = points[pt_id_1]
        r2 = points[pt_id_2]
        mdpt = (
            (r1[0]+r2[0]) / 2, (r1[1]+r2[1]) / 2, (r1[2]+r2[2]) / 2)

#  Calc edge weight by Euc. distance between adjacent pixels
        pids = tvtk.to_tvtk(vtkdata.GetCell(i).GetPointIds())
        shift_pids = np.roll(pids, 1)
        edw = np.sum(np.linalg.norm(points[shift_pids][1:] -
                                    points[pids][1:], axis=1))

        n1, n2 = findedges(r1, r2, G)
        G.add_edge(n1, n2, weight=edw, cellID=i, midpoint=mdpt)

    deg = nx.degree(G)
    for i in G.nodes_iter():
        G.node[i]['degree'] = deg[i]
    nds = G.nodes(data=True)
    edgs = G.edges(data=True)
    return nds, edgs, G

# ===========================================================================
if __name__ == '__main__':
    vtkF = wr.ddwalk(op.join(rawdir, 'normalizedVTK'),
                     '*skeleton.vtk', start=5, stop=-13)

    for mediatype in sorted(vtkF.keys())[:]:
        nlist = []
        elist = []
        glist = []
        print 'creating edge node lists for %s' % mediatype
        print'number of files = %-3d\n' % len(vtkF[mediatype])
        for files in sorted(vtkF[mediatype].keys())[:]:
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
