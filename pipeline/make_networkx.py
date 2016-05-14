"""
Module to create mitonetwork in graph representation, with nodes, edges and
their respective attributes
"""
import os
import os.path as op
import cPickle as pickle
import numpy as np
import networkx as nx
import vtk
from tvtk.api import tvtk
import wrappers as wr
# pylint: disable=C0103


def find_edges(X, Y, M):
    """
    find connections between line segment endpoints in list X and Y
    """
    n1 = X[:, np.newaxis]-M
    n2 = Y[:, np.newaxis]-M
    test1 = np.ravel([np.nonzero(np.all(n == 0, axis=1)) for n in n1])
    test2 = np.ravel([np.nonzero(np.all(n == 0, axis=1)) for n in n2])
    e_list = [tuple(sorted(i)) for i in zip(test1, test2)]
    return e_list


def make_ebunch(e_list, vtkdat, pnts, X, Y):
    """
    returns the edge attributes (edge length, etc.) from the edgelist e_list
    returned by find_edges
    """
    ebunch = {}
    for i, el in enumerate(e_list):
        p = vtkdat.GetCell(i).GetPointIds()
        pids = [p.GetId(k) for k in range(p.GetNumberOfIds())]

        # calc length of line by shifting line by one pixel and taking diff.
        # of every pixel of the pair of lines
        length = np.sum(np.linalg.norm(pnts[np.r_[0, pids[:-1]]][1:] -
                                       pnts[pids][1:],
                                       axis=1))
        ebunch[i] = (el[0], el[1],
                     {'weight': length,
                      'cellID': i,
                      'midpoint': tuple((X[i]+Y[i])/2)})
    return ebunch


def check_exist(bounds, l1, l2, ex_dic=None):
    """
    Check for existence of a connection between nodes.
    `a, b, c` and `d` are all the possible connections between the endpoints
    of lines `p1` and `p2`. Check each point for the connectiongs if they
    coincide (same coordinates).

    Returns:
    --------
    Dictionary of type (eg. 'a') as `True` if a connection exist for the type
    """
    if ex_dic is None:
        ex_dic = {}
    a = l1[bounds] - l2[:bounds]
    b = l1[bounds] - l1[:bounds]
    c = l2[bounds] - l1[:bounds]
    d = l2[bounds] - l2[:bounds]
    dics = {'a': a, 'b': b, 'c': c, 'd': d}
    for key in sorted(dics):
        ex_dic[key] = np.any(np.all(dics[key] == 0, axis=1))
    return ex_dic


def makegraph(vtkdata, graphname, scalartype='DY_raw'):
    """
    Return networkX graph object from vtk skel

    Parameters
    ----------
    vtkdata: vtkPolyData
        must use VTK, not tvtk (conversion is handled in function)

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
        `NetworkX` graph object
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

    # Create node list of graph
    for i in range(vtkdata.GetNumberOfLines()-1, -1, -1):
        exist = check_exist(i, fp, lp)
        # test for a single pixel
        exist['single_pt'] = np.all((fp[i] - lp[i]) == 0)
        inten1 = scalars[first[i]]
        inten2 = scalars[last[i]]
        # only add a node if there is no connection to another endpoint
        if not (exist['a'] or exist['b']):
            G.add_node(nnodes, coord=tuple(fp[i]), inten=inten1)
            nnodes += 1
        if not (exist['c'] or exist['d'] or exist['single_pt']):
            G.add_node(nnodes, coord=tuple(lp[i]), inten=inten2)
            nnodes += 1

    # Create edgelist of graph
    Ncoords = np.array([nat['coord'] for n, nat in G.nodes(data=True)])
    edges = find_edges(fp, lp, Ncoords)
    ebunch = make_ebunch(edges, vtkdata, points, fp, lp)
    G.add_edges_from(ebunch.values())

    # Degree connectivity of nodes
    for n in G.nodes():
        G.node[n]['degree'] = G.degree(n)
    return G.nodes(data=True), G.edges(data=True), G

# ===========================================================================
if __name__ == '__main__':
    # writes out a pickle file containing the graph list of every file for
    # for each mediatype
    datadir = op.join(os.getcwd(), 'mutants')
    rawdir = op.join(os.getcwd(), 'mutants')
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
