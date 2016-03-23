"""
     make sure VTK file exist, if not run MapZisoPtcldRAD.py first!
     this module returns edgeslist,nodelist and graph objects fromm VTK files
     outputs the graph objects as a pickle file *_grph.pkl
"""
import os
import os.path as op
import math
import cPickle as pickle
from collections import defaultdict
import numpy as np
import networkx as nx
import vtk
import wrappers as wr
# pylint: disable=C0103
vtkF = defaultdict(dict)
datadir = op.join(os.getcwd(), 'data')
rawdir = op.join(os.getcwd(), 'output')

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

#    vtkF = wr.ddwalk(op.join(rawdir, 'normalizedVTK'),
#                     '*skeleton.vtk', start=5, stop=-13)
    vtkF=   vtkF = defaultdict(dict,
            {'YPD': {'YPD_042515_003_RFPstack_001': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_003_RFPstack_001_skeleton.vtk',
              'YPD_042515_003_RFPstack_002': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_003_RFPstack_002_skeleton.vtk',
              'YPD_042515_003_RFPstack_003': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_003_RFPstack_003_skeleton.vtk',
              'YPD_042515_005_RFPstack_004': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_005_RFPstack_004_skeleton.vtk',
              'YPD_042515_007_RFPstack_006': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_007_RFPstack_006_skeleton.vtk',
              'YPD_042515_009_RFPstack_007': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_009_RFPstack_007_skeleton.vtk',
              'YPD_042515_009_RFPstack_008': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_009_RFPstack_008_skeleton.vtk',
              'YPD_042515_009_RFPstack_010': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_009_RFPstack_010_skeleton.vtk',
              'YPD_042515_011_RFPstack_011': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_011_RFPstack_011_skeleton.vtk',
              'YPD_042515_011_RFPstack_012': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_011_RFPstack_012_skeleton.vtk',
              'YPD_042515_011_RFPstack_013': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_011_RFPstack_013_skeleton.vtk',
              'YPD_042515_011_RFPstack_014': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_011_RFPstack_014_skeleton.vtk',
              'YPD_042515_013_RFPstack_015': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_013_RFPstack_015_skeleton.vtk',
              'YPD_042515_013_RFPstack_016': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_013_RFPstack_016_skeleton.vtk',
              'YPD_042515_015_RFPstack_017': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_015_RFPstack_017_skeleton.vtk',
              'YPD_042515_015_RFPstack_018': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_015_RFPstack_018_skeleton.vtk',
              'YPD_042515_015_RFPstack_019': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_015_RFPstack_019_skeleton.vtk',
              'YPD_042515_017_RFPstack_020': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_017_RFPstack_020_skeleton.vtk',
              'YPD_042515_019_RFPstack_023': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_019_RFPstack_023_skeleton.vtk',
              'YPD_042515_019_RFPstack_024': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_019_RFPstack_024_skeleton.vtk',
              'YPD_042515_019_RFPstack_025': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042515_019_RFPstack_025_skeleton.vtk',
              'YPD_042715_001_RFPstack_027': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_001_RFPstack_027_skeleton.vtk',
              'YPD_042715_001_RFPstack_028': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_001_RFPstack_028_skeleton.vtk',
              'YPD_042715_001_RFPstack_029': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_001_RFPstack_029_skeleton.vtk',
              'YPD_042715_003_RFPstack_031': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_003_RFPstack_031_skeleton.vtk',
              'YPD_042715_005_RFPstack_033': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_005_RFPstack_033_skeleton.vtk',
              'YPD_042715_005_RFPstack_034': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_005_RFPstack_034_skeleton.vtk',
              'YPD_042715_007_RFPstack_035': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_007_RFPstack_035_skeleton.vtk',
              'YPD_042715_007_RFPstack_036': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_007_RFPstack_036_skeleton.vtk',
              'YPD_042715_009_RFPstack_038': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_009_RFPstack_038_skeleton.vtk',
              'YPD_042715_009_RFPstack_039': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_009_RFPstack_039_skeleton.vtk',
              'YPD_042715_009_RFPstack_040': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_009_RFPstack_040_skeleton.vtk',
              'YPD_042715_009_RFPstack_041': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_009_RFPstack_041_skeleton.vtk',
              'YPD_042715_009_RFPstack_042': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_009_RFPstack_042_skeleton.vtk',
              'YPD_042715_011_RFPstack_045': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_011_RFPstack_045_skeleton.vtk',
              'YPD_042715_013_RFPstack_048': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_013_RFPstack_048_skeleton.vtk',
              'YPD_042715_013_RFPstack_050': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_013_RFPstack_050_skeleton.vtk',
              'YPD_042715_015_RFPstack_051': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_015_RFPstack_051_skeleton.vtk',
              'YPD_042715_015_RFPstack_052': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_015_RFPstack_052_skeleton.vtk',
              'YPD_042715_015_RFPstack_053': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_015_RFPstack_053_skeleton.vtk',
              'YPD_042715_015_RFPstack_054': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_015_RFPstack_054_skeleton.vtk',
              'YPD_042715_017_RFPstack_056': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_017_RFPstack_056_skeleton.vtk',
              'YPD_042715_017_RFPstack_057': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042715_017_RFPstack_057_skeleton.vtk',
              'YPD_042915_003_RFPstack_059': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_003_RFPstack_059_skeleton.vtk',
              'YPD_042915_005_RFPstack_060': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_005_RFPstack_060_skeleton.vtk',
              'YPD_042915_005_RFPstack_061': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_005_RFPstack_061_skeleton.vtk',
              'YPD_042915_007_RFPstack_062': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_007_RFPstack_062_skeleton.vtk',
              'YPD_042915_007_RFPstack_063': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_007_RFPstack_063_skeleton.vtk',
              'YPD_042915_007_RFPstack_064': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_007_RFPstack_064_skeleton.vtk',
              'YPD_042915_011_RFPstack_066': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_011_RFPstack_066_skeleton.vtk',
              'YPD_042915_013_RFPstack_067': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_013_RFPstack_067_skeleton.vtk',
              'YPD_042915_013_RFPstack_068': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_013_RFPstack_068_skeleton.vtk',
              'YPD_042915_017_RFPstack_071': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_017_RFPstack_071_skeleton.vtk',
              'YPD_042915_017_RFPstack_072': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_017_RFPstack_072_skeleton.vtk',
              'YPD_042915_017_RFPstack_073': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_017_RFPstack_073_skeleton.vtk',
              'YPD_042915_017_RFPstack_074': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_017_RFPstack_074_skeleton.vtk',
              'YPD_042915_017_RFPstack_075': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_017_RFPstack_075_skeleton.vtk',
              'YPD_042915_019_RFPstack_077': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_019_RFPstack_077_skeleton.vtk',
              'YPD_042915_021_RFPstack_078': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_021_RFPstack_078_skeleton.vtk',
              'YPD_042915_021_RFPstack_079': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_021_RFPstack_079_skeleton.vtk',
              'YPD_042915_021_RFPstack_080': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_021_RFPstack_080_skeleton.vtk',
              'YPD_042915_024_RFPstack_082': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_024_RFPstack_082_skeleton.vtk',
              'YPD_042915_024_RFPstack_083': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_042915_024_RFPstack_083_skeleton.vtk',
              'YPD_052315_001_RFPstack_001': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_001_RFPstack_001_skeleton.vtk',
              'YPD_052315_001_RFPstack_002': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_001_RFPstack_002_skeleton.vtk',
              'YPD_052315_003_RFPstack_003': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_003_RFPstack_003_skeleton.vtk',
              'YPD_052315_003_RFPstack_004': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_003_RFPstack_004_skeleton.vtk',
              'YPD_052315_003_RFPstack_005': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_003_RFPstack_005_skeleton.vtk',
              'YPD_052315_005_RFPstack_006': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_005_RFPstack_006_skeleton.vtk',
              'YPD_052315_005_RFPstack_007': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_005_RFPstack_007_skeleton.vtk',
              'YPD_052315_005_RFPstack_008': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_005_RFPstack_008_skeleton.vtk',
              'YPD_052315_007_RFPstack_010': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_007_RFPstack_010_skeleton.vtk',
              'YPD_052315_009_RFPstack_011': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_009_RFPstack_011_skeleton.vtk',
              'YPD_052315_009_RFPstack_012': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_009_RFPstack_012_skeleton.vtk',
              'YPD_052315_011_RFPstack_013': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_011_RFPstack_013_skeleton.vtk',
              'YPD_052315_011_RFPstack_014': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_011_RFPstack_014_skeleton.vtk',
              'YPD_052315_013_RFPstack_015': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_013_RFPstack_015_skeleton.vtk',
              'YPD_052315_013_RFPstack_016': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_013_RFPstack_016_skeleton.vtk',
              'YPD_052315_015_RFPstack_017': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_015_RFPstack_017_skeleton.vtk',
              'YPD_052315_015_RFPstack_018': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_015_RFPstack_018_skeleton.vtk',
              'YPD_052315_017_RFPstack_021': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_017_RFPstack_021_skeleton.vtk',
              'YPD_052315_019_RFPstack_022': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_019_RFPstack_022_skeleton.vtk',
              'YPD_052315_019_RFPstack_023': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_019_RFPstack_023_skeleton.vtk',
              'YPD_052315_019_RFPstack_024': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_019_RFPstack_024_skeleton.vtk',
              'YPD_052315_019_RFPstack_025': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_019_RFPstack_025_skeleton.vtk',
              'YPD_052315_021_RFPstack_026': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_021_RFPstack_026_skeleton.vtk',
              'YPD_052315_021_RFPstack_027': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_021_RFPstack_027_skeleton.vtk',
              'YPD_052315_023_RFPstack_029': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_023_RFPstack_029_skeleton.vtk',
              'YPD_052315_023_RFPstack_030': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_023_RFPstack_030_skeleton.vtk',
              'YPD_052315_025_RFPstack_031': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_025_RFPstack_031_skeleton.vtk',
              'YPD_052315_027_RFPstack_032': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_027_RFPstack_032_skeleton.vtk',
              'YPD_052315_031_RFPstack_034': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_031_RFPstack_034_skeleton.vtk',
              'YPD_052315_031_RFPstack_035': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_031_RFPstack_035_skeleton.vtk',
              'YPD_052315_031_RFPstack_036': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_031_RFPstack_036_skeleton.vtk',
              'YPD_052315_033_RFPstack_037': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_033_RFPstack_037_skeleton.vtk',
              'YPD_052315_033_RFPstack_038': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPD\\Norm_YPD_052315_033_RFPstack_038_skeleton.vtk'},
             'YPE': {'YPE_042515_001_RFPstack_000': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042515_001_RFPstack_000_skeleton.vtk',
              'YPE_042515_005_RFPstack_001': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042515_005_RFPstack_001_skeleton.vtk',
              'YPE_042515_007_RFPstack_004': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042515_007_RFPstack_004_skeleton.vtk',
              'YPE_042515_009_RFPstack_007': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042515_009_RFPstack_007_skeleton.vtk',
              'YPE_042515_011_RFPstack_009': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042515_011_RFPstack_009_skeleton.vtk',
              'YPE_042515_013_RFPstack_010': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042515_013_RFPstack_010_skeleton.vtk',
              'YPE_042515_013_RFPstack_011': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042515_013_RFPstack_011_skeleton.vtk',
              'YPE_042515_015_RFPstack_012': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042515_015_RFPstack_012_skeleton.vtk',
              'YPE_042515_017_RFPstack_013': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042515_017_RFPstack_013_skeleton.vtk',
              'YPE_042515_017_RFPstack_014': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042515_017_RFPstack_014_skeleton.vtk',
              'YPE_042515_017_RFPstack_015': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042515_017_RFPstack_015_skeleton.vtk',
              'YPE_042515_019_RFPstack_016': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042515_019_RFPstack_016_skeleton.vtk',
              'YPE_042515_019_RFPstack_017': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042515_019_RFPstack_017_skeleton.vtk',
              'YPE_042715_001_RFPstack_018': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_001_RFPstack_018_skeleton.vtk',
              'YPE_042715_001_RFPstack_019': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_001_RFPstack_019_skeleton.vtk',
              'YPE_042715_001_RFPstack_020': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_001_RFPstack_020_skeleton.vtk',
              'YPE_042715_001_RFPstack_022': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_001_RFPstack_022_skeleton.vtk',
              'YPE_042715_003_RFPstack_023': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_003_RFPstack_023_skeleton.vtk',
              'YPE_042715_003_RFPstack_024': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_003_RFPstack_024_skeleton.vtk',
              'YPE_042715_003_RFPstack_026': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_003_RFPstack_026_skeleton.vtk',
              'YPE_042715_005_RFPstack_029': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_005_RFPstack_029_skeleton.vtk',
              'YPE_042715_005_RFPstack_030': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_005_RFPstack_030_skeleton.vtk',
              'YPE_042715_005_RFPstack_032': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_005_RFPstack_032_skeleton.vtk',
              'YPE_042715_007_RFPstack_033': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_007_RFPstack_033_skeleton.vtk',
              'YPE_042715_007_RFPstack_034': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_007_RFPstack_034_skeleton.vtk',
              'YPE_042715_007_RFPstack_035': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_007_RFPstack_035_skeleton.vtk',
              'YPE_042715_007_RFPstack_036': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_007_RFPstack_036_skeleton.vtk',
              'YPE_042715_009_RFPstack_038': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_009_RFPstack_038_skeleton.vtk',
              'YPE_042715_009_RFPstack_039': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_009_RFPstack_039_skeleton.vtk',
              'YPE_042715_009_RFPstack_040': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_009_RFPstack_040_skeleton.vtk',
              'YPE_042715_009_RFPstack_041': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_009_RFPstack_041_skeleton.vtk',
              'YPE_042715_011_RFPstack_043': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_011_RFPstack_043_skeleton.vtk',
              'YPE_042715_011_RFPstack_044': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_011_RFPstack_044_skeleton.vtk',
              'YPE_042715_013_RFPstack_045': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_013_RFPstack_045_skeleton.vtk',
              'YPE_042715_013_RFPstack_046': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_013_RFPstack_046_skeleton.vtk',
              'YPE_042715_013_RFPstack_047': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_013_RFPstack_047_skeleton.vtk',
              'YPE_042715_013_RFPstack_048': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_013_RFPstack_048_skeleton.vtk',
              'YPE_042715_013_RFPstack_049': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_013_RFPstack_049_skeleton.vtk',
              'YPE_042715_015_RFPstack_050': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_015_RFPstack_050_skeleton.vtk',
              'YPE_042715_015_RFPstack_051': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_015_RFPstack_051_skeleton.vtk',
              'YPE_042715_020_RFPstack_055': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_020_RFPstack_055_skeleton.vtk',
              'YPE_042715_023_RFPstack_056': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_023_RFPstack_056_skeleton.vtk',
              'YPE_042715_023_RFPstack_057': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_023_RFPstack_057_skeleton.vtk',
              'YPE_042715_023_RFPstack_058': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_023_RFPstack_058_skeleton.vtk',
              'YPE_042715_023_RFPstack_059': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_023_RFPstack_059_skeleton.vtk',
              'YPE_042715_023_RFPstack_060': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_042715_023_RFPstack_060_skeleton.vtk',
              'YPE_052315_001_RFPstack_000': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_001_RFPstack_000_skeleton.vtk',
              'YPE_052315_001_RFPstack_001': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_001_RFPstack_001_skeleton.vtk',
              'YPE_052315_003_RFPstack_002': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_003_RFPstack_002_skeleton.vtk',
              'YPE_052315_003_RFPstack_003': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_003_RFPstack_003_skeleton.vtk',
              'YPE_052315_003_RFPstack_004': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_003_RFPstack_004_skeleton.vtk',
              'YPE_052315_003_RFPstack_005': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_003_RFPstack_005_skeleton.vtk',
              'YPE_052315_005_RFPstack_006': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_005_RFPstack_006_skeleton.vtk',
              'YPE_052315_005_RFPstack_007': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_005_RFPstack_007_skeleton.vtk',
              'YPE_052315_005_RFPstack_008': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_005_RFPstack_008_skeleton.vtk',
              'YPE_052315_007_RFPstack_009': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_007_RFPstack_009_skeleton.vtk',
              'YPE_052315_007_RFPstack_010': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_007_RFPstack_010_skeleton.vtk',
              'YPE_052315_009_RFPstack_011': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_009_RFPstack_011_skeleton.vtk',
              'YPE_052315_011_RFPstack_012': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_011_RFPstack_012_skeleton.vtk',
              'YPE_052315_011_RFPstack_013': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_011_RFPstack_013_skeleton.vtk',
              'YPE_052315_011_RFPstack_014': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_011_RFPstack_014_skeleton.vtk',
              'YPE_052315_011_RFPstack_015': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_011_RFPstack_015_skeleton.vtk',
              'YPE_052315_013_RFPstack_016': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_013_RFPstack_016_skeleton.vtk',
              'YPE_052315_013_RFPstack_017': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_013_RFPstack_017_skeleton.vtk',
              'YPE_052315_013_RFPstack_018': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_013_RFPstack_018_skeleton.vtk',
              'YPE_052315_013_RFPstack_019': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_013_RFPstack_019_skeleton.vtk',
              'YPE_052315_015_RFPstack_020': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_015_RFPstack_020_skeleton.vtk',
              'YPE_052315_015_RFPstack_021': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_015_RFPstack_021_skeleton.vtk',
              'YPE_052315_017_RFPstack_022': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_017_RFPstack_022_skeleton.vtk',
              'YPE_052315_019_RFPstack_023': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_019_RFPstack_023_skeleton.vtk',
              'YPE_052315_019_RFPstack_024': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_019_RFPstack_024_skeleton.vtk',
              'YPE_052315_019_RFPstack_025': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_019_RFPstack_025_skeleton.vtk',
              'YPE_052315_019_RFPstack_026': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_019_RFPstack_026_skeleton.vtk',
              'YPE_052315_019_RFPstack_027': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_019_RFPstack_027_skeleton.vtk',
              'YPE_052315_021_RFPstack_029': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_021_RFPstack_029_skeleton.vtk',
              'YPE_052315_023_RFPstack_030': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_023_RFPstack_030_skeleton.vtk',
              'YPE_052315_023_RFPstack_031': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_023_RFPstack_031_skeleton.vtk',
              'YPE_052315_025_RFPstack_033': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_025_RFPstack_033_skeleton.vtk',
              'YPE_052315_025_RFPstack_034': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_025_RFPstack_034_skeleton.vtk',
              'YPE_052315_027_RFPstack_035': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_027_RFPstack_035_skeleton.vtk',
              'YPE_052315_027_RFPstack_036': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_027_RFPstack_036_skeleton.vtk',
              'YPE_052315_030_RFPstack_037': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_030_RFPstack_037_skeleton.vtk',
              'YPE_052315_030_RFPstack_038': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_030_RFPstack_038_skeleton.vtk',
              'YPE_052315_032_RFPstack_039': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_032_RFPstack_039_skeleton.vtk',
              'YPE_052315_032_RFPstack_041': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_032_RFPstack_041_skeleton.vtk',
              'YPE_052315_032_RFPstack_042': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_032_RFPstack_042_skeleton.vtk',
              'YPE_052315_032_RFPstack_043': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_032_RFPstack_043_skeleton.vtk',
              'YPE_052315_034_RFPstack_044': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_034_RFPstack_044_skeleton.vtk',
              'YPE_052315_034_RFPstack_046': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_034_RFPstack_046_skeleton.vtk',
              'YPE_052315_034_RFPstack_047': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_052315_034_RFPstack_047_skeleton.vtk',
              'YPE_c42515_001_RFPstack_062': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_001_RFPstack_062_skeleton.vtk',
              'YPE_c42515_001_RFPstack_063': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_001_RFPstack_063_skeleton.vtk',
              'YPE_c42515_001_RFPstack_064': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_001_RFPstack_064_skeleton.vtk',
              'YPE_c42515_001_RFPstack_065': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_001_RFPstack_065_skeleton.vtk',
              'YPE_c42515_001_RFPstack_066': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_001_RFPstack_066_skeleton.vtk',
              'YPE_c42515_003_RFPstack_067': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_003_RFPstack_067_skeleton.vtk',
              'YPE_c42515_005_RFPstack_069': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_005_RFPstack_069_skeleton.vtk',
              'YPE_c42515_007_RFPstack_070': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_007_RFPstack_070_skeleton.vtk',
              'YPE_c42515_007_RFPstack_071': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_007_RFPstack_071_skeleton.vtk',
              'YPE_c42515_007_RFPstack_072': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_007_RFPstack_072_skeleton.vtk',
              'YPE_c42515_009_RFPstack_074': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_009_RFPstack_074_skeleton.vtk',
              'YPE_c42515_009_RFPstack_075': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_009_RFPstack_075_skeleton.vtk',
              'YPE_c42515_011_RFPstack_076': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_011_RFPstack_076_skeleton.vtk',
              'YPE_c42515_011_RFPstack_077': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_011_RFPstack_077_skeleton.vtk',
              'YPE_c42515_011_RFPstack_078': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_011_RFPstack_078_skeleton.vtk',
              'YPE_c42515_011_RFPstack_079': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_011_RFPstack_079_skeleton.vtk',
              'YPE_c42515_013_RFPstack_081': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_013_RFPstack_081_skeleton.vtk',
              'YPE_c42515_013_RFPstack_082': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_013_RFPstack_082_skeleton.vtk',
              'YPE_c42515_017_RFPstack_084': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_017_RFPstack_084_skeleton.vtk',
              'YPE_c42515_019_RFPstack_087': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_019_RFPstack_087_skeleton.vtk',
              'YPE_c42515_019_RFPstack_088': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPE\\Norm_YPE_c42515_019_RFPstack_088_skeleton.vtk'},
             'YPL': {'YPL_042515_001_RFPstack_002': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_001_RFPstack_002_skeleton.vtk',
              'YPL_042515_001_RFPstack_004': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_001_RFPstack_004_skeleton.vtk',
              'YPL_042515_001_RFPstack_005': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_001_RFPstack_005_skeleton.vtk',
              'YPL_042515_008_RFPstack_009': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_008_RFPstack_009_skeleton.vtk',
              'YPL_042515_008_RFPstack_011': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_008_RFPstack_011_skeleton.vtk',
              'YPL_042515_011_RFPstack_012': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_011_RFPstack_012_skeleton.vtk',
              'YPL_042515_011_RFPstack_013': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_011_RFPstack_013_skeleton.vtk',
              'YPL_042515_011_RFPstack_014': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_011_RFPstack_014_skeleton.vtk',
              'YPL_042515_015_RFPstack_016': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_015_RFPstack_016_skeleton.vtk',
              'YPL_042515_015_RFPstack_018': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_015_RFPstack_018_skeleton.vtk',
              'YPL_042515_017_RFPstack_020': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_017_RFPstack_020_skeleton.vtk',
              'YPL_042515_017_RFPstack_021': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_017_RFPstack_021_skeleton.vtk',
              'YPL_042515_017_RFPstack_022': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_017_RFPstack_022_skeleton.vtk',
              'YPL_042515_019_RFPstack_024': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_019_RFPstack_024_skeleton.vtk',
              'YPL_042515_019_RFPstack_025': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_019_RFPstack_025_skeleton.vtk',
              'YPL_042515_021_RFPstack_028': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_021_RFPstack_028_skeleton.vtk',
              'YPL_042515_021_RFPstack_029': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_021_RFPstack_029_skeleton.vtk',
              'YPL_042515_021_RFPstack_030': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_021_RFPstack_030_skeleton.vtk',
              'YPL_042515_021_RFPstack_031': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_021_RFPstack_031_skeleton.vtk',
              'YPL_042515_021_RFPstack_032': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_021_RFPstack_032_skeleton.vtk',
              'YPL_042515_023_RFPstack_033': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_023_RFPstack_033_skeleton.vtk',
              'YPL_042515_023_RFPstack_034': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_023_RFPstack_034_skeleton.vtk',
              'YPL_042515_023_RFPstack_035': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_023_RFPstack_035_skeleton.vtk',
              'YPL_042515_030_RFPstack_039': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_030_RFPstack_039_skeleton.vtk',
              'YPL_042515_030_RFPstack_040': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_030_RFPstack_040_skeleton.vtk',
              'YPL_042515_030_RFPstack_041': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_030_RFPstack_041_skeleton.vtk',
              'YPL_042515_030_RFPstack_042': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_030_RFPstack_042_skeleton.vtk',
              'YPL_042515_030_RFPstack_043': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_030_RFPstack_043_skeleton.vtk',
              'YPL_042515_032_RFPstack_044': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_032_RFPstack_044_skeleton.vtk',
              'YPL_042515_032_RFPstack_045': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042515_032_RFPstack_045_skeleton.vtk',
              'YPL_042715_004_RFPstack_048': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_004_RFPstack_048_skeleton.vtk',
              'YPL_042715_004_RFPstack_050': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_004_RFPstack_050_skeleton.vtk',
              'YPL_042715_004_RFPstack_051': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_004_RFPstack_051_skeleton.vtk',
              'YPL_042715_006_RFPstack_052': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_006_RFPstack_052_skeleton.vtk',
              'YPL_042715_006_RFPstack_053': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_006_RFPstack_053_skeleton.vtk',
              'YPL_042715_008_RFPstack_054': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_008_RFPstack_054_skeleton.vtk',
              'YPL_042715_008_RFPstack_055': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_008_RFPstack_055_skeleton.vtk',
              'YPL_042715_008_RFPstack_056': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_008_RFPstack_056_skeleton.vtk',
              'YPL_042715_010_RFPstack_057': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_010_RFPstack_057_skeleton.vtk',
              'YPL_042715_010_RFPstack_058': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_010_RFPstack_058_skeleton.vtk',
              'YPL_042715_012_olp14RFPstack_059': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_012_olp14RFPstack_059_skeleton.vtk',
              'YPL_042715_016_RFPstack_063': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_016_RFPstack_063_skeleton.vtk',
              'YPL_042715_016_RFPstack_064': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_016_RFPstack_064_skeleton.vtk',
              'YPL_042715_016_RFPstack_065': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_016_RFPstack_065_skeleton.vtk',
              'YPL_042715_018_RFPstack_066': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_018_RFPstack_066_skeleton.vtk',
              'YPL_042715_018_RFPstack_067': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_018_RFPstack_067_skeleton.vtk',
              'YPL_042715_020_RFPstack_068': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_020_RFPstack_068_skeleton.vtk',
              'YPL_042715_022_RFPstack_069': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_022_RFPstack_069_skeleton.vtk',
              'YPL_042715_022_RFPstack_070': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_022_RFPstack_070_skeleton.vtk',
              'YPL_042715_024_RFPstack_071': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_024_RFPstack_071_skeleton.vtk',
              'YPL_042715_024_RFPstack_072': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_024_RFPstack_072_skeleton.vtk',
              'YPL_042715_024_RFPstack_073': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_024_RFPstack_073_skeleton.vtk',
              'YPL_042715_026_RFPstack_074': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_026_RFPstack_074_skeleton.vtk',
              'YPL_042715_026_RFPstack_075': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_026_RFPstack_075_skeleton.vtk',
              'YPL_042715_028_RFPstack_077': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_028_RFPstack_077_skeleton.vtk',
              'YPL_042715_028_RFPstack_078': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_028_RFPstack_078_skeleton.vtk',
              'YPL_042715_028_RFPstack_079': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_028_RFPstack_079_skeleton.vtk',
              'YPL_042715_030_RFPstack_080': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_030_RFPstack_080_skeleton.vtk',
              'YPL_042715_030_RFPstack_081': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042715_030_RFPstack_081_skeleton.vtk',
              'YPL_042915_001_RFPstack_082': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_001_RFPstack_082_skeleton.vtk',
              'YPL_042915_001_RFPstack_083': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_001_RFPstack_083_skeleton.vtk',
              'YPL_042915_003_RFPstack_084': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_003_RFPstack_084_skeleton.vtk',
              'YPL_042915_003_RFPstack_085': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_003_RFPstack_085_skeleton.vtk',
              'YPL_042915_003_RFPstack_086': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_003_RFPstack_086_skeleton.vtk',
              'YPL_042915_005_RFPstack_087': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_005_RFPstack_087_skeleton.vtk',
              'YPL_042915_005_RFPstack_088': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_005_RFPstack_088_skeleton.vtk',
              'YPL_042915_005_RFPstack_089': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_005_RFPstack_089_skeleton.vtk',
              'YPL_042915_007_RFPstack_090': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_007_RFPstack_090_skeleton.vtk',
              'YPL_042915_007_RFPstack_091': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_007_RFPstack_091_skeleton.vtk',
              'YPL_042915_007_RFPstack_092': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_007_RFPstack_092_skeleton.vtk',
              'YPL_042915_009_RFPstack_093': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_009_RFPstack_093_skeleton.vtk',
              'YPL_042915_011_RFPstack_094': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_011_RFPstack_094_skeleton.vtk',
              'YPL_042915_013_RFPstack_097': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_013_RFPstack_097_skeleton.vtk',
              'YPL_042915_109_RFPStack_MM_098': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_109_RFPStack_MM_098_skeleton.vtk',
              'YPL_042915_109_RFPStack_MM_099': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_042915_109_RFPStack_MM_099_skeleton.vtk',
              'YPL_052315_001_RFPstack_000': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_001_RFPstack_000_skeleton.vtk',
              'YPL_052315_003_RFPstack_001': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_003_RFPstack_001_skeleton.vtk',
              'YPL_052315_003_RFPstack_002': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_003_RFPstack_002_skeleton.vtk',
              'YPL_052315_005_RFPstack_004': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_005_RFPstack_004_skeleton.vtk',
              'YPL_052315_005_RFPstack_005': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_005_RFPstack_005_skeleton.vtk',
              'YPL_052315_005_RFPstack_006': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_005_RFPstack_006_skeleton.vtk',
              'YPL_052315_007_RFPstack_007': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_007_RFPstack_007_skeleton.vtk',
              'YPL_052315_007_RFPstack_008': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_007_RFPstack_008_skeleton.vtk',
              'YPL_052315_009_RFPstack_010': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_009_RFPstack_010_skeleton.vtk',
              'YPL_052315_009_RFPstack_011': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_009_RFPstack_011_skeleton.vtk',
              'YPL_052315_009_RFPstack_012': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_009_RFPstack_012_skeleton.vtk',
              'YPL_052315_011_RFPstack_013': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_011_RFPstack_013_skeleton.vtk',
              'YPL_052315_011_RFPstack_014': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_011_RFPstack_014_skeleton.vtk',
              'YPL_052315_013_RFPstack_016': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_013_RFPstack_016_skeleton.vtk',
              'YPL_052315_013_RFPstack_017': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_013_RFPstack_017_skeleton.vtk',
              'YPL_052315_015_RFPstack_018': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_015_RFPstack_018_skeleton.vtk',
              'YPL_052315_017_RFPstack_021': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_017_RFPstack_021_skeleton.vtk',
              'YPL_052315_017_RFPstack_022': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_017_RFPstack_022_skeleton.vtk',
              'YPL_052315_017_RFPstack_023': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_017_RFPstack_023_skeleton.vtk',
              'YPL_052315_019_RFPstack_024': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_019_RFPstack_024_skeleton.vtk',
              'YPL_052315_019_RFPstack_025': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_019_RFPstack_025_skeleton.vtk',
              'YPL_052315_021_RFPstack_027': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_021_RFPstack_027_skeleton.vtk',
              'YPL_052315_021_RFPstack_028': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_021_RFPstack_028_skeleton.vtk',
              'YPL_052315_021_RFPstack_029': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_021_RFPstack_029_skeleton.vtk',
              'YPL_052315_023_RFPstack_030': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_023_RFPstack_030_skeleton.vtk',
              'YPL_052315_023_RFPstack_031': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_023_RFPstack_031_skeleton.vtk',
              'YPL_052315_023_RFPstack_032': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_023_RFPstack_032_skeleton.vtk',
              'YPL_052315_025_RFPstack_033': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_025_RFPstack_033_skeleton.vtk',
              'YPL_052315_025_RFPstack_034': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_025_RFPstack_034_skeleton.vtk',
              'YPL_052315_025_RFPstack_035': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_025_RFPstack_035_skeleton.vtk',
              'YPL_052315_025_RFPstack_036': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_025_RFPstack_036_skeleton.vtk',
              'YPL_052315_027_RFPstack_037': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_027_RFPstack_037_skeleton.vtk',
              'YPL_052315_027_RFPstack_038': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_027_RFPstack_038_skeleton.vtk',
              'YPL_052315_027_RFPstack_039': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_027_RFPstack_039_skeleton.vtk',
              'YPL_052315_029_RFPstack_040': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_029_RFPstack_040_skeleton.vtk',
              'YPL_052315_029_RFPstack_041': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_029_RFPstack_041_skeleton.vtk',
              'YPL_052315_029_RFPstack_042': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_029_RFPstack_042_skeleton.vtk',
              'YPL_052315_031_RFPstack_043': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_031_RFPstack_043_skeleton.vtk',
              'YPL_052315_031_RFPstack_044': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_031_RFPstack_044_skeleton.vtk',
              'YPL_052315_031_RFPstack_045': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_031_RFPstack_045_skeleton.vtk',
              'YPL_052315_033_RFPstack_046': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_033_RFPstack_046_skeleton.vtk',
              'YPL_052315_033_RFPstack_047': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPL\\Norm_YPL_052315_033_RFPstack_047_skeleton.vtk'},
             'YPR': {'YPR_042715_003_RFPstack_004': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_003_RFPstack_004_skeleton.vtk',
              'YPR_042715_003_RFPstack_006': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_003_RFPstack_006_skeleton.vtk',
              'YPR_042715_005_RFPstack_008': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_005_RFPstack_008_skeleton.vtk',
              'YPR_042715_007_RFPstack_010': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_007_RFPstack_010_skeleton.vtk',
              'YPR_042715_007_RFPstack_011': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_007_RFPstack_011_skeleton.vtk',
              'YPR_042715_007_RFPstack_012': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_007_RFPstack_012_skeleton.vtk',
              'YPR_042715_007_RFPstack_013': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_007_RFPstack_013_skeleton.vtk',
              'YPR_042715_009_RFPstack_014': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_009_RFPstack_014_skeleton.vtk',
              'YPR_042715_009_RFPstack_015': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_009_RFPstack_015_skeleton.vtk',
              'YPR_042715_011_RFPstack_017': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_011_RFPstack_017_skeleton.vtk',
              'YPR_042715_011_RFPstack_018': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_011_RFPstack_018_skeleton.vtk',
              'YPR_042715_011_RFPstack_019': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_011_RFPstack_019_skeleton.vtk',
              'YPR_042715_011_RFPstack_020': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_011_RFPstack_020_skeleton.vtk',
              'YPR_042715_011_RFPstack_021': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_011_RFPstack_021_skeleton.vtk',
              'YPR_042715_013_RFPstack_022': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_013_RFPstack_022_skeleton.vtk',
              'YPR_042715_013_RFPstack_023': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_013_RFPstack_023_skeleton.vtk',
              'YPR_042715_013_RFPstack_024': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_013_RFPstack_024_skeleton.vtk',
              'YPR_042715_015_RFPstack_026': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_015_RFPstack_026_skeleton.vtk',
              'YPR_042715_021_RFPstack_034': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_021_RFPstack_034_skeleton.vtk',
              'YPR_042715_021_RFPstack_035': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_021_RFPstack_035_skeleton.vtk',
              'YPR_042715_021_RFPstack_036': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_021_RFPstack_036_skeleton.vtk',
              'YPR_042715_023_RFPstack_037': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_023_RFPstack_037_skeleton.vtk',
              'YPR_042715_028_RFPstack_040': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_028_RFPstack_040_skeleton.vtk',
              'YPR_042715_028_RFPstack_041': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_028_RFPstack_041_skeleton.vtk',
              'YPR_042715_032_RFPstack_043': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_032_RFPstack_043_skeleton.vtk',
              'YPR_042715_032_RFPstack_044': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_032_RFPstack_044_skeleton.vtk',
              'YPR_042715_032_RFPstack_045': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_032_RFPstack_045_skeleton.vtk',
              'YPR_042715_032_RFPstack_047': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_032_RFPstack_047_skeleton.vtk',
              'YPR_042715_034_RFPstack_048': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_034_RFPstack_048_skeleton.vtk',
              'YPR_042715_034_RFPstack_049': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_034_RFPstack_049_skeleton.vtk',
              'YPR_042715_034_RFPstack_050': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_034_RFPstack_050_skeleton.vtk',
              'YPR_042715_036_RFPstack_051': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_036_RFPstack_051_skeleton.vtk',
              'YPR_042715_036_RFPstack_052': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_036_RFPstack_052_skeleton.vtk',
              'YPR_042715_036_RFPstack_053': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_036_RFPstack_053_skeleton.vtk',
              'YPR_042715_038_RFPstack_055': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_038_RFPstack_055_skeleton.vtk',
              'YPR_042715_038_RFPstack_056': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042715_038_RFPstack_056_skeleton.vtk',
              'YPR_042915_005_RFPstack_058': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_005_RFPstack_058_skeleton.vtk',
              'YPR_042915_007_RFPstack_062': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_007_RFPstack_062_skeleton.vtk',
              'YPR_042915_007_RFPstack_064': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_007_RFPstack_064_skeleton.vtk',
              'YPR_042915_009_RFPstack_067': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_009_RFPstack_067_skeleton.vtk',
              'YPR_042915_012_RFPstack_071': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_012_RFPstack_071_skeleton.vtk',
              'YPR_042915_012_RFPstack_072': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_012_RFPstack_072_skeleton.vtk',
              'YPR_042915_016_RFPstack_079': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_016_RFPstack_079_skeleton.vtk',
              'YPR_042915_016_RFPstack_080': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_016_RFPstack_080_skeleton.vtk',
              'YPR_042915_016_RFPstack_081': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_016_RFPstack_081_skeleton.vtk',
              'YPR_042915_016_RFPstack_082': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_016_RFPstack_082_skeleton.vtk',
              'YPR_042915_016_RFPstack_083': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_016_RFPstack_083_skeleton.vtk',
              'YPR_042915_020_RFPstack_086': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_020_RFPstack_086_skeleton.vtk',
              'YPR_042915_020_RFPstack_087': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_020_RFPstack_087_skeleton.vtk',
              'YPR_042915_020_RFPstack_088': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_020_RFPstack_088_skeleton.vtk',
              'YPR_042915_022_RFPstack_089': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_022_RFPstack_089_skeleton.vtk',
              'YPR_042915_022_RFPstack_090': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_022_RFPstack_090_skeleton.vtk',
              'YPR_042915_022_RFPstack_091': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_022_RFPstack_091_skeleton.vtk',
              'YPR_042915_022_RFPstack_092': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_022_RFPstack_092_skeleton.vtk',
              'YPR_042915_022_RFPstack_093': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_022_RFPstack_093_skeleton.vtk',
              'YPR_042915_026_RFPstack_097': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_026_RFPstack_097_skeleton.vtk',
              'YPR_042915_028_RFPstack_099': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_042915_028_RFPstack_099_skeleton.vtk',
              'YPR_052315_001_RFPstack_000': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_001_RFPstack_000_skeleton.vtk',
              'YPR_052315_001_RFPstack_001': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_001_RFPstack_001_skeleton.vtk',
              'YPR_052315_004_RFPstack_002': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_004_RFPstack_002_skeleton.vtk',
              'YPR_052315_004_RFPstack_003': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_004_RFPstack_003_skeleton.vtk',
              'YPR_052315_004_RFPstack_004': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_004_RFPstack_004_skeleton.vtk',
              'YPR_052315_004_RFPstack_005': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_004_RFPstack_005_skeleton.vtk',
              'YPR_052315_006_RFPstack_006': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_006_RFPstack_006_skeleton.vtk',
              'YPR_052315_006_RFPstack_007': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_006_RFPstack_007_skeleton.vtk',
              'YPR_052315_006_RFPstack_008': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_006_RFPstack_008_skeleton.vtk',
              'YPR_052315_006_RFPstack_009': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_006_RFPstack_009_skeleton.vtk',
              'YPR_052315_008_RFPstack_010': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_008_RFPstack_010_skeleton.vtk',
              'YPR_052315_008_RFPstack_011': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_008_RFPstack_011_skeleton.vtk',
              'YPR_052315_008_RFPstack_012': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_008_RFPstack_012_skeleton.vtk',
              'YPR_052315_008_RFPstack_013': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_008_RFPstack_013_skeleton.vtk',
              'YPR_052315_010_RFPstack_014': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_010_RFPstack_014_skeleton.vtk',
              'YPR_052315_012_RFPstack_015': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_012_RFPstack_015_skeleton.vtk',
              'YPR_052315_012_RFPstack_016': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_012_RFPstack_016_skeleton.vtk',
              'YPR_052315_014_RFPstack_017': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_014_RFPstack_017_skeleton.vtk',
              'YPR_052315_014_RFPstack_018': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_014_RFPstack_018_skeleton.vtk',
              'YPR_052315_014_RFPstack_019': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_014_RFPstack_019_skeleton.vtk',
              'YPR_052315_016_RFPstack_020': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_016_RFPstack_020_skeleton.vtk',
              'YPR_052315_018_RFPstack_021': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_018_RFPstack_021_skeleton.vtk',
              'YPR_052315_018_RFPstack_022': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_018_RFPstack_022_skeleton.vtk',
              'YPR_052315_020_RFPstack_023': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_020_RFPstack_023_skeleton.vtk',
              'YPR_052315_020_RFPstack_024': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_020_RFPstack_024_skeleton.vtk',
              'YPR_052315_020_RFPstack_025': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_020_RFPstack_025_skeleton.vtk',
              'YPR_052315_022_RFPstack_026': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_022_RFPstack_026_skeleton.vtk',
              'YPR_052315_022_RFPstack_027': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_022_RFPstack_027_skeleton.vtk',
              'YPR_052315_024_RFPstack_029': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_024_RFPstack_029_skeleton.vtk',
              'YPR_052315_024_RFPstack_030': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_024_RFPstack_030_skeleton.vtk',
              'YPR_052315_027_RFPstack_032': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_027_RFPstack_032_skeleton.vtk',
              'YPR_052315_029_RFPstack_033': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_029_RFPstack_033_skeleton.vtk',
              'YPR_052315_031_RFPstack_034': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_031_RFPstack_034_skeleton.vtk',
              'YPR_052315_033_RFPstack_035': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_033_RFPstack_035_skeleton.vtk',
              'YPR_052315_033_RFPstack_036': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_033_RFPstack_036_skeleton.vtk',
              'YPR_052315_033_RFPstack_037': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_033_RFPstack_037_skeleton.vtk',
              'YPR_052315_033_RFPstack_038': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_033_RFPstack_038_skeleton.vtk',
              'YPR_052315_033_RFPstack_039': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_033_RFPstack_039_skeleton.vtk',
              'YPR_052315_033_RFPstack_040': 'C:\\Users\\sweel_rafelski\\Desktop\\WorkingDir\\output\\normalizedVTK\\YPR\\Norm_YPR_052315_033_RFPstack_040_skeleton.vtk'}})

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

#            filename = op.join(datadir, '%s_grph.pkl' % mediatype)
#            with open(filename, 'wb') as output:
#                pickle.dump((nlist, elist, glist), output)
