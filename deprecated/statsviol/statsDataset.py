"""
    Calculate the statistical and topological measures of interest of cells
    in various carbon sources
"""
import os
import fnmatch
import numpy as np
import cPickle as pickle
from tvtk.api import tvtk
import networkx as nx
import scipy.stats as sp
import seaborn as sns
sns.set(style="white")
# =============================================================================
#               init vars and get vtk, graph data
# =============================================================================
# pylint: disable=C0103
vtkF = {}
G = {}
backgroundGFP = {}
backgroundRFP = {}
parDir = os.path.dirname(os.getcwd())
for root, dirs, files in os.walk(os.getcwd()):
    for i in files:
        if fnmatch.fnmatch(i, '*vtk'):
            vtkF.setdefault(root.rsplit('\\', 1)[1], []).append(
                os.path.join(root, i))
        if fnmatch.fnmatch(i, '*grph.pkl'):
            G.setdefault(root.rsplit('\\', 1)[1], []).append(
                os.path.join(root, i))

media = sorted(vtkF.keys())
# pylint: enable=C0103

#    networkX graph objects of mitograph
for i in G:
    with open(G[i][0], 'rb') as inpt:
        temp = pickle.load(inpt)[2]
        G[i].append(temp)

#    get metadatas
with open(parDir+'\\'+'fileMetas.pkl', 'rb') as inpt:
    METAS = pickle.load(inpt)
for i in METAS:
    backgroundRFP[i] = METAS[i][0]
    backgroundGFP[i] = METAS[i][1]

# =============================================================================
#               begin pipeline
# =============================================================================
# pylint: disable=C0103
mitoLength = {}
mitoNodes = {}
mitoEdgeNorm = {}
mitoStdNorm = {}
mitoEdgeNorm_unscaled = {}
mitoStdNorm_unscaled = {}
cellNormDY = {}
cellNormDY_raw = {}
mitoEdges = {}
mitoNodesAll = {}
centralityBP = {}
closeCentralityBP = {}
rfpWidthCorCoef = {}
kNN3 = {}
ClstCoef = {}
charPathL = {}
AvgDeg = {}
phi = {}
pk3 = {}
A = {}
B = {}
centralityBP_W = {}
closeCentralityBP_W = {}
kNN3_W = {}
charPathL_W = {}
ClstCoef_W = {}
# pylint: enable=C0103
for mem in media:
    print'%s\n' % mem + "-" * 79
    for n, a in enumerate(vtkF[mem]):
        Norm = []
        NormRaw = []
        GFP = []
        RFP = []
        arrPts = []
        W = []
        W2 = []
        rGFP = []
        lineId = {}

        curGrph = G[mem][1][n]
        reader = tvtk.PolyDataReader()
        reader.set(file_name=a)
        reader.update()
        data = reader.output
        scalarsNorm = data.point_data.scalars
        temp = data.point_data
        dyRaw = np.ravel(temp.get_array('DY_raw'))
        rawGFP = np.ravel(temp.get_array('rGFP'))
        rawRFP = np.ravel(temp.get_array('rRFP'))
        WidthEq = np.ravel(temp.get_array('WidthEq'))
        tubeWidth = np.ravel(temp.get_array('tubeWidth'))
        fileKey = a.rsplit('\\', 1)[1][5:][:-13]

        if backgroundRFP[fileKey] > min(rawRFP):
            minA = backgroundRFP[fileKey]-1
        else:
            minA = backgroundRFP[fileKey]-1
        minB = min(backgroundGFP[fileKey], min(rawGFP))

#    Convert multigraph to graph for clustering coef calc (choose long edges)
        GG = nx.Graph()
        for n, nbrs in curGrph.adjacency_iter():
            for nbr, attr in nbrs.items():
                maxvalue = max([d['weight'] for d in attr.values()])
                GG.add_edge(n, nbr, weight=maxvalue)

#    get the line of interest
        bpts = [curGrph.node[i]['coord'] for i in curGrph.nodes()
                if curGrph.node[i]['degree'] > 1]
        bptsId = [data.find_point(el) for el in bpts]

        epts = [curGrph.node[i]['coord'] for i in curGrph.nodes()
                if curGrph.node[i]['degree'] == 1]
        eptsId = [data.find_point(el) for el in epts]

        for j in range(data.number_of_lines):
            cellIds = list(data.get_cell(j).point_ids)
            Norm.append([scalarsNorm[i] for i in cellIds])
            NormRaw.append([dyRaw[i] for i in cellIds])
            GFP.append([rawGFP[i] for i in cellIds])
            RFP.append([rawRFP[i] for i in cellIds])
            W.append([WidthEq[i] for i in cellIds])

        t = curGrph.edges(data=True)
        x = curGrph.nodes(data=True)
        y = nx.betweenness_centrality(curGrph)
        z = nx.closeness_centrality(curGrph)
        knn3 = nx.k_nearest_neighbors(curGrph)[3]  # bpts kNN
        cc = nx.average_clustering(GG)
        k3 = len([i[0] for i in x if i[1]['degree'] > 1])
        # weighted by edgelens versions
        yW = nx.betweenness_centrality(curGrph, weight='weight')
        zW = nx.closeness_centrality(curGrph, distance='weight')
        knn3W = nx.k_nearest_neighbors(curGrph, weight='weight')[3]
        ccW = nx.average_clustering(GG, weight='weight')
# =============================================================================
#               these are the stats of interest
# =============================================================================
#       Topology
        mitoNodes.setdefault(mem, []).append(
            [i[1]['degree'] for i in x if i[1]['degree'] > 2])

        mitoNodesAll.setdefault(mem, []).append(
            [i[1]['degree'] for i in x])

        mitoEdges.setdefault(mem, []).append(
            [i[2]['weight'] for i in t])

        mitoLength.setdefault(mem, []).append(
            np.sum([i[2]['weight'] for i in t]))

#       Function
        mitoEdgeNorm.setdefault(mem, []).append(  # per cell!
            [np.mean(el) for el in Norm])

        mitoStdNorm.setdefault(mem, []).append(  # per cell!
            [np.std(el) for el in Norm])

        cellNormDY.setdefault(mem, []).append(np.mean(scalarsNorm))

        cellNormDY_raw.setdefault(mem, []).append(np.mean(dyRaw))

        mitoEdgeNorm_unscaled.setdefault(mem, []).append(  # per cell!
            [np.mean(el) for el in NormRaw])

        mitoStdNorm_unscaled.setdefault(mem, []).append(  # per cell!
            [np.std(el) for el in NormRaw])

#       Connectivity
        rfpWidthCorCoef.setdefault(mem, []).append(
            sp.pearsonr(rawRFP, tubeWidth))

        centralityBP.setdefault(mem, []).append(
            [y[i[0]] for i in x if i[1]['degree'] > 2])  # per cell!

        closeCentralityBP.setdefault(mem, []).append(
            [z[i[0]] for i in x if i[1]['degree'] > 2])  # per cell

        kNN3.setdefault(mem, []).append(knn3)  # nearest neighb conn

        ClstCoef.setdefault(mem, []).append(cc)

        charPathL.setdefault(mem, []).append(
            np.max(
                [nx.average_shortest_path_length(g) for g
                 in nx.connected_component_subgraphs(curGrph)
                 if len(g) > 1]))

#       Connectivity weighted by length
        centralityBP_W.setdefault(mem, []).append(
            [yW[i[0]] for i in x if i[1]['degree'] > 2])  # per cell!

        closeCentralityBP_W.setdefault(mem, []).append(
            [zW[i[0]] for i in x if i[1]['degree'] > 2])  # per cell

        kNN3_W.setdefault(mem, []).append(knn3W)  # nearest neighb conn

        ClstCoef_W.setdefault(mem, []).append(ccW)

        charPathL_W.setdefault(mem, []).append(
            np.max(
                [nx.average_shortest_path_length(g, weight='weight') for g
                 in nx.connected_component_subgraphs(curGrph)
                 if len(g) > 1]))

        phi.setdefault(mem, []).append(
            1. * np.max(
                [g.number_of_edges() for g
                 in nx.connected_component_subgraphs(curGrph)]) /
            curGrph.number_of_edges())

        pk3.setdefault(mem, []).append(
            1. * k3 / curGrph.number_of_nodes())

    #   for average degree
    B.setdefault(mem, []).append(
        [len(i) for i in mitoNodesAll[mem]])
    A.setdefault(mem, []).append(
        [2.*len(j) for j in mitoEdges[mem]])
    AvgDeg.setdefault(mem, []).append(
        [A[mem][0][k] / B[mem][0][k] for k in range(len(A[mem][0]))])

OUT = (media, mitoEdges, mitoLength, rfpWidthCorCoef,
       mitoEdgeNorm, mitoStdNorm, cellNormDY, cellNormDY_raw,
       mitoEdgeNorm_unscaled, mitoStdNorm_unscaled,
       centralityBP, closeCentralityBP, ClstCoef, charPathL,
       phi, pk3, kNN3, AvgDeg,
       centralityBP_W, closeCentralityBP_W, ClstCoef_W, charPathL_W, kNN3_W)
with open('statsResults.pkl', 'wb') as OUTPUT:
    pickle.dump(OUT, OUTPUT)
