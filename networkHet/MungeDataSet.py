"""
    Calculate the statistical and topological measures of interest of cells
    in various carbon sources and munge it into dataframe
"""
import os
import fnmatch
import numpy as np
from tvtk.api import tvtk
import networkx as nx
import scipy.stats as sp
import matplotlib.pyplot as plt
from collections import defaultdict
from mungedata import MungeDataFuncs as md
import pandas as pd
import cPickle as pickle
from numpy.random import choice as samp_no_rep
plt.close('all')

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
        if fnmatch.fnmatch(i, 'N*vtk'):
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
mito_avgdeg = defaultdict(dict)
mito_bpts_dy = defaultdict(dict)
mito_bpts_dyraw = defaultdict(dict)
mito_bptcoefvar_raw = defaultdict(dict)
mito_btwcntr_uw = defaultdict(dict)
mito_btwcntr_w = defaultdict(dict)
mito_cell_avedy = defaultdict(dict)
mito_cell_avedyr = defaultdict(dict)
mito_cell_stddy = defaultdict(dict)
mito_cell_stddyr = defaultdict(dict)
mito_charpl_uw = defaultdict(dict)
mito_charpl_w = defaultdict(dict)
mito_clscntr_uw = defaultdict(dict)
mito_clscntr_w = defaultdict(dict)
mito_clstcf_uw = defaultdict(dict)
mito_clstcf_w = defaultdict(dict)
mito_edge_avedy = defaultdict(dict)
mito_edge_avedyr = defaultdict(dict)
mito_edge_coefvar = defaultdict(dict)
mito_edge_coefvarr = defaultdict(dict)
mito_edge_stddy = defaultdict(dict)
mito_edge_stddyr = defaultdict(dict)
mito_edgelen = defaultdict(dict)
mito_edgenum = defaultdict(dict)
mito_knn_uw = defaultdict(dict)
mito_knn_w = defaultdict(dict)
mito_beta_top = defaultdict(dict)
mito_beta_geo = defaultdict(dict)
mito_phi = defaultdict(dict)
mito_pk3 = defaultdict(dict)
mito_totlen = defaultdict(dict)
mito_widcoef = defaultdict(dict)
mito_widcoefDY = defaultdict(dict)
mito_cell_ave_gfp = defaultdict(dict)
mito_cell_ave_rfp = defaultdict(dict)
mito_cell_w = defaultdict(dict)
mito_iso_dyr = defaultdict(dict)
mito_bootbpts_dyraw = defaultdict(dict)
mito_tubew = defaultdict(dict)
mito_nodenumbers = defaultdict(dict)

cnnsub = nx.connected_component_subgraphs
avg_shpthl = nx.average_shortest_path_length
avg_nnd = nx.average_neighbor_degree
# pylint: enable=C0103
for mem in media:
    print'\nNow on %s\n' % mem + "=" * 79
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
        filekey = a.rsplit('\\', 1)[1][5:][:-13]

        if backgroundRFP[filekey] > min(rawRFP):
            minA = backgroundRFP[filekey]-1
        else:
            minA = backgroundRFP[filekey]-1
        minB = min(backgroundGFP[filekey], min(rawGFP))

#    Convert multigraph to graph for clustering coef calc (choose long edges)
        GG = nx.Graph()
        for n, nbrs in curGrph.adjacency_iter():
            for nbr, attr in nbrs.items():
                maxvalue = max([d['weight'] for d in attr.values()])
                GG.add_edge(n, nbr, weight=maxvalue)

#    get btps intensity value within radofinfluence
        branchpoints = {j: attr['coord'] for j, attr
                        in curGrph.nodes(data=True)
                        if attr['degree'] > 2}

        bptpid = md.bpts_inten(data, branchpoints)
        bptdy = {key: np.mean([scalarsNorm[el] for el in vals])
                 for key, vals in sorted(bptpid.iteritems())}

        bptdy_raw = {key: np.mean([dyRaw[el] for el in vals])
                     for key, vals in sorted(bptpid.iteritems())}

        bptcoefvar_raw = {key: sp.variation([dyRaw[el] for el in vals])
                          for key, vals in sorted(bptpid.iteritems())}

#    make bootstrapped btps and bootstrap dyraw around rad of influence
        bpids = np.unique([el for lis in bptpid.values() for el in lis])
        allids = {i: np.array(data.get_cell(i).point_ids) for i
                  in range(data.number_of_lines)}
        nonbpids = np.unique([el for lis in allids.values() for el
                              in lis if el not in bpids])

        bootbp = defaultdict(dict)
        nboot = 100  # number of replicates for bootstrap
        for n in range(nboot):
            bootbp[n] = {}
            for k, v in bptpid.items():
                bootbp[n][k] = samp_no_rep(nonbpids,
                                           size=len(v),
                                           replace=False)
        mean_bs = {}
        for key in bptpid.keys():
            mean_bs[key] = [np.mean([dyRaw[el] for el
                                     in bootbp[n][key]]) for n
                            in range(nboot)]
        bootbpdy_raw = {key: np.mean(vals) for key, vals
                        in mean_bs.iteritems()}

#    get the line of interest
        for line in range(data.number_of_lines):
            cellIds = list(data.get_cell(line).point_ids)
            Norm.append([scalarsNorm[cid] for cid in cellIds])
            NormRaw.append([dyRaw[cid] for cid in cellIds])
            GFP.append([rawGFP[cid] for cid in cellIds])
            RFP.append([rawRFP[cid] for cid in cellIds])
            W.append([WidthEq[cid] for cid in cellIds])

        dedges = {(a, b, eattr['cellID']): eattr['weight'] for a, b, eattr
                  in curGrph.edges(data=True)}  # edgelists
        node_btwcent = nx.betweenness_centrality(curGrph)
        node_clscent = nx.closeness_centrality(curGrph)
        anndeg = avg_nnd(curGrph, nodes=branchpoints)  # bpts nearest ngbr deg
        lc = nx.clustering(GG)
        lcW = nx.clustering(GG, weight='weight')
        cc = nx.average_clustering(GG)  # clus. coef of graph
        k3 = len(branchpoints)  # number of nodes deg == 3
        # weighted by edgelens versions
        node_btwcentW = nx.betweenness_centrality(curGrph,
                                                  weight='weight')
        node_clscentW = nx.closeness_centrality(curGrph,
                                                distance='weight')
        anndeg_w = avg_nnd(curGrph, nodes=branchpoints,
                           weight='weight')
        ccW = nx.average_clustering(GG, weight='weight')

        conncomps = [cnn for cnn in cnnsub(curGrph)]
        largest_cnn = conncomps[np.argmax([g.number_of_edges() for g
                                           in cnnsub(curGrph)])]
#        smallest_cnn = conncomps[np.argmin([g.number_of_edges() for g
#                                            in cnnsub(curGrph)])]
        isocnn = [a for a in conncomps if len(a.edges()) == 1]

#        isoedgecid = [eattr['cellID'] for n1, n2, eattr
#                      in smallest_cnn.edges(data=True)][0]
#        isoedgepid = data.get_cell(isoedgecid).point_ids
        isoedgecid = [eattr['cellID'] for subg
                      in isocnn for n1, n2, eattr
                      in subg.edges(data=True)]
        isoedgepid = {cid: np.array(data.get_cell(cid).point_ids) for cid
                      in isoedgecid}

# =============================================================================
#              these are the stats of interest
# =============================================================================
#       Topology
        # make sure the edges are returned in the same order as DY (by cellid)

        x = {i[2]: dedges[i] for i in dedges.keys()}

        mito_edgelen[filekey] = [x[key] for key in sorted(x.keys())]

        mito_totlen[filekey] = np.sum(mito_edgelen[filekey])

        mito_edgenum[filekey] = curGrph.number_of_edges()
        mito_nodenumbers[filekey]= curGrph.number_of_nodes()

        for line in isoedgepid:
            mito_iso_dyr[filekey][line] = [dyRaw[pid] for pid
                                           in isoedgepid[line]]
#       Function
        mito_edge_avedy[filekey] = [
            np.mean(el) for el in Norm]  # per cell

        mito_edge_stddy[filekey] = [
            np.std(el) for el in Norm]  # per cell!

        mito_edge_avedyr[filekey] = [
            np.mean(el) for el in NormRaw]

        mito_edge_stddyr[filekey] = [
            np.std(el) for el in NormRaw]

        mito_edge_coefvar[filekey] = [
            np.std(el)/np.mean(el) for el in Norm]

        mito_edge_coefvarr[filekey] = [
            np.std(el)/np.mean(el) for el in NormRaw]

        mito_cell_avedy[filekey] = np.mean(scalarsNorm)

        mito_cell_avedyr[filekey] = np.mean(dyRaw)

        mito_cell_stddy[filekey] = np.std(scalarsNorm)

        mito_cell_stddyr[filekey] = np.std(dyRaw)

        mito_bpts_dy[filekey] = [bptdy[key] for key
                                 in sorted(branchpoints)]

        mito_bpts_dyraw[filekey] = [bptdy_raw[key] for key
                                    in sorted(branchpoints)]
        mito_bptcoefvar_raw[filekey] = [bptcoefvar_raw[key] for key
                                        in sorted(branchpoints)]
        mito_bootbpts_dyraw[filekey] = [bootbpdy_raw[key] for key
                                        in sorted(branchpoints)]

        mito_cell_ave_gfp[filekey] = np.mean(rawGFP-minB)
        mito_cell_ave_rfp[filekey] = np.mean(rawRFP-minA)
        mito_cell_w[filekey] = np.mean(WidthEq)
        mito_tubew[filekey] = np.mean(tubeWidth)
#       Connectivity
        mito_widcoef[filekey] = sp.pearsonr(rawRFP, tubeWidth)
        mito_widcoefDY[filekey] = sp.pearsonr(dyRaw, tubeWidth)

        mito_btwcntr_uw[filekey] = [node_btwcent[key] for key
                                    in sorted(branchpoints)]  # per cell!

        mito_clscntr_uw[filekey] = [node_clscent[key] for key
                                    in sorted(branchpoints)]  # per cell

        mito_knn_uw[filekey] = [anndeg[key] for key  # nearest neighb conn
                                in sorted(branchpoints)]

        mito_clstcf_uw[filekey] = [lc[key] for key
                                   in sorted(branchpoints)]

        mito_charpl_uw[filekey] = np.max(
            [avg_shpthl(g) for g in cnnsub(curGrph)
             if len(g) > 1])

#       Connectivity weighted by length
        mito_btwcntr_w[filekey] = [node_btwcentW[key] for key
                                   in sorted(branchpoints)]  # per cell!

        mito_clscntr_w[filekey] = [node_clscentW[key] for key
                                   in sorted(branchpoints)]  # per cell

        mito_knn_w[filekey] = [anndeg_w[key] for key
                               in sorted(branchpoints)]  # nrst nb. conn.

        mito_clstcf_w[filekey] = [lcW[key] for key
                                  in sorted(branchpoints)]

        mito_charpl_w[filekey] = np.max(
            [avg_shpthl(g, weight='weight') for g in cnnsub(curGrph)
             if len(g) > 1])

        mito_phi[filekey] = (1. * largest_cnn.number_of_nodes() /
                             curGrph.number_of_nodes())

        mito_beta_top[filekey] = (1. * largest_cnn.number_of_edges() /
                                  curGrph.number_of_edges())

        mito_beta_geo[filekey] = (np.sum([eattr['weight'] for e, f, eattr
                                          in largest_cnn.edges(data=True)]) /
                                  mito_totlen[filekey])

        mito_pk3[filekey] = 1. * k3 / curGrph.number_of_nodes()

        mito_avgdeg[filekey] = (2. * curGrph.number_of_edges() /
                                curGrph.number_of_nodes())


OUT = ('mito_avgdeg',
       'mito_beta_geo',
       'mito_beta_top',
       'mito_bpts_dy',
       'mito_bpts_dyraw',
       'mito_bptcoefvar_raw',
       'mito_btwcntr_uw',
       'mito_btwcntr_w',
       'mito_cell_avedy',
       'mito_cell_avedyr',
       'mito_cell_stddy',
       'mito_cell_stddyr',
       'mito_charpl_uw',
       'mito_charpl_w',
       'mito_clscntr_uw',
       'mito_clscntr_w',
       'mito_clstcf_uw',
       'mito_clstcf_w',
       'mito_edge_avedy',
       'mito_edge_avedyr',
       'mito_edge_coefvar',
       'mito_edge_coefvarr',
       'mito_edge_stddy',
       'mito_edge_stddyr',
       'mito_edgelen',
       'mito_edgenum',
       'mito_knn_uw',
       'mito_knn_w',
       'mito_phi',
       'mito_pk3',
       'mito_totlen',
       'mito_widcoef',
       'mito_widcoefDY',
       'mito_cell_ave_gfp',
       'mito_cell_ave_rfp',
       'mito_cell_w',
       'mito_iso_dyr',
       'mito_bootbpts_dyraw',
       'mito_tubew',
       'mito_nodenumbers')

FRAMES = []
for i in OUT:
    FRAMES.append(pd.Series(globals()[i], name=i))
df = pd.concat(FRAMES, axis=1)
df['cat'] = df.index
df['media'] = df['cat'].apply(lambda x: x.split('_', 1)[0])
df = df.drop('cat', axis=1)
df['charpl_norm_len'] = df.mito_charpl_uw / df.mito_totlen
df['charpl_norm_numedge'] = df.mito_charpl_w / df.mito_edgenum
df['cell_coefvar'] = df['mito_edge_avedy'].apply(sp.variation)
df['cell_coefvar_r'] = df['mito_edge_avedyr'].apply(sp.variation)

#with open('munged_dataframe.pkl', 'wb') as OUTPUT:
#    pickle.dump(df, OUTPUT)
