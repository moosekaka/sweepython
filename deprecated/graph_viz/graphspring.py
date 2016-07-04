# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 13:21:38 2015
@author: sweel
    Calculate the statistical and topological measures of interest of cells
    in various carbon sources
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import cPickle as pickle
import networkx as nx
import seaborn as sns
from collections import defaultdict
import matplotlib.gridspec as gridspec
sns.set(style="white")
plt.close('all')


def drawneato(grph, graphname):
    """Draw a graph using pygraphwiz fructerman diagram style
    """
    n = int(np.ceil(np.sqrt(len(grph))))
    gs = gridspec.GridSpec(n, n)
    plt.figure(figsize=(11, 17))

    for ind, grh in enumerate(grph):
        cols = []
        sizes = []
        plt.suptitle(graphname, fontsize=16)
        for n, attr in grh.nodes(data=True):
            if attr['degree'] == 1:
                cols.append('#3366FF')
                sizes.append(10)
            else:
                cols.append('#FF5050')
                sizes.append(18)
        axes = plt.subplot(gs[ind])
        nx.draw_graphviz(grh,
                         # ‘neato’|’dot’|’twopi’|’circo’|’fdp’|
                         prog='neato',
                         ax=axes,
                         node_size=sizes,
                         node_color=cols)
    plt.show()
    plt.savefig('%s.png' % graphname)
    return (axes, n)


def retgrplist(graphs, gmea, lower, upper):
    """Return a subset of graphs whose indices fall within the lower and upper
    percentiles of a graph connectivity measure contained in gmea
    """
#    filt = np.percentile(gmea, [lower, upper])
#    ind_filt = [ind for ind, ele in enumerate(gmea) if
#                ele > filt[0] and
#                ele < filt[1]]
#    graphs_filt = [graphs[el] for el in ind_filt]
    gdic = {i.graph['cell']: i for i in graphs}

    filt = gmea[(gmea > gmea.quantile(lower)) &
                (gmea < gmea.quantile(upper))]
    graphs_filt = [gdic[el] for el in filt.index]

    return graphs_filt

# =============================================================================
#       init vars and get vtk, graph data
# =============================================================================
# pylint: disable=C0103
#with open('statsResults.pkl', 'rb') as inpt:
#    (media, mitoEdges, mitoLength, rfpWidthCorCoef,
#     mitoEdgeNorm, mitoStdNorm, cellNormDY, cellNormDY_raw,
#     mitoEdgeNorm_unscaled, mitoStdNorm_unscaled,
#     centralityBP, closeCentralityBP, ClstCoef, charPathL,
#     phi, pk3, kNN3, AvgDeg,
#     centralityBP_W, closeCentralityBP_W, ClstCoef_W, charPathL_W, kNN3_W) =\
#     pickle.load(inpt)
dirList = []
Graphs = defaultdict(dict)
# pylint: enable=C0103

for root, dirs, files in os.walk(os.getcwd()):
    for f in dirs:
        if f.startswith('YP'):
            dirList.append(
                os.path.join(root, f))


with open('munged_dataframe.pkl', 'rb') as INPUT:
    df = pickle.load(INPUT)
with open('lagedges.pkl', 'rb') as INPT:
    dflags = pickle.load(INPT)

df['lags_1'] = dflags

# =============================================================================
#       Main
# =============================================================================
for media in dirList:
    labs = media[-3:]
    print'Now on %s' % labs+"\n"+"="*79
    Graphs[labs] = []
    with open(os.path.join(media,
                           '%s_grph.pkl' % labs), 'rb') as inpt:
        G = pickle.load(inpt)[2]
        for h, i in enumerate(G):
            Graphs[labs].append(G[h])

filt = df.loc[:, ['mito_knn_uw', 'media']]

YPD = retgrplist(Graphs['YPD'],
                 filt[filt['media'] == 'YPD'].ix[:, 0].apply(np.mean),
                 .2, .22)
YPE = retgrplist(Graphs['YPE'],
                 filt[filt['media'] == 'YPE'].ix[:, 0].apply(np.mean),
                 .90, .91)
YPL = retgrplist(Graphs['YPL'],
                 filt[filt['media'] == 'YPL'].ix[:, 0].apply(np.mean),
                 .90, .913)
YPR = retgrplist(Graphs['YPR'],
                 filt[filt['media'] == 'YPR'].ix[:, 0].apply(np.mean),
                 .52, .53)

gs = gridspec.GridSpec(6, 4)
G = [YPD, YPE, YPL, YPR]
labels = ['YPD', 'YPE', 'YPL', 'YPR']
for ind, grh in enumerate(G):

    for idx, graph in enumerate(grh):

        cols = []
        sizes = []

        for n, attr in graph.nodes(data=True):
            if attr['degree'] == 1:
                cols.append('#3366FF')
                sizes.append(10)
            else:
                cols.append('#FF5050')
                sizes.append(18)
        axes = plt.subplot(gs[idx, ind])
        if idx == 0:
            axes.set_title(labels[ind])

        nx.draw_graphviz(graph,
                         # ‘neato’|’dot’|’twopi’|’circo’|’fdp’|
                         prog='neato',
                         ax=axes,
                         node_size=sizes,
                         node_color=cols)
        plt.suptitle(filt.columns[0])
plt.show()

#==============================================================================
# draw single graph
#==============================================================================
cols = []
sizes = []
fig, axes = plt.subplots(1, 1, figsize=(8, 6))
GYPL = {i.graph['cell']: i for i in Graphs['YPL']}
graph = GYPL['YPL_042915_005_RFPstack_089']
subgr = nx.Graph()
for n, nbrs in graph.adjacency_iter():
    for nbr, attr in nbrs.items():
        maxvalue = max([d['weight'] for d in attr.values()])
        subgr.add_edge(n, nbr, weight=maxvalue)

for n, attr in graph.nodes(data=True):
    if attr['degree'] == 1:
        cols.append('#3366FF')
        sizes.append(30)
    else:
        cols.append('#FF5050')
        sizes.append(50)

nx.draw_graphviz(subgr,
                 # ‘neato’|’dot’|’twopi’|’circo’|’fdp’|
                 prog='neato',
                 ax=axes,
                 node_size=sizes,
                 node_color=cols)
pos = nx.graphviz_layout(subgr)
labels = nx.average_neighbor_degree(subgr)
labs = {i: '%4.3f' % labels[i] for i in labels}
#nx.draw_networkx_labels(subgr, pos=pos, labels=labs)
