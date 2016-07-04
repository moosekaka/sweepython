# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 00:25:17 2015
Draw a toygraph with labels for illustration and example
@author: sweel
"""
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('dark')
plt.close('all')

#==============================================================================
# elist = [('a', 'b', 1.0),
#          ('b', 'c', 1.0),
#          ('a', 'c', 2.0),
#          # ('d', 'e', 2.0),
#          # ('d', 'f', 2.0),
#          # ('b', 'e', 2.0),
#          ('a', 'd', 2.0)]

elist =[
         (2867, 2869, 1.5),
         (2868, 2869, 0.9),
         (2870, 2880, 2.5),
         (2871, 2881, 2.1),
         (2871, 2875, 1.2),
         (2874, 2875, 0.3),
         (2875, 2885, 0.9),
         (2876, 2876, 2.4),
         (2876, 2878, 1.9),
         (2877, 2878, 0.2),
         (2878, 2880, 0.4),
         (2879, 2883, 0.5),
         (2880, 2883, 0.8),
         (2881, 2884, 1.7),
         (2881, 2885, 0.6),
         (2882, 2885, 0.4),
         (2883, 2884, 1.3),
         (2884, 2885, 1.6) ]

#==============================================================================

G = nx.Graph()
G.add_weighted_edges_from(elist)

elab = {}
for a, b, c in G.edges(data=True):
    elab[(a, b)] = '%s' % c['weight']

nlab = {}
for key, val in nx.betweenness_centrality(G,
                                          weight=False).items():
    nlab[key] = val

nlab2 = {}
for key, val in nx.closeness_centrality(G,
                                        distance=False).items():
    nlab2[key] = val

nlab3 = {}
for key, val in nx.betweenness_centrality(G,
                                          normalized=True,
                                          weight='weight').items():
    nlab3[key] = val

nlab4 = {}
for key, val in nx.closeness_centrality(G,
                                        distance='weight').items():
    nlab4[key] = val

fig, ax1 = plt.subplots(1, 1)
pose = nx.graphviz_layout(G)

nlab5 = {}
for key, val in nx.clustering(G,
                                        weight=False).items():
    nlab5[key] = val


nlab6 = {}
for key, val in nx.clustering(G,
                                        weight='weight').items():
    nlab6[key] = val

nx.draw_networkx(G,
                 pos=pose,
                 ax=ax1,
                 node_size=300,
                 with_labels=True)

nx.draw_networkx_edge_labels(G,
                             pos=pose,
                             edge_labels=elab,
                             ax=ax1)

for k, v in pose.items():

#    plt.text(v[0]+2,
#             v[1]+18,
#             s='betw cent: %4.2f' % nlab[k],
#             bbox={'facecolor': 'red', 'alpha': 0.5},
#             horizontalalignment='center')

    plt.text(v[0]+2,
             v[1]+18,
             s='close cent: %4.2f' % nlab2[k],
             bbox={'facecolor': 'g', 'alpha': 0.5},
             horizontalalignment='center')
#
#    plt.text(v[0]+2,
#             v[1]-18,
#             s='betw cent_w: %4.2f' % nlab3[k],
#             bbox={'facecolor': 'c', 'alpha': 0.5},
#             horizontalalignment='center')
#
    plt.text(v[0]+2,
             v[1]-20,
             s='close cent_w: %4.2f' % nlab4[k],
             bbox={'facecolor': 'm', 'alpha': 0.5},
             horizontalalignment='center')

#    plt.text(v[0]+2,
#             v[1]+15,
#             s='loc_cluster: %4.2f' % nlab5[k],
#             bbox={'facecolor': 'b', 'alpha': 0.5},
#             horizontalalignment='center')
#
#    plt.text(v[0]+2,
#             v[1]-20,
#             s='loc_cluster_w: %4.2f' % nlab6[k],
#             bbox={'facecolor': 'orange', 'alpha': 0.5},
#             horizontalalignment='center')
