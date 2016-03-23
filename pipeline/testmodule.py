# -*- coding: utf-8 -*-
"""
test modules for vectorized implementation of make_networkx
@author: sweel_rafelski
"""
import os
import os.path as op
import cPickle as pickle
import numpy as np
import wrappers as wr
# pylint: disable=C0103
datadir = op.join(os.getcwd(), 'data')
rawdir = op.join(os.getcwd(), 'output')

glist = wr.swalk(datadir, '*grph.pkl')
testlist = wr.swalk(datadir, '*grph_b.pkl')

# List from vectorized method to be tested
for k in sorted(glist)[:]:
    with open(glist[k], 'rb') as ipt:
        nlist, elist, nx = pickle.load(ipt)
    tot = 0
    tot2 = 0
    tot3 = 0
    tot4 = 0
    for i in nlist:
        for j in i:
            tot = tot + j[1]['inten']
            tot2 = tot2 + j[1]['degree']
    for i in elist:
        for j in i:
            tot3 = tot3 + j[2]['weight']
            tot4 = tot4 + j[2]['cellID']
        test = [i[0] for i in elist]
        sumt = [t[2]['weight'] for t in test]

    ntest1 = []
    for ind, el in enumerate(elist):
        nlist_d = {}
        for node, coord in nlist[ind]:
            nlist_d[node] = coord

        ntest1.append([(nlist_d[n1]['coord'], nlist_d[n2]['coord'])
                       for n1, n2, eatt in el if eatt['cellID'] ==
                       len(el)%2][0])
    x1 = np.array(ntest1)
    x1 = np.sum(x1, axis=(2, 1, 0))
    print ("total intensities for {} = {:,.2f}\ndegree= {:,}, "
           "edw= {:,.2f}, # of edges= {:,}\n"
           "total for coords = {:,.4f}\n".format(k, tot, tot2, tot3, tot4, x1))

# List data from unvectorized method to test against
for k in sorted(testlist)[:]:
    with open(testlist[k], 'rb') as ipt:
        nlist_t, elist_t, nx_t = pickle.load(ipt)
    tot = 0
    tot2 = 0
    tot3 = 0
    tot4 = 0
    for i in nlist_t:
        for j in i:
            tot = tot + j[1]['inten']
            tot2 = tot2 + j[1]['degree']
    for i in elist_t:
        for j in i:
            tot3 = tot3 + j[2]['weight']
            tot4 = tot4 + j[2]['cellID']
        test2 = [i[0] for i in elist_t]
        sumt2 = [t[2]['weight'] for t in test2]

    ntest2 = []
    for ind, el in enumerate(elist_t):
        nlist_td = {}
        for node, coord in nlist_t[ind]:
            nlist_td[node] = coord
        ntest2.append([(nlist_td[n1]['coord'], nlist_td[n2]['coord'])
                       for n1, n2, eatt in el if eatt['cellID'] ==
                       len(el)%2][0])
    x2 = np.array(ntest2)
    x2 = np.sum(x2, axis=(2, 1, 0))
    print ("total intensities for {} = {:,.2f}\ndegree= {:,}, "
           "edw= {:,.2f}, # of edges= {:,}\n"
           "total for coords = {:,.4f}\n".format(k, tot, tot2, tot3, tot4, x2))
