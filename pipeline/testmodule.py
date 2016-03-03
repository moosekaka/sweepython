# -*- coding: utf-8 -*-
"""
test modules
@author: sweel_rafelski
"""
import os
import os.path as op
import cPickle as pickle
import wrappers as wr
# pylint: disable=C0103
datadir = op.join(os.getcwd(), 'data')
rawdir = op.join(os.getcwd(), 'output')

glist = wr.swalk(datadir, '*grph.pkl')
testlist = wr.swalk(datadir, '*grph_b.pkl')


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
            tot4 = tot4 + len(i)
    print "total inten for %s is %.2f, degree is %d, edw %.2f, #edges %d\n" % (k, tot, tot2, tot3, tot4)

for k in sorted(testlist):
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
            tot4 = tot4 + len(i)
    print "total inten for %s is %.2f, degree is %d, edw %.2f, , #edges %d\n" % (k, tot, tot2, tot3, tot4)
