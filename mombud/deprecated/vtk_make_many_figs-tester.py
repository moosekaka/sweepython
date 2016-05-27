# -*- coding: utf-8 -*-
"""
Batch visualize skel and surface vtk files
"""
import sys
import os
import os.path as op
import matplotlib.pyplot as plt
from mayavi import mlab
from pipeline.make_networkx import makegraph as mg
from mombud.functions import vtkvizfuncs as vf
import wrappers as wr
import cPickle as pickle
# pylint: disable=C0103
#plt.close('all')
#mlab.close(all=True)
inputdir = op.join(os.getcwd(), 'input')
rawdir = op.join(os.getcwd(), 'output')
direc = op.join(os.getcwd(), 'data')
with open(op.join(direc, 'YPL_grph_b.pkl'), 'rb') as inpt:
    node_data, edge_data, nxgrph =pickle.load(inpt)

nxgrph = nxgrph[78]
# filelist and graph list
if __name__ == '__main__':
    try:
        vtkF = wr.ddwalk(op.join(rawdir, 'normSkel'),
                         '*skeleton.vtk', start=5, stop=-13)
        vtkS = wr.ddwalk(op.join(inputdir, 'surfaceFiles'),
                         '*surface.vtk', stop=-12)

    except Exception:
            print "Error: check your filepaths"
            sys.exit()

    filekeys = {item: vtkF[media][item] for media
                in sorted(vtkF.keys()) for item
                in sorted(vtkF[media].keys())}

    for key in sorted(filekeys.keys())[400:401]:
        key = 'YPL_052315_005_RFPstack_004'
#        data = vf.callreader(vtkF[key[:3]][key])
        data = vf.callreader(vtkF['YPL'][key])
#        node_data, edge_data, nxgrph = mg(data, key)

        figone = mlab.figure(figure=key+"old",
                             size=(1200, 800),
                             bgcolor=(.15, .15, .15))
        vtkobj, _ = vf.cellplot(figone, filekeys[key])
        vf.rendsurf(vtkS[key[:3]][key[4:]])
        vf.labelbpoints(nxgrph, bsize=0.1, esize=0.08)
        vf.labellines(vtkobj)