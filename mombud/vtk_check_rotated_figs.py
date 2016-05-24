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
from mombud.vtk_viz import vtkvizfuncs as vf
import wrappers as wr

# pylint: disable=C0103
plt.close('all')
mlab.close(all=True)
inputdir = op.join(os.getcwd(), 'mutants', 'test')
rawdir = op.join(os.getcwd(), 'mutants', 'test')


# filelist and graph list
if __name__ == '__main__':
    try:
        vtkF = wr.swalk(op.join(rawdir),
                         '*.vtk', start=0, stop=-4)

    except Exception:
            print "Error: check your filepaths"
            sys.exit()

    for key in vtkF.keys():
        temp = key.partition("_")
        etype = temp[0]
        cellkey = temp[-1]
        data = vf.callreader(vtkF[key])
        node_data, edge_data, nxgrph = mg(data, key)
        figone = mlab.figure(figure=key,
                             size=(800, 600),
                             bgcolor=(.15, .15, .15))
        vtkobj, _ = vf.cellplot(figone, vtkF[key])
        mlab.savefig(op.join(rawdir, '%s.png' % key))
        mlab.close()
