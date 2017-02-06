# -*- coding: utf-8 -*-
"""
Plot mitoskel network in with various scalar values
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
datadir = op.join(os.getcwd(), 'data')
inptdir = op.join(os.getcwd(), 'input')

# filelist and graph list
if __name__ == '__main__':

    filekey = 'YPE_042715_018_RFPstack_052'
    try:
        vtkF = wr.swalk(op.join(inptdir, 'pipelineFigs'),
                        'N*Skeleton.vtk', start=5, stop=-13)
        vtkS = wr.swalk(op.join(inptdir, 'surfaceFiles'),
                        '*surface.vtk', stop=-12)

    except Exception:
        print ("Check your filepaths\nSearch directory is %s\n" % inptdir)
        sys.exit()

    data = vf.callreader(vtkF[filekey])
    node_data, edge_data, nxgrph = mg(data, filekey)

    figone = mlab.figure(figure=filekey,
                         size=(1200, 800),
                         bgcolor=(.086, .086, .086))
    dic = {'DY_minmax',
           'WidthEq',
           'DY_raw',
           'rRFP',
           'rGFP',
           'bkstRFP',
           'bkstGFP'}
    for i in dic:
        vtkobj, vtktube = vf.cellplot(figone,
                                      vtkF[filekey],
                                      scalartype=i,
                                      rad=.08)
        vtktube.actor.mapper.scalar_visibility = True  # False for no heatmap
        #    vf.rendsurf(vtkS[filekey[:3]][filekey[4:]])
        vf.labelbpoints(nxgrph, esize=.12)
        mlab.savefig(op.join(datadir, 'pipelineFigs', i + '.png'))
