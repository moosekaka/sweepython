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
from mombud.functions import vtkvizfuncs as vz
import wrappers as wr
# pylint: disable=C0103
plt.close('all')
mlab.close(all=True)
datadir = op.join(os.getcwd(), 'mutants')
inptdir = op.join(os.getcwd(), 'mutants')
outputdir = op.join(os.getcwd(), 'mutants')

# filelist and graph list
if __name__ == '__main__':

    filekey = 'NUM1_032016_011_RFPstack_030'
    try:
 #        vtkF = wr.swalk(op.join(outputdir, 'normalizedVTK'),
 #                        'N*Skeleton.vtk', start=5, stop=-13)
        vtkF = wr.swalk(os.getcwd(),
                        '*.vtk', start=0, stop=-4)
        vtkS = wr.swalk(op.join(inptdir, 'surfaceFiles'),
                        '*surface.vtk', stop=-4)

    except Exception:
        print "Check your filepaths\nSearch directory is %s\n" % inptdir
        sys.exit()

    data = vz.callreader(vtkF[filekey])
    node_data, edge_data, nxgrph = mg(data, filekey)

    figone = mlab.figure(figure=filekey,
                         size=(1200, 800),
                         bgcolor=(.0, .0, .0))
    scaltypes = ('DY_raw', 'DY_minmax', 'WidthEq', 'rRFP',
                 'rGFP', 'bkstRFP', 'bkstGFP')
    for i in scaltypes[:1]:
        mlab.clf(figure=figone)
        vtkobj, vtktube = vz.cellplot(figone,
                                      vtkF[filekey],
                                      scalartype=i,
                                      scalar_visible=True,
                                      rad=0.08,
                                      legend=True)
        vz.labelbpoints(nxgrph,
                        bsize=.12, esize=0.06,
                        bcol=vz.rgbcol('light magenta'))
        vact = vz.rendsurf(vtkS['%s_surface' % filekey],
                    color=vz.rgbcol('yellow orange'))
#        vact.actor.actor.user_tranform =
    #    mlab.savefig(op.join(datadir, 'pipelineFigs', i + '.png'))

    f, ax = plt.subplots()
    vz.nicegrph(nxgrph, ax)
