# -*- coding: utf-8 -*-
"""
Plot mitoskel network in with various scalar values
"""
import sys
import os
import os.path as op
import matplotlib.pyplot as plt
from mayavi import mlab
from mayavi.sources.vtk_file_reader import VTKFileReader
from pipeline.make_networkx import makegraph as mg
from mombud.vtk_viz import vtkvizfuncs as vf
import wrappers as wr
# pylint: disable=C0103
plt.close('all')
mlab.close(all=True)
datadir = op.join(os.getcwd(), 'mutants')
inptdir = op.join(os.getcwd(), 'mutants')
outputdir =  op.join(os.getcwd(), 'mutants')


def setup_data(fname):
    """Given a VTK file name `fname`, this creates a mayavi2 reader
    for it and adds it to the pipeline.  It returns the reader
    created.
    """
    # 'mayavi' is always defined on the interpreter.
    d = VTKFileReader()
    d.initialize(fname)
    return d


# filelist and graph list
if __name__ == '__main__':

    filekey = 'WT_032716_007_RFPstack_060'
    try:
        vtkF = wr.swalk(op.join(outputdir, 'normalizedVTK'),
                        'N*Skeleton.vtk', start=5, stop=-13)
        vtkS = wr.swalk(op.join(inptdir, 'surfaceFiles'),
                        '*surface.vtk', stop=-12)

    except Exception:
        print ("Check your filepaths\nSearch directory is %s\n" % inptdir)
        sys.exit()

    source = setup_data(vtkF[filekey])
    data = vf.callreader(vtkF[filekey])
    node_data, edge_data, nxgrph = mg(data, filekey)

    figone = mlab.figure(figure=filekey,
                         size=(1200, 800),
                         bgcolor=(.0, .0, .0))
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
                                      rad=.06)
        vtktube.actor.mapper.scalar_visibility = True  # False for no heatmap
        #    vf.rendsurf(vtkS[filekey[:3]][filekey[4:]])
        vf.labelbpoints(nxgrph, bsize=.12, esize=0.06)
        mlab.savefig(op.join(datadir, 'pipelineFigs', i + '.png'))

    f, ax = plt.subplots()
    vf.nicegrph(nxgrph, ax)
