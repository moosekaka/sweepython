# -*- coding: utf-8 -*-
"""
Plot mitoskel network in with various scalar values
"""
import sys
import os
import os.path as op
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sp
from mayavi import mlab
from tvtk.api import tvtk
from pipeline.make_networkx import makegraph as mg
from mombud.vtk_viz import vtkvizfuncs as vf
import wrappers as wr

# pylint: disable=C0103
plt.close('all')
mlab.close(all=True)
datadir = op.join(os.getcwd(), 'data')
inptdir = op.join(os.getcwd(), 'input')


def create_vtk(pt, lns, scl):
    """
    wrapper to make tvtk polydata from given pts, lines , scalars
    """
    src = tvtk.PolyData()
    src.points = pt
    src.lines = lns
    src.point_data.scalars = scl
    return src


def vtkarray(scl, arrname):
    """
    wrapper to create array and name it
    """
    temp = tvtk.DoubleArray()
    temp.from_array(scl)
    temp.set(name=arrname)
    return temp

# filelist and graph list
if __name__ == '__main__':

    filekey = 'YPE_042515_001_RFPstack_000'
    try:
        vtkF = wr.swalk(op.join(inptdir, 'tubule'),
                        'N*Skeleton.vtk', start=5, stop=-13)
        vtkS = wr.swalk(op.join(inptdir, 'surfaceFiles'),
                        '*surface.vtk', stop=-12)

    except Exception:
        print "Check your filepaths\nSearch directory is %s\n" % inptdir
        sys.exit()

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
                                      rad=.08)
        vtktube.actor.mapper.scalar_visibility = True  # False for no heatmap
        vf.labellines(vtkobj)
        mlab.savefig(op.join(datadir, 'tubule', i + '.png'))

    # Edgeplots for real scalars
    data = tvtk.to_tvtk(vtkobj.outputs[0])
    pts = data.points.to_array()
    scals = data.point_data.get_array('DY_raw')

    fig2 = mlab.figure(figure=filekey + ' real',
                       size=(1200, 800),
                       bgcolor=(.086, .086, .086))
    vf.edgeplot(fig2, vtkobj.outputs[0], 14)

    # Edgeplot for shuffled scalars
    shufscals = np.random.permutation(scals)
    shufarr = vtkarray(shufscals, 'DY_raw')
    pshuf = create_vtk(pts, data.lines, shufarr)
    fig3 = mlab.figure(figure=filekey + ' shuf',
                       size=(1200, 800),
                       bgcolor=(.0, .0, .0))
    vf.edgeplot(fig3, pshuf, 14, scalartype='DY_raw')

    # Edgeplot for norm dist scalars
    nrmscals = sp.norm(np.mean(scals), np.std(scals))
    nrmarr = vtkarray(nrmscals.rvs(size=len(data.points)), 'DY_raw')
    pnrm = create_vtk(pts, data.lines, nrmarr)
    fig4 = mlab.figure(figure=filekey + ' norm',
                       size=(1200, 800),
                       bgcolor=(.086, .086, .086))
    vf.edgeplot(fig4, pnrm, 14, scalartype='DY_raw')

    # Edgeplot for norm dist scalars
    mn = np.mean(scals)
    dv = np.std(scals)
    a = mn - 1.5 * dv
    urfscals = sp.uniform(a.clip(min=0), 1.5 * mn + dv)
    nrmarr = vtkarray(nrmscals.rvs(size=len(data.points)), 'DY_raw')
    purf = create_vtk(pts, data.lines, nrmarr)
    fig5 = mlab.figure(figure=filekey + ' uniform',
                       size=(1200, 800),
                       bgcolor=(.086, .086, .086))
    vf.edgeplot(fig5, purf, 14, scalartype='DY_raw')
