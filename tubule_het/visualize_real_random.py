# -*- coding: utf-8 -*-
"""
Plot mitoskel network in with various scalar values
"""
import sys
import os
import os.path as op
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as sp
import seaborn as sns
from mayavi import mlab
from tvtk.api import tvtk
from pipeline.make_networkx import makegraph as mg
from mombud.functions import vtkvizfuncs as vf
import wrappers as wr

# pylint: disable=C0103
plt.close('all')
mlab.close(all=True)
datadir = op.join(os.getcwd(), 'data')
inptdir = op.join(os.getcwd(), 'input')
savefolder = op.expanduser(os.sep.join(['~', 'Dropbox', 'SusanneSweeShared',
                                        'paper']))

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
                         bgcolor=(.05, .05, .05))

    vtkobj, vtktube = vf.cellplot(figone,
                                  vtkF[filekey],
                                  scalartype='DY_raw',
                                  rad=.08)
    vtktube.actor.mapper.scalar_visibility = True  # False for no heatmap
    vf.labellines(vtkobj)

    # Edgeplots for real scalars
    data = tvtk.to_tvtk(vtkobj.outputs[0])
    pts = data.points.to_array()
    scals = data.point_data.get_array('DY_raw')

    fig2 = mlab.figure(figure=filekey + ' real',
                       size=(1200, 800),
                       bgcolor=(.05, .05, .05))
    real, mmr = vf.edgeplot(fig2, vtkobj.outputs[0], 14, scalartype='DY_raw')

    # Edgeplot for shuffled scalars
    shufscals = np.random.permutation(scals)
    shufarr = vtkarray(shufscals, 'DY_raw')
    pshuf = create_vtk(pts, data.lines, shufarr)
    fig3 = mlab.figure(figure=filekey + ' shuf',
                       size=(1200, 800),
                       bgcolor=(.05, .05, .05))
    shuf, mms = vf.edgeplot(fig3, pshuf, 14, scalartype='DY_raw')

    # Edgeplot for norm dist scalars
    nrmscals = sp.norm(np.mean(scals), np.std(scals))
    nrmarr = vtkarray(nrmscals.rvs(size=len(data.points)), 'DY_raw')
    pnrm = create_vtk(pts, data.lines, nrmarr)
    fig4 = mlab.figure(figure=filekey + ' norm',
                       size=(1200, 800),
                       bgcolor=(.05, .05, .05))
    norm, mmn = vf.edgeplot(fig4, pnrm, 14, scalartype='DY_raw')

    # Edgeplot for norm dist scalars
    mn = np.mean(scals)
    dv = np.std(scals)
    a = mn - 1.5 * dv
    urfscals = sp.uniform(a.clip(min=0), 1.5 * mn + dv)
    nrmarr = vtkarray(nrmscals.rvs(size=len(data.points)), 'DY_raw')
    purf = create_vtk(pts, data.lines, nrmarr)
    fig5 = mlab.figure(figure=filekey + ' uniform',
                       size=(1200, 800),
                       bgcolor=(.05, .05, .05))
    unif, mmu = vf.edgeplot(fig5, purf, 14, scalartype='DY_raw')

    # dataframe long form of various scalars in simulation for ΔΨ
    df = pd.DataFrame(data={'shuffled': shuf,
                            'real': real,
                            'normal': norm,
                            'uniform': unif})
    df = df.reset_index()
    df2 = pd.melt(df, id_vars='index', var_name='type')

    # reset scalebars to be the same for edgeplots
    for i in [mmr, mms, mmn, mmu]:
        i.lut.set(range=[1000, 4000])
#==============================================================================
# Plots
#==============================================================================
    with sns.plotting_context('talk', 1.2):
        sns.set_style('whitegrid')
        g = sns.FacetGrid(df2, col_wrap=2, col='type', hue='type',
                          col_order=['real', 'shuffled', 'normal', 'uniform'])
        g.map(plt.plot, 'index', 'value')
        g.set(ylim=[0, 4000],
              yticks=np.arange(0, 4500, 1000),
              ylabel='', xlabel='')
        plt.savefig(op.join(savefolder,
                            'panel_simul.svg'))
        plt.savefig(op.join(savefolder,
                            'panel_simul.png'))

    dic = {'real': fig2,
           'shuf': fig3,
           'norm': fig4,
           'unif': fig5}

    for key, val in dic.iteritems():
        mlab.savefig(op.join(savefolder,
                             ('{}.png').format(key)),
                     figure=val)
