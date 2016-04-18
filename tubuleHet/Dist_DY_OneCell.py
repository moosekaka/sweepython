# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 16:12:01 2015

@author: sweel
"""
# pylint: disable=C0103
import os
import os.path as op
import cPickle as pickle
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import vtk
from tvtk.api import tvtk
from tubuleHet.autoCor.fitDistr import fitDist
import wrappers as wr
from pipeline.make_networkx import makegraph
# pylint: disable=C0103

sns.set_context("talk")
sns.set(style="darkgrid")
colors = ["medium green",
          "greyish blue",
          "yellowy brown",
          "reddish grey"]
labs = {'actual': r'$\Delta \Psi$ Actual ',
        'normal': r'$\Delta \Psi$ Normal Dist. Fit',
        'permute': r'$\Delta \Psi$ Shuffled Dist. Fit',
        'uniform': r'$\Delta \Psi$ Uniform Dist. Fit'}

sns.set_palette(sns.xkcd_palette(colors))
color = sns.color_palette()
plt.close('all')

# =============================================================================
#           Data initialization
# =============================================================================
datadir = op.join(os.getcwd(), 'data')
rawdir = op.join(os.getcwd(), 'output')
vtkF = wr.ddwalk(op.join(rawdir, 'normalizedVTK'),
                 '*skeleton.vtk', start=5, stop=-13)


def plotcelldist(axeshand, dstr, labels, **kwargs):
    """Plot one cell with diff distributions

    Parameters
    ----------
    axeshand:
        handle to the axes figure
    dstr:
        the distribution to be plotted
    labels:
        value from a dict of labels to identify distsr

    """
    axeshand.text(0.4, 0.15, labels,
                  fontsize=24, ha='center', va='center',
                  transform=axeshand.transAxes)
    axeshand.plot(np.hstack(dstr), **kwargs)
    _, yl1 = axeshand.get_ybound()
    axeshand.set_yticks(np.arange(500, yl1, 1000))
    plt.setp(axeshand.get_xticklabels(), visible=True)
    count = 0
    for _, i in enumerate(dstr):
        count += len(i)
        axeshand.axvline(count,
                         color='k',
                         lw='1.5',
                         alpha=.6,
                         ls='-.')
    return axeshand

#   main block
if __name__ == "__main__":
    filekey = "YPE_042515_001_RFPstack_000"
    mediatype = "YPE"
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(vtkF[mediatype][filekey])
    reader.Update()
    vtkdata = reader.GetOutput()
    tvtkdata = tvtk.to_tvtk(vtkdata)
    _, _, nxgrph = makegraph(vtkdata, filekey)
# ONLY RUN THESE LINES IF WANT TO REFIT NEW DISTRIBUTIONS!!
#    output = fitDist(tvtkdata, nxgrph)
#    normal, uniform, actual, permute = output[4:8]
#    cell = tvtkdata  # YPE_042515_001_RFPstack_000'
#    with open(op.join(rawdir, 'DataForOneCellDist.pkl'), 'wb') as output:
#        pickle.dump((actual, permute, normal, uniform), output)

    with open(op.join(rawdir, 'DataForOneCellDist.pkl'), 'rb') as INPT:
        actual, permute, normal, uniform = pickle.load(INPT)
    with sns.plotting_context('talk', font_scale=1.5):
        f, axarr = plt.subplots(4, sharex=True, figsize=(11, 8.25))
        plotcelldist(axarr[0], actual, labs['actual'], color=color[0])
        plotcelldist(axarr[1], normal, labs['normal'], color=color[1])
        plotcelldist(axarr[2], permute, labs['permute'], color=color[2])
        plotcelldist(axarr[3], uniform, labs['uniform'], color=color[3])
    plt.show()
