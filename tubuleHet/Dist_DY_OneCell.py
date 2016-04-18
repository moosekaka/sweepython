# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 16:12:01 2015

@author: sweel
"""
# pylint: disable=C0103
import glob
import os
import os.path as op
import cPickle as pickle
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from tubuleHet.autoCor.fitDistr import fitDist
import wrappers as wr
import vtk
from tvtk.api import tvtk
from pipeline.make_networkx import makegraph
sns.set_context("talk")
sns.set(style="darkgrid")
colors = ["medium green",
          "greyish blue",
          "yellowy brown",
          "reddish grey"]
labs = {'real': r'$\Delta \Psi$ Actual ',
        'norm': r'$\Delta \Psi$ Normal Dist. Fit',
        'shuf': r'$\Delta \Psi$ Shuffled Dist. Fit',
        'unif': r'$\Delta \Psi$ Uniform Dist. Fit'}

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
#    axeshand.text(0.4, 0.15, labels,
#                  fontsize=24, ha='center', va='center',
#                  transform=axeshand.transAxes)
    axeshand.plot([el for lis in dstr for el in lis], **kwargs)
    _, yl1 = axeshand.get_ybound()
    axeshand.set_yticks(np.arange(500, yl1, 1000))
    plt.setp(axeshand.get_xticklabels(), visible=True)
    count = 0
    for h, i in enumerate(dstr):
        count += len(i)
        axeshand.axvline(count,
                         color='k',
                         lw='1.5',
                         alpha=.6,
                         ls='-.')
    return axeshand

# =============================================================================
#   main block
# =============================================================================
if __name__ == "__main__":
# ONLY RUN THESE LINES IF WANT TO REFIT NEW DISTRIBUTIONS!!
    filekey = "YPE_042515_001_RFPstack_000"
    mediatype = "YPE"

    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(vtkF[mediatype][filekey])
    reader.Update()
    vtkdata = reader.GetOutput()
    tvtkdata = tvtk.to_tvtk(vtkdata)

    _, _, nxgrph = makegraph(vtkdata, filekey)
    output = fitDist(tvtkdata, nxgrph)
    sampN, sampU, Norm, NormPermute = output[4:8]
#    for root, dirs, fls in os.walk(os.getcwd()):
#        for f in dirs:
#            if f.startswith('YPE'):
#                dirlist[f[:3]] = (os.path.join(root, f))
#
#    files = glob.glob(dirlist['YPE']+r'\Norm*vtk')
#    with open(
#        os.path.join(
#            dirlist['YPE'], 'YPE_grph.pkl'), 'rb') as inpt:
#        G = pickle.load(inpt)[2]
#        Graphs = {i.graph['cell']: i for i in G}
#        output = fitDist(files, Graphs)

    data = output[0]
    sampN, sampU, Norm, NormPermute = output[5:9]

    cell = sorted(data.keys())[0]  # YPE_042515_001_RFPstack_000'

    for line in range(data[cell].number_of_lines):
        M = data[cell].get_cell(line).number_of_points
        randNDY.append(sampN[cell].rvs(size=M))
        randUDY.append(sampU[cell].rvs(size=M))
    with open('DataForOneCellDist.pkl', 'wb') as output:
        pickle.dump(
            (cell,randNDY,randUDY,Norm[cell],NormPermute[cell]), output)
# =============================================================================

    with open('DataForOneCellDist.pkl', 'rb') as INPT:
        cell, randNDY, randUDY, Norm, NormPermute = pickle.load(INPT)
    with sns.plotting_context('talk', font_scale=1.5):
        f, axarr = plt.subplots(4, sharex=True, figsize=(11, 8.25))
        plotcelldist(axarr[0], Norm, labs['real'], color=color[0])
        plotcelldist(axarr[1], randNDY, labs['norm'], color=color[1])
        plotcelldist(axarr[2], NormPermute, labs['shuf'], color=color[2])
        plotcelldist(axarr[3], randUDY, labs['unif'], color=color[3])
    plt.show()
