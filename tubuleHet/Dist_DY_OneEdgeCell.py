# -*- coding: utf-8 -*-
"""
Created on Sun Jul 05 22:14:07 2015
plot ONE EDGE of ONE CELL, run the commented lines below main block to
generate the distributions first
@author: sweel
"""
# pylint: disable=C0103
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import cPickle as pickle
sns.set_context("talk")

sns.set(rc={"legend.markerscale": 3})
colors = ["medium green",
          "greyish blue",
          "yellowy brown",
          "reddish grey"]
labs = {'real': r'$\Delta \Psi$ Actual ',
        'norm': r'$\Delta \Psi$ Normal Dist Fit',
        'shuf': r'$\Delta \Psi$ Shuffled Dist Fit',
        'unif': r'$\Delta \Psi$ Uniform Dist Fit'}

sns.set_palette(sns.xkcd_palette(colors))
color = sns.color_palette()
#plt.close('all')


def plotedgedistr(axeshand, dstr, labels, lineid, length=None, **kwargs):
    """Plot edge segments of a distribution

    Parameters
    ----------
        axeshand:
            handle to the axes figure
        dstr:
            the distribution to be plotted
        labels:
            value from a dict of labels to identify distsr
        lineid:
            cellid of line to be plotted
        length:
            length of line to be plotted,default is full length of line

    """
    if length is None:
        length = len(dstr[lineid])
    A = np.hstack(dstr)
    y = dstr[lineid]
    bounds = np.argwhere(A==y[0])
#    axeshand.text(0.25, 0.15, labels,
#                  fontsize=24, ha='center', va='center',
#                  transform=axeshand.transAxes)
    x=np.arange(250, 294)
#    axeshand.plot([el for el in dstr[lineid][0:length]], **kwargs)
    axeshand.plot(x,y,**kwargs)
    _, yl1 = axeshand.get_ybound()
    axeshand.set_yticks(np.arange(500, yl1, 1000))
    plt.setp(axeshand.get_xticklabels(), visible=True)
    return axeshand

# =============================================================================
#   Main block
# =============================================================================
if __name__ == "__main__":
    with open('DataForOneCellDist.pkl', 'rb') as INPT:
        cell, randNDY, randUDY, Norm, NormPermute = pickle.load(INPT)
    sns.set(style="whitegrid")
    with sns.plotting_context('talk', font_scale=1.5):

        fig, axarr = plt.subplots(4, sharex=True, figsize=(11, 8.25))
        plotedgedistr(axarr[0], Norm, labs['real'],
                      16, color=color[0])
        plotedgedistr(axarr[1], randNDY, labs['norm'],
                      16, color=color[1])
        plotedgedistr(axarr[2], NormPermute, labs['shuf'],
                      16, color=color[2])
        plotedgedistr(axarr[3], randUDY, labs['unif'],
                      16, color=color[3])
    plt.show()

# =============================================================================
# ONLY RUN THESE LINES IF WANT TO REFIT NEW DISTRIBUTIONS!!
#    dirlist = {}
#    files = []
#    randNDY = []
#    randUDY = []
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
#        output = fitd(files, Graphs)
#
#    data = output[0]
#    sampN, sampU, Norm, NormPermute = output[5:9]
#
#    cell = sorted(data.keys())[0]  # YPE_042515_001_RFPstack_000'
#
#    for line in range(data[cell].number_of_lines):
#        M = data[cell].get_cell(line).number_of_points
#        randNDY.append(sampN[cell].rvs(size=M))
#        randUDY.append(sampU[cell].rvs(size=M))
#    with open('DataForOneCellDist.pkl', 'wb') as output:
#        pickle.dump(
#            (cell,randNDY,randUDY,Norm[cell],NormPermute[cell]), output)
#    plotHet(cell, randNDY, randUDY, Norm[cell], NormPermute[cell])
# =============================================================================
