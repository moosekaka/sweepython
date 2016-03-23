import networkx as nx
import math
import vtk
import glob
import string
import os
import scipy.stats as sp
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
from itertools import cycle
import pandas as pd
import matplotlib.gridspec as gridspec

r = pd.read_csv('BackgroundRFPGFP.csv')
g = pd.read_csv('BackgroundRFPGFP.csv')
backgroundGFP = {i[0]: i[1] for i in zip(g['Label'], g['MeanGFP'])}
backgroundRFP = {i[0]: i[1] for i in zip(r['Label'], r['MeanRFP'])}
with open('fileDates.pkl', 'rb') as inpt:
    s = pickle.load(inpt)

plt.close('all')
b = glob.glob(os.getcwd() + '\[0-9][0-9]Y*')  # directory names
labs = [i.rsplit('\\', 1)[1].split('_', 1)[0][2:] for i in b]

minmaxRFP = [[] for k in range(len(b))]
minmaxGFP = [[] for k in range(len(b))]
meanRFP = [[] for k in range(len(b))]
meanGFP = [[] for k in range(len(b))]
dates = [[] for k in range(len(b))]
bckgrndGFP = [[] for k in range(len(b))]
bckgrndRFP = [[] for k in range(len(b))]
colors = ['b', 'g', 'r', 'm', 'c']
cycler = cycle(colors)

for k in range(len(b)):

    fileLoc = b[k] + '\\2.5_raw*vtk'
    files = glob.glob(fileLoc)

    for a in range(len(files)):
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(files[a])
        reader.Update()
        data = reader.GetOutput()
        rawGFP = data.GetPointData().GetScalars("rGFP")
        rawRFP = data.GetPointData().GetScalars("rRFP")
############ the min, max and mean values of  the skeletons are here #####
        RFPtemp = [
            rawRFP.GetTuple(i) for i in range(
                rawRFP.GetNumberOfTuples())]
        GFPtemp = [
            rawGFP.GetTuple(i) for i in range(
                rawGFP.GetNumberOfTuples())]
        meanRFP[k].append(np.mean(RFPtemp))
        meanGFP[k].append(np.mean(GFPtemp))

        fileKey = files[a].rsplit('\\', 1)[1][8:][:-13]
        minmaxRFP[k].append(rawRFP.GetRange())
        minmaxGFP[k].append(rawGFP.GetRange())
        if(backgroundRFP[fileKey] > minmaxRFP[k][a][0]):
            minA = backgroundRFP[fileKey] - 55
        else:
            minA = backgroundRFP[fileKey]
        minB = min(backgroundGFP[fileKey], minmaxGFP[k][a][0])
        bckgrndGFP[k].append(minB)
        bckgrndRFP[k].append(minA)
        dates[k].append(s[fileKey])

##################################   PLOT  GFP   #########################
gs = gridspec.GridSpec(2, 3)
fig = plt.figure(figsize=(17, 11))
fig.text(0.5, 0.95, 'GFP intensities', fontsize=22, ha='center')

for k in range(len(b)):
    clr = cycler.next()
    ax = plt.subplot(gs[k])
    ax.set_ylim(1000, 14000), ax.set_xlim(0, 115)
    data = minmaxGFP[k]
    data2 = meanGFP[k]
    dateX = dates[k]
    ax.plot([i for i in data2], 'o', color=clr, ms=5)
    ax.plot([i[0] for i in data], '_', mew=2, ms=5, color=clr, alpha=.8)
    ax.plot([i[1] for i in data], '_', mew=2, ms=5, color=clr)
    ax.plot([i for i in bckgrndGFP[k]], 'v', mew=.5, ms=5, color=clr, alpha=.8)
    ax.vlines(
        range(
            len(data)), [
            i[0] for i in data], [
                i[1] for i in data], colors=clr, linewidth=.5)
    ###### mark out different dates with and 'x' #############################
    ax.plot([h for h,
             i in enumerate(data2) if dateX[h] == '062914'],
            [i for h,
             i in enumerate(data2) if dateX[h] == '062914'],
            'x',
            color='y',
            mew=1,
            ms=5)
    ##########################################################################
    ax.text(
        0.65,
        0.9,
        labs[k],
        transform=ax.transAxes,
        fontsize=14,
        va='top',
        ha='center')
plt.show()

###################################### PLOT RFP ##########################
gs2 = gridspec.GridSpec(2, 3)
fig2 = plt.figure(figsize=(17, 11))
fig2.text(0.5, 0.95, 'RFP intensities', fontsize=22, ha='center')

for k in range(len(b)):
    clr = cycler.next()
    ax2 = plt.subplot(gs2[k])
    ax2.set_ylim(1000, 14000), ax2.set_xlim(0, 115)
    data = minmaxRFP[k]
    data2 = meanRFP[k]
    dateX = dates[k]
    ax2.plot([i for i in data2], 'o', color=clr, ms=5)
    ax2.plot([i[0] for i in data], '_', mew=2, ms=5, color=clr)
    ax2.plot([i[1] for i in data], '_', mew=2, ms=5, color=clr)
    ax2.plot([i for i in bckgrndRFP[k]],
             'v',
             mew=.5,
             ms=5,
             color=clr,
             alpha=.8)
    ax2.vlines(
        range(
            len(data)), [
            i[0] for i in data], [
                i[1] for i in data], colors=clr, linewidth=.5)
    ###### mark out different dates with and 'x' #############################
    ax2.plot([h for h,
              i in enumerate(data2) if dateX[h] == '062914'],
             [i for h,
              i in enumerate(data2) if dateX[h] == '062914'],
             'x',
             color='y',
             mew=1,
             ms=5)
    ##########################################################################
    ax2.text(
        0.65,
        0.9,
        labs[k],
        transform=ax2.transAxes,
        fontsize=14,
        va='top',
        ha='center')
plt.show()

###################################### PLOT GFP vs RFP ###################

gs3 = gridspec.GridSpec(2, 3)
fig3 = plt.figure(figsize=(17, 11))
fig3.text(0.5, 0.95, 'GFP vs RFP intensities', fontsize=22, ha='center')
ax4 = plt.subplot(gs3[-1])
ax4.text(
    0.5,
    0.9,
    'All Groups',
    transform=ax4.transAxes,
    fontsize=14,
    va='top',
    ha='center')

for k in range(len(b)):
    clr = cycler.next()
    ax3 = plt.subplot(gs3[k])
    ax3.set_ylim(2000, 8500), ax3.set_xlim(2000, 8500)
    dataG = [i for i in meanGFP[k]]
    dataR = [i for i in meanRFP[k]]
    dateX = dates[k]
    ax3.plot(dataR, dataG, 's', color=clr, ms=5, alpha=0.85)
    ax4.plot(dataR, dataG, 's', color=clr, ms=5, alpha=0.75)
    ####### mark out different dates with and 'x' ############################
    ax3.plot([i for h,
              i in enumerate(dataR) if dateX[h] == '062914'],
             [i for h,
              i in enumerate(dataG) if dateX[h] == '062914'],
             'x',
             color='y',
             mew=1,
             ms=5)
    ax4.plot([i for h,
              i in enumerate(dataR) if dateX[h] == '062914'],
             [i for h,
              i in enumerate(dataG) if dateX[h] == '062914'],
             'x',
             color='y',
             mew=1,
             ms=5)
    ##########################################################################
    ax3.text(
        0.5,
        0.9,
        labs[k],
        transform=ax3.transAxes,
        fontsize=14,
        va='top',
        ha='center')
plt.show()
