# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 17:44:48 2015
Script to create a distribution of lags for Intensities along edges of cells
@author: sweel
"""
import glob
import os
import os.path as op
import cPickle as pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from tubuleHet.autoCor.fitDistr import fitDist as fitd
datadir = op.join(os.getcwd(), 'data')
sns.set_context("talk")
sns.set(style="dark")

# =============================================================================
#           Data initialization
# =============================================================================
plt.close('all')
temp = []
for root, dirs, files in os.walk(os.getcwd()):
    for f in dirs:
        if f.startswith('YP'):
            temp.append(
                os.path.join(root, f))

for mem in temp:
    print'Now on %s' % mem+"\n"+"="*79
    randNDY = {}  # normal dist fit
    randUDY = {}  # uniform dist fit
    files = glob.glob(mem+r'\Norm*vtk')
    labs = mem[-3:]

    with open(os.path.join(mem,
                           '%s_grph.pkl' % labs), 'rb') as inpt:
        G = pickle.load(inpt)[2]
        Graphs = {i.graph['cell']: i for i in G}

    output = fitd(files, Graphs)
    data = output[0]
#    1:5 for scaled, 5:9 for unscaled
    lineId = output[9]
#
    sampN, sampU, Norm, NormPermute = output[1:5]
    for cell in data.keys():
        for line in range(data[cell].number_of_lines):
            M = data[cell].get_cell(line).number_of_points
            randNDY.setdefault(cell, []).append(sampN[cell].rvs(size=M))
            randUDY.setdefault(cell, []).append(sampU[cell].rvs(size=M))

    out = (randNDY, randUDY, Norm, NormPermute, data, lineId)
    with open('%s_lagsUnscaled.pkl' % labs, 'wb') as OUT:
        pickle.dump(out, OUT)

#    1:5 for scaled, 5:9 for unscaled
#    sampN, sampU, Norm, NormPermute = output[1:5]
#    for cell in data.keys():
#        for line in range(data[cell].number_of_lines):
#            M = data[cell].get_cell(line).number_of_points
#            randNDY.setdefault(cell, []).append(sampN[cell].rvs(size=M))
#            randUDY.setdefault(cell, []).append(sampU[cell].rvs(size=M))
#    out2 = (randNDY, randUDY, Norm, NormPermute)
#    with open('%s_lagscaled.pkl' % labs, 'wb') as OUT2:
#        pickle.dump(out2, OUT2)

Y=np.hstack(Norm[cell])
X=np.arange(len(Y))


#    plt.plot(freqs[idx], ps[idx])
