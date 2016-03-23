"""
Created on Wed Jun 17 21:41:36 2015
cell by cell autocor coef averaged by pixel length along edge
@author: sweel
"""
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import islice, izip_longest

colors = {
    'YPD': 'b',
    'YPE': 'g',
    'YPR': 'r',
    'YPL': 'm'}

with open('autoCors.pkl', 'rb') as inpt:
    s = pickle.load(inpt)

media = sorted(s.keys())

plt.close('all')
gs = gridspec.GridSpec(2, 2)
fig = plt.figure(figsize=(11, 8.25))

for d, f in enumerate(media):
    celltemp = {}
    ax = plt.subplot(gs[d])
    for g, h in enumerate(s[f]):
        temp = []
        for i in islice(
                        izip_longest(*h, fillvalue=.0), 0, 15):
            temp.append(np.mean(i))
        celltemp[g] = temp

#    Plot
    plt.ylim(-.8, 1.)
    for i in range(len(celltemp)):
        plt.plot(celltemp[i],
                 alpha=.25*(i % 2)+.25, color=colors[f])
    ax.text(5, .5, f)
    ax.axhline(color='k', alpha=0.7)


plt.suptitle("Cell. AutoCor Coef for all edge lengths",
             fontsize=13)
s = None
celltemp = None
