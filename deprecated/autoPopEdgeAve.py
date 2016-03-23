"""
Created on Wed Jun 17 23:43:41 2015
population autocor coef averaged by pixel length
@author: sweel
"""
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import islice, izip_longest

minL = 10
maxL = 45
celltemp = {}
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
fig = plt.figure("autoCors", figsize=(11, 8.25))

for d, f in enumerate(media):
    ax = plt.subplot(gs[d])
    plt.ylim(-.5, 1.)
    ax.text(7.5, .5, f)
    ax.axhline(color='k', alpha=0.5)

    for a, b in enumerate(range(minL, maxL, 5)):
        X = np.ravel(
            [el for lis in s[f] for el in lis
             if len(el) > b])
        celltemp.setdefault(f, [])
        temp = []
        for i in islice(izip_longest(*X, fillvalue=.0), 0, 15):
            temp.append(np.mean(i))
        celltemp[f].append(temp)

# =============================================================================
#        plot
# =============================================================================
        plt.plot(celltemp[f][a],
                 alpha=.35 * np.log(a+2), color=colors[f])

plt.suptitle("Pop. AutoCor Coef for various minimum "
             "edge lengths of %4.2f - %4.2f $\mu$m in "
             " increments of %4.3f $\mu$m" %
             (minL*.055, maxL*.055, .055), fontsize=13)
