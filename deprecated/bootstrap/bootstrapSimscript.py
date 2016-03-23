"""
   Bootstrap a simulation of randomly selected branchpoints
"""
import matplotlib.pyplot as plt
import glob
import os
import scipy.stats as sp
import numpy as np
import cPickle as pickle
import seaborn as sns
from bootstrapEnBptsFunc import bootStrp as bstrp
sns.set_context("talk")
sns.set(style="dark")

NUM_RUNS = 5  # number of simulations per bootstrap
plt.close('all')
# =============================================================================
#           Data initialization
# =============================================================================
temp = []
for root, dirs, files in os.walk(os.getcwd()):
    for f in dirs:
        temp.append(
            os.path.join(root, f))

for mem in temp:
    print('Now on %s' % mem+"\n"+"="*79)
    graphPkl = os.path.join(mem, mem[-3:]+'_grph.pkl')
    with open(graphPkl, 'rb') as input:
        nodesEdges = pickle.load(input)
    files = glob.glob(mem+'\Norm*vtk')
    labs = mem[-3:]
    C = len(files)
    nodes = [nodesEdges[0][i] for i in range(C)]
    edges = [nodesEdges[1][i] for i in range(C)]
    Grphs = [nodesEdges[2][i] for i in range(C)]

# =============================================================================
#           Main Function block
# =============================================================================
    temp2 = []
    for runs in range(NUM_RUNS):
        enrchd, difft = bstrp(0.3, files, Grphs)[-2:]
        temp2.append([len(enrchd), len(difft)])

    enrichBP = np.mean([i[0] for i in temp2])
    diffBP = np.mean([i[1] for i in temp2])
    deenrichBP = diffBP-enrichBP

    print("%s mean # of enrich BP : %4.3f\n"
          "%s mean # of dennrich BP : %4.3f\n" %
          (labs, enrichBP, labs, deenrichBP))
# =============================================================================
#           PLOT
# =============================================================================
    fig = plt.figure(figsize=(7, 11))
    kernel = sp.gaussian_kde([i[0] for i in temp2])
    kernel.set_bandwidth(bw_method='silverman')
    kernel2 = sp.gaussian_kde([i[1]-i[0] for i in temp2])
    kernel2.set_bandwidth(bw_method='silverman')
    positions = np.linspace(0, 8, num=50)

    Z = kernel(positions)
    Z2 = kernel2(positions)

    ax1 = plt.subplot(211)
    ax1.fill_between(positions, Z, color='b', alpha=.3)
    ax1.axvline(enrichBP, color='b', alpha=.5)
    ax1.text(
        .3, .5, 'Number of enriched branchpoints: %4.2f'
        % enrichBP,
        transform=ax1.transAxes, verticalalignment='top')

    ax2 = plt.subplot(212)
    ax2.axvline(deenrichBP, color='g', alpha=0.5)
    ax2.fill_between(positions, Z2, color='g', alpha=.3)
    ax2.text(
         .3, .5, 'Number of de-enriched branchpoints: %4.2f'
         % (deenrichBP),
         transform=ax2.transAxes, verticalalignment='top')

    ax1.text(
         .5, .85, labs, transform=ax1.transAxes,
         fontsize=14, verticalalignment='top')
    plt.show()
