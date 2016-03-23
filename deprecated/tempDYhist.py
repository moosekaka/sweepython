"""
    Plot the KDE means channel distributions of cell
    and populations in different carbon sources
"""
import os
import fnmatch
import scipy.stats as sp
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from tvtk.api import tvtk

# =============================================================================
#               init vars and get vtk, graph data
# =============================================================================
G = {}
Bpoints = {}
Graphs = {}
vtkF = {}
BranchPtMeanFP = {}
notBPMeanFP = {}
backgroundGFP = {}
backgroundRFP = {}
bckgrndGFP = {}
bckgrndRFP = {}
Type = {}
Dates = {}
meanFP = {}
parDir = os.path.dirname(os.getcwd())
colors = {
    'YPD': 'b',
    'YPE': 'g',
    'YPR': 'r',
    'YPL': 'm'}
condition = {
    'All Dates': ['042515', '042715', '042915', 'c42515', '052315'],
    '04-27-15': ['042715'],
    '04-25-15': ['042515', 'c42515'],
    '04-29-15': ['042915'],
    'c4-25-15': ['c42515'],
    '05-23-15': ['052315']}
conditionKey = 'All Dates'

with open(parDir + '\\fileMetas.pkl', 'rb') as inpt:
    s = pickle.load(inpt)

for root, dirs, files in os.walk(os.getcwd()):
    for i in files:
        if fnmatch.fnmatch(i, '*grph.pkl'):
            G.setdefault(
                root.rsplit('\\', 1)[1], []).append(os.path.join(root, i))
        if fnmatch.fnmatch(i, '*vtk'):
            vtkF.setdefault(
                root.rsplit('\\', 1)[1], []).append(os.path.join(root, i))
media = sorted(G.keys())

#   get the graphs objects and hash into Graphs
for i in media:
    with open(G[i][0], 'rb') as inpt:
        temp = pickle.load(inpt)[2]
        Graphs.setdefault(i, []).append(temp)

#    iterate thru every type of media, find and hash the bpoints coords ##
for mem in media:
    for n, m in enumerate(Graphs[mem][0]):
        graphName = m.graph['cell']
        Bpoints.setdefault(graphName, [])
        B = m.nodes(data=True)
        Bpoints[graphName] = {
            i[0]: i[1]['coord'] for i in B if i[1]['degree'] > 2}

#    hash the metaDatas for the file/graph cells
with open(parDir+'\\'+'fileMetas.pkl', 'rb') as inpt:
    s = pickle.load(inpt)
for i in s:
    backgroundRFP[i] = s[i][0]
    backgroundGFP[i] = s[i][1]
    Type[i] = i[:4]
    Dates[i] = i.split('_')[1][:]

# =============================================================================
#               Begin Routine
# =============================================================================
for mem in media:
    for n, a in enumerate(vtkF[mem]):
        placeHold = []
        notBP = []
        fileKey = a.rsplit('\\', 1)[1][5:-13]
        reader = tvtk.PolyDataReader()
        reader.set(file_name=a)
        reader.update()
        data = reader.output
        temp = data.point_data
        scalarsNorm = np.ravel(data.point_data.scalars)
        dyRaw = np.ravel(temp.get_array('DY_raw'))
        rawGFP = np.ravel(temp.get_array('rGFP'))
        rawRFP = np.ravel(temp.get_array('rRFP'))
        indices = np.indices(rawRFP.shape)[0]

#   list of branchpoints  in cell
        BP = Bpoints[fileKey]
        BranchPtMeanFP.setdefault(mem, [])
        notBPMeanFP.setdefault(mem, [])

        if Dates[fileKey] in condition[conditionKey]:

            for i in BP:
                polyData = tvtk.PolyData()
                polyData.points = data.points
                result = tvtk.IdList()
                loc = tvtk.PointLocator()
                loc.set(data_set=polyData)
                loc.build_locator()
                loc.find_points_within_radius(0.3, BP[i], result)
                placeHold.append(result)

#   branchpoints

            placeHold = [el for lis in placeHold for el in lis]
            BranchPtMeanFP[mem].append([(scalarsNorm[j]) for j in placeHold])
#   not branchpoints
            notBP = np.setdiff1d(indices, placeHold)
            notBPMeanFP[mem].append([(scalarsNorm[j]) for j in notBP])

# =============================================================================
#             PLOT Population Intensities
# =============================================================================
plt.close('all')
gs = gridspec.GridSpec(4, 1)
fig = plt.figure(figsize=(11, 8.5))
plt.show()
temp = []

for h, i in enumerate(notBPMeanFP):
    temp2 = [el for ls in notBPMeanFP[i] for el in ls]  # Pop NON BRANCHPOINT
    temp4 = [el for ls in BranchPtMeanFP[i] for el in ls]  # BRANCHPOINTS
#   Population BRANCHPOINTS ##########
    ax1 = plt.subplot(gs[h])
    ax2 = plt.subplot(gs[h])
    kernel = sp.gaussian_kde(temp2)
    kernel2 = sp.gaussian_kde(temp4)
    positions = np.linspace(0, 1, num=100)

    Z = kernel(positions)
    ax1.fill_between(positions, Z, color=colors[i], alpha=0.2)

    ax1.text(
        0.5, 0.9, i, transform=ax1.transAxes, fontsize=14,
        verticalalignment='top')

    ax1.axvline(
        np.percentile(temp2, 50), color='blue', lw=1.5, alpha=0.5)

    Z2 = kernel2(positions)
    ax2.fill_between(positions, Z2, color=colors[i], alpha=0.2, hatch="/")

    ax2.text(
        0.5, 0.7, "Non Branchpoints median=%6.3f" % np.percentile(temp2, 50),
        transform=ax1.transAxes, fontsize=12, verticalalignment='top')

    ax2.text(
        0.5, 0.5, "Branchpoints median=%6.3f" % np.percentile(temp4, 50),
        transform=ax1.transAxes, fontsize=12, verticalalignment='top')

    ax2.axvline(np.percentile(temp4, 50), color='red', lw=1.5, alpha=0.5)

plt.suptitle(
    'Scaled $\Delta \Psi$ Intensities Distributions for Branchpoints '
    '(hatched) and non branchpoints ', fontsize=14)
plt.show()

# =============================================================================
#            PLOT Cells Intensities
# =============================================================================

for mem in media:
    N = int(np.ceil(np.sqrt(len(BranchPtMeanFP[mem]))))  # for grid size
    gs = gridspec.GridSpec(N, N)
    fig = plt.figure(figsize=(17, 11))

    temp = []
    for h, i in enumerate(notBPMeanFP[mem]):
        n_bpoints = [el for el in i]  # Non BraNCHPOINTS #
        b_points = [el for el in BranchPtMeanFP[mem][h]]  # BRANCHPOINTS

        #   BRANCHPOINTS
        ax1 = plt.subplot(gs[h])
        ax2 = plt.subplot(gs[h])
        kernel = sp.gaussian_kde(n_bpoints)
        kernel2 = sp.gaussian_kde(b_points)
        positions = np.linspace(0, 1, num=100)

        #   non brpts
        Z = kernel(positions)
        ax1.fill_between(positions, Z, color='green', alpha=0.2)
        ax1.axvline(
            np.percentile(n_bpoints, 50), color='blue', lw=1.5, alpha=0.5)

        ax1.get_xaxis().set_ticklabels([])
        ax1.get_yaxis().set_ticklabels([])
        #   bpts
        Z2 = kernel2(positions)
        ax2.fill_between(positions, Z2, color='magenta', alpha=0.2)
        ax2.axvline(
            np.percentile(b_points, 50), color='red', lw=1.5, alpha=0.5)

    plt.suptitle(
        'Scaled Normed $\Delta \Psi$ Intensities Distributions '
        'for Branchpoints (magenta) and non branchpoints (green) in ' +
        mem, fontsize=14)

plt.show()
