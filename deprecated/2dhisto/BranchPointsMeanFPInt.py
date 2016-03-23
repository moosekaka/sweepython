import os
import fnmatch
import scipy.stats as sp
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from tvtk.api import tvtk
import seaborn as sns
import pandas as pd

sns.set_style("dark")
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

colors = {'YPD': 'b', 'YPE': 'g', 'YPR': 'r', 'YPL': 'm'}

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

parDir = os.path.dirname(os.getcwd())

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
# =============================================================================
#  iterate thru every type of media, find and hash the bpoints coords
# =============================================================================
for mem in media:
    for n, m in enumerate(Graphs[mem][0]):
        graphName = m.graph['cell']
        Bpoints.setdefault(graphName, [])
        B = m.nodes(data=True)
        Bpoints[graphName] = {
            i[0]: i[1]['coord'] for i in B if i[1]['degree'] > 2}

#   hash the metaDatas for the file/graph cells
with open(parDir + '\\'+'fileMetas.pkl', 'rb') as inpt:
    s = pickle.load(inpt)
for i in s:
    backgroundRFP[i] = s[i][0]
    backgroundGFP[i] = s[i][1]
    Type[i] = i[:4]
    Dates[i] = i.split('_')[1][:]
plt.close('all')

# =============================================================================
#       Main part of the routine
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
        rawGFP = np.ravel(temp.get_array('rGFP'))
        rawRFP = np.ravel(temp.get_array('rRFP'))
        indices = np.indices(rawRFP.shape)[0]

#   Establish Background levels of channel
        if(backgroundRFP[fileKey] > min(rawRFP)):
            minA = backgroundRFP[fileKey]-1
        else:
            minA = backgroundRFP[fileKey]-1
            minB = min(backgroundGFP[fileKey], min(rawGFP))
        bckgrndGFP[fileKey] = minB
        bckgrndRFP[fileKey] = minA

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
                loc.find_points_within_radius(.3, BP[i], result)
                placeHold.append(result)

#            without background subtraction
#               BranchPtMeanFP[mem].append(
#                   [(rawRFP[j],rawGFP[j]) for j in result])
#            with background subtraction
                BranchPtMeanFP[mem].append(
                    [(rawRFP[j]-bckgrndRFP[fileKey], rawGFP[j] -
                     bckgrndGFP[fileKey]) for j in result])

            placeHold = [el for lis in placeHold for el in lis]

            notBP = np.setdiff1d(indices, placeHold)  # indices to find nonbp
            notBPMeanFP[mem].append(
                [(rawRFP[j]-bckgrndRFP[fileKey], rawGFP[j] -
                 bckgrndGFP[fileKey]) for j in notBP])


# =============================================================================
#    Gaussian KDEs  GFP-RFP for BRANCHPOINTS
# =============================================================================
gs2 = gridspec.GridSpec(2, 2)
fig2 = plt.figure(figsize=(17, 11))
xmin = 0
ymin = 0
xmax = 5500
ymax = 5500

for h, ind in enumerate(media):
    if len(BranchPtMeanFP[ind]):
        temp = [el for lis in BranchPtMeanFP[ind] for el in lis]
        y = [j[1] for j in temp]
        x = [j[0] for j in temp]

        X, Y = np.mgrid[xmin:xmax:55j, ymin:ymax:55j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        kernel = sp.gaussian_kde([x, y])
        Z = np.reshape(kernel(positions).T, X.shape)

        axes = plt.subplot(gs2[h])
        axes.imshow(
            np.rot90(Z), cmap=plt.cm.coolwarm,
            extent=[xmin, xmax, ymin, ymax])
        axes.set_title(ind)
        axes.set_aspect('equal')

plt.suptitle(
    '2d GFP-RFP Pixel Intensities Distributions for Branchpoints on ' +
    conditionKey)
plt.show()

# =============================================================================
#        boxplot
# =============================================================================
data = []
colors2 = []
plt.figure(figsize=(11, 8.25))
for h, ind in enumerate(media):
    temp = [el for ls in BranchPtMeanFP[ind] for el in ls]
    data.append([j[1] for j in temp])
    colors2.append(colors[ind])
bplot = plt.boxplot(data, labels=media, notch=1, patch_artist=True)

plt.ylim(0, 5500)
for patch, color in zip(bplot['boxes'], colors2[:]):
    patch.set_facecolor(color)
    patch.set_alpha(.7)
    patch.set_linewidth(1)
plt.suptitle('Branchpoints GFP intensities on ' + conditionKey)
plt.show()

# =============================================================================
#       Gaussian KDEs  GFP-RFP for NON BPTS
# =============================================================================
gs2 = gridspec.GridSpec(2, 2)
fig2 = plt.figure(figsize=(17, 11))
xmin = 0
ymin = 0
xmax = 5500
ymax = 5500

for h, ind in enumerate(media):
    if len(notBPMeanFP[ind]):
        temp = [el for lis in notBPMeanFP[ind] for el in lis]
        y = [j[1] for j in temp]
        x = [j[0] for j in temp]
        x1 = pd.Series(x, name="RFP")
        x2 = pd.Series(y, name="GFP")

#        X, Y = np.mgrid[xmin:xmax:55j, ymin:ymax:55j]
#        positions = np.vstack([X.ravel(), Y.ravel()])
#        kernel = sp.gaussian_kde([x, y])
#        Z = np.reshape(kernel(positions).T, X.shape)

        axes = plt.subplot(gs2[h])
        axes.imshow(
            np.rot90(Z), cmap=plt.cm.coolwarm,
            extent=[xmin, xmax, ymin, ymax])
#        sns.jointplot(x1, x2, kind="kde", size=7, space=0)
        axes.set_title(ind)
        axes.set_aspect('equal')

plt.suptitle(
    '2d GFP-RFP Pixel Intensities Distributions for non Branchpoints on ' +
    conditionKey)
plt.show()

# =============================================================================
#        boxplot
# =============================================================================
data = []
colors2 = []
plt.figure(figsize=(11, 8.25))
plt.show()
for h, ind in enumerate(media):
    temp = [el for ls in notBPMeanFP[ind] for el in ls]
    data.append([j[1] for j in temp])
    colors2.append(colors[ind])
bplot = plt.boxplot(data, labels=media, notch=1, patch_artist=True)

plt.ylim(0, 5500)
for patch, color in zip(bplot['boxes'], colors2[:]):
    patch.set_facecolor(color)
    patch.set_alpha(.7)
    patch.set_linewidth(1)
plt.suptitle('Non Branchpoints GFP intensities on ' + conditionKey)

# =============================================================================
#  Plot distribution for bpts and non bpts
# =============================================================================
gs = gridspec.GridSpec(4, 1)
fig = plt.figure(figsize=(11, 8.5))
plt.show()

temp = []
for h, i in enumerate(media):
    if len(notBPMeanFP[i]):
        temp = [el for ls in notBPMeanFP[i] for el in ls]
        temp2 = [j[1] for j in temp]   # NON BRANCHPOINTS
        temp3 = [el for ls in BranchPtMeanFP[i] for el in ls]
        temp4 = [j[1] for j in temp3]  # BRANCHPOINTS

        ax1 = plt.subplot(gs[h])
        ax2 = plt.subplot(gs[h])
        ax1.set_ylim(0, .0012)
        ax2.set_ylim(0, .0012)
        kernel = sp.gaussian_kde(temp2)
        kernel2 = sp.gaussian_kde(temp4)
        positions = np.linspace(0, 6000, num=120)

        Z = kernel(positions)
        ax1.fill_between(positions, Z, color=colors[i], alpha=0.2)

        ax1.text(
            .5, .9, i, transform=ax1.transAxes, fontsize=14,
            verticalalignment='top')

        ax1.axvline(np.percentile(temp2, 50), color='blue', lw=1.5, alpha=.5)

        Z2 = kernel2(positions)
        ax2.fill_between(
            positions, Z2, color=colors[i], alpha=0.2, hatch="/")

        ax2.text(
            .5, .7, "Non-Branchpoints median=%6d" % np.percentile(temp2, 50),
            transform=ax1.transAxes, fontsize=12, verticalalignment='top')
        ax2.text(
            .5, .5, "Branchpoints median=%6d" % np.percentile(temp4, 50),
            transform=ax1.transAxes, fontsize=12, verticalalignment='top')

        ax2.axvline(np.percentile(temp4, 50), color='red', lw=1.5, alpha=.5)

plt.suptitle(
    'GFP-RFP Pixel Intensities Distributions for Branchpoints (hatched)'
    ' and non branchpoints', fontsize=14)
plt.show()
