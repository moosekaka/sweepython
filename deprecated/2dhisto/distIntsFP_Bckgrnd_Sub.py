import os
import fnmatch
import scipy.stats as sp
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from tvtk.api import tvtk

G = {}
Bpoints = {}
Graphs = {}
vtkF = {}
meanFP = {}
Type = {}
Dates = {}
meanFP = {}
backgroundGFP = {}
backgroundRFP = {}
bckgrndGFP = {}
bckgrndRFP = {}
parDir = os.path.dirname(os.getcwd())
conditionKey = 'All Dates'
plt.close('all')

colors = {'YPD': 'b', 'YPE': 'g', 'YPR': 'r', 'YPL': 'm'}
condition = {'All Dates': ['042515', '042715', '042915', 'c42515', '052315'],
             '04-27-15': ['042715'],
             '04-25-15': ['042515', 'c42515'],
             '04-29-15': ['042915'],
             'c4-25-15': ['c42515'],
             '05-23-15': ['052315']}

with open(parDir + '\\fileMetas.pkl', 'rb') as inpt:
    s = pickle.load(inpt)

parDir = os.path.dirname(os.getcwd())

for root, dirs, files in os.walk(os.getcwd()):
    for i in files:
        if fnmatch.fnmatch(i, '*grph.pkl'):
            G.setdefault(
                root.rsplit('\\', 1)[1], []).append(
                os.path.join(root, i))
        if fnmatch.fnmatch(i, '*vtk'):
            vtkF.setdefault(
                root.rsplit('\\', 1)[1], []).append(
                os.path.join(root, i))

media = sorted(G.keys())

for i in s:
    backgroundRFP[i] = s[i][0]
    backgroundGFP[i] = s[i][1]
    Dates[i] = i.split('_')[1][:]

for mem in media:
    for n, a in enumerate(vtkF[mem]):
        fileKey = a.rsplit('\\', 1)[1][5:-13]
        reader = tvtk.PolyDataReader()
        reader.set(file_name=a)
        reader.update()
        data = reader.output
        temp = data.point_data
        rawGFP = np.ravel(temp.get_array('rGFP'))
        rawRFP = np.ravel(temp.get_array('rRFP'))

#      Establish background of channel#########

        if(backgroundRFP[fileKey] > min(rawRFP)):
            minA = backgroundRFP[fileKey] - 1
        else:
            minA = backgroundRFP[fileKey] - 1
            minB = min(backgroundGFP[fileKey], min(rawGFP))

        bckgrndGFP[fileKey] = minB
        bckgrndRFP[fileKey] = minA

#      subtract the raw FP channels with the background of channel#

        meanFP.setdefault(mem, [])
        conditionKey = 'All Dates'
        if Dates[fileKey] in condition[conditionKey]:
            # without background subtraction
            #RFPtemp=[el for el in rawRFP]
            #GFPtemp=[el for el in rawGFP]
            # with background subtraction
            RFPtemp = [el - bckgrndRFP[fileKey] for el in rawRFP]
            GFPtemp = [el - bckgrndGFP[fileKey] for el in rawGFP]
            temp = zip(RFPtemp, GFPtemp)
            meanFP[mem].append(temp)

#      individual channel GFP histograms
gs = gridspec.GridSpec(4, 1)
fig = plt.figure(figsize=(17, 11))
plt.show()

temp = []
for h, i in enumerate(media):
    if len(meanFP[i]):

        temp = [el for ls in meanFP[i] for el in ls]
        tempFilt2 = [j[1] for j in temp]
        ax2 = plt.subplot(gs[h])
        ax2.set_ylim(0, 0.0010)
        kernel2 = sp.gaussian_kde(tempFilt2)
        positions = np.linspace(2000, 7500, num=55)

        Z2 = kernel2(positions)
        ax2.fill_between(positions, Z2, color=colors[i], alpha=0.3)
        ax2.text(0.5, 0.9, i, transform=ax2.transAxes, fontsize=14,
                 verticalalignment='top')

plt.suptitle('2d GFP-RFP Pixel Intensities Distributions on ' + conditionKey)

#    Gaussian KDEs  GFP-RFP
gs2 = gridspec.GridSpec(2, 2)
fig2 = plt.figure(figsize=(17, 11))

xmin = 0
ymin = 0
xmax = 5500
ymax = 5500
for h, i in enumerate(media):
    if len(meanFP[i]):
        temp = [el for lis in meanFP[i] for el in lis]
        y = [j[1] for j in temp]
        x = [j[0] for j in temp]

        X, Y = np.mgrid[xmin:xmax:55j, ymin:ymax:55j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        kernel = sp.gaussian_kde([x, y])
        Z = np.reshape(kernel(positions).T, X.shape)

        axes = plt.subplot(gs2[h])
        axes.imshow(
            np.rot90(Z),
            cmap=plt.cm.coolwarm,
            extent=[xmin, xmax, ymin, ymax])
        axes.set_title(i)
        axes.set_aspect('equal')
plt.suptitle('2d GFP-RFP Pixel Intensities Distributions on ' + conditionKey)
plt.show()

#    boxplot
data = []
plt.figure(figsize=(11, 8.25))
colors2 = []
for h, ind in enumerate(media):
    temp = [el for ls in meanFP[ind] for el in ls]
    data.append([j[1] for j in temp])
    colors2.append(colors[ind])
bplot = plt.boxplot(
    data,
    labels=media,
    notch=1,
    patch_artist=True,
    bootstrap=100,
    showmeans=True,
    meanline=0)
plt.show()
plt.ylim(2000, 7500)
plt.title(conditionKey)
for patch, color in zip(bplot['boxes'], colors2[:]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
    patch.set_linewidth(1)
