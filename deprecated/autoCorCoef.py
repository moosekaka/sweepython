import os
import fnmatch
import numpy as np
import cPickle as pickle
from tvtk.api import tvtk


def estimated_autocorrelation(x):
    X = np.array(x)
    n = len(X)
    variance = X.var()
    X = X-X.mean()
    r = np.correlate(X, X, mode='full')[-n:]
    result = r/(variance*(np.arange(n, 0, -1)))
    return (result, n)

# =============================================================================
#               init vars and get vtk,graph data
# =============================================================================
vtkF = {}
G = {}
backgroundGFP = {}
backgroundRFP = {}
Type = {}
Dates = {}
rfpWidthCorCoef = {}
autoCorCoef = {}
parDir = os.path.dirname(os.getcwd())
for root, dirs, files in os.walk(os.getcwd()):
    for i in files:
        if fnmatch.fnmatch(i, '*vtk'):
            vtkF.setdefault(root.rsplit('\\', 1)[1], []).append(
                os.path.join(root, i))
        if fnmatch.fnmatch(i, '*grph.pkl'):
            G.setdefault(root.rsplit('\\', 1)[1], []).append(
                os.path.join(root, i))

media = sorted(vtkF.keys())

#    get metadatas
with open(parDir+'\\'+'fileMetas.pkl', 'rb') as inpt:
    s = pickle.load(inpt)
for i in s:
    backgroundRFP[i] = s[i][0]
    backgroundGFP[i] = s[i][1]
    Type[i] = i[:4]
    Dates[i] = i.split('_')[1][:]

# =============================================================================
#               begin pipeline
# =============================================================================
for mem in media:
    for n, a in enumerate(vtkF[mem]):
        Norm = []
        GFP = []
        RFP = []
        arrPts = []
        W = []
        W2 = []
        rGFP = []
        lineId = {}
        R = []

        reader = tvtk.PolyDataReader()
        reader.set(file_name=a)
        reader.update()
        data = reader.output
        scalarsNorm = data.point_data.scalars
        temp = data.point_data
        rawGFP = np.ravel(temp.get_array('rGFP'))
        rawRFP = np.ravel(temp.get_array('rRFP'))
        dyRaw = np.ravel(temp.get_array('DY_raw'))
        WidthEq = np.ravel(temp.get_array('WidthEq'))
        fileKey = a.rsplit('\\', 1)[1][5:][:-13]

        if backgroundRFP[fileKey] > min(rawRFP):
            minA = backgroundRFP[fileKey]-1
        else:
            minA = backgroundRFP[fileKey]-1
        minB = min(backgroundGFP[fileKey], min(rawGFP))

#   Estimate edge autocorrelation Rtemp
        for j in range(data.number_of_lines):
            cellIds = list(data.get_cell(j).point_ids)
            temp = [cellIds[0], cellIds[-1]]  # first and last point on cell
            Norm.append([scalarsNorm[i] for i in cellIds])
            GFP.append([rawGFP[i] for i in cellIds])
            RFP.append([rawRFP[i] for i in cellIds])
            W.append([WidthEq[i] for i in cellIds])
            Rtemp = estimated_autocorrelation(
                [dyRaw[i] for i in cellIds])
            R.append(Rtemp[0])

        autoCorCoef.setdefault(mem, []).append(
            [i for i in R if len(i)])

with open('autoCors.pkl', 'wb') as output:
    pickle.dump(autoCorCoef, output)
# =============================================================================
#               Plot
# =============================================================================
# plt.close('all')
# gs = gridspec.GridSpec(2, 2)
# fig = plt.figure(figsize=(11, 8.25))
#
# for h, i in enumerate(media):
#     X = np.ravel([el for lis in autoCorCoef[i] for el in lis])
#     for i in X:
#         ax = plt.subplot(gs[h])
#         n = len(i)
#         if n > 5:
#             ax.plot(i[:15])
# plt.show()
