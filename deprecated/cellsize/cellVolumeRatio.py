import glob
import os
import math
import matplotlib.pyplot as plt
import pandas as pd
import cPickle as pickle


def pltBx(X, grphTitle, labs, *args, **kwargs):
    plt.figure()
    plt.boxplot(X, labels=labs, notch=1)
    plt.title(grphTitle)
    plt.show()

b = glob.glob(os.getcwd() + '\[0-9][0-9]Y*')
cellMitoVol = []
input = open('CellVol.pkl', 'rb')
cellVol = pickle.load(input)
input.close()
Vols = [i[1] for i in cellVol]
labels = [i[0] for i in cellVol]

for k in range(len(b)):
    media = b[k].rsplit('\\', 2)
    mediaCur = media[-1].split('_')[0]
    A = glob.glob(b[k] + '\\' + '*RFP\*gnet')

    SkelVol = []

    for i in range(len(A)):
        x = pd.read_table(A[i], skiprows=1)
        temp = x.sum()[2]
        temp2 = x.keys()
        temp3 = temp + float(temp2[-1])
        SkelVol.append(temp3 * math.pi * .15 * .15)

    cellMitoVol.append(SkelVol)

cellVolRatio = []
for k in range(len(b)):
    temp = [cellMitoVol[k][m] / cellVol[k][1][m]
            for m in range(len(cellMitoVol[k]))]

    cellVolRatio.append(temp)

pltBx(Vols, 'cellVol', labels)
pltBx(cellMitoVol, 'mitolVol', labels)
pltBx(cellVolRatio, 'VolRatio', labels)
