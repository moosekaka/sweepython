# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as sp
import math,vtk,glob,string,os
import numpy as np
import cPickle as pickle
from itertools import cycle


b=glob.glob(os.getcwd()+'\[0-9][0-9]Y*')

averageDeg=[]
for k in range(len(b)):
    temp=glob.glob(b[k])
    with open(temp[0],'rb') as inpt:
        nodes,edges,G=pickle.load(inpt)
    temp2=[2*float(len(edges[i]))/len(nodes[i]) for i in range(len(G))]
    #print(len(temp2))
    averageDeg.append((k,temp2))

with open('CellVolFileNames.pkl') as inpt:
   VolName=pickle.load(inpt)

with open('MitoGraphFileNames.pkl') as inpt:
   GraphName=pickle.load(inpt)

with open('CellSurfArea.pkl') as input:
    cellSurfArea=pickle.load(input)
    area=[i[1] for i in cellSurfArea]

with open('degBranchpointsCt.pkl') as inpt:
   degM=pickle.load(inpt)

with open('MitoLengths.pkl') as inpt:
    sumED2=pickle.load(inpt)
#degBranchpointsCt is obtained from GRaphStats.py
degMCount=[[len(i) for i in degM[j]] for j in range(5)]

#order the names of the common mitograph and cellVol filenames as a
#dictionary to make referencing common cells easier
GrphFilt=[{i:h for h,i in enumerate(GraphName[j]) if i in VolName[j]} for j in range(5)]
VolFilt=[{i:h for h,i in enumerate(VolName[j]) if i in GraphName[j]} for j in range(5)]

#quasiDen = mitolenght/cellSurfArea
quasiDen=[];aveDegFilt=[];
for i in range(5):
    temp=[]
    tempD=[]
    for key,value in GrphFilt[i].iteritems():
        temp.append(sumED2[i][value]/area[i][VolFilt[i][key]])
        tempD.append(averageDeg[i][1][value])
    quasiDen.append(temp)
    aveDegFilt.append(tempD)

plt.figure();labs=['YPD','YPE','YPR','YPL','YPLatg32']
plt.boxplot([i for i in quasiDen],labels=labs,notch='true');
plt.title('quasiDen');plt.show()

plt.figure();labs=['YPD','YPE','YPR','YPL','YPLatg32']
plt.title('avgDeg');plt.show()
plt.boxplot([i for i in aveDegFilt],labels=labs,notch='true');

colors = ['b','g','r','m','c']
cycler =cycle(colors)
plt.figure();
plt.xlim(0.5,4.25);plt.ylim(0.5,3.5)
for i in range(5):



    plt.plot(quasiDen[i],aveDegFilt[i],'o',color=cycler.next(),label=labs[i])

    plt.xlabel('quasiDensity');plt.ylabel('averageDegree')
    plt.title(str(labs[i]));plt.show()
    plt.legend(labs,loc='lower right',numpoints=1);plt.show()

print('Media ,Mean <k>/q,CI_low,CI_high,Rvalue')
for i in range(5):
    z = np.polyfit(quasiDen[i],aveDegFilt[i], 1)
    p = np.poly1d(z)
    plt.plot(quasiDen[i],p(quasiDen[i]),'-',color=cycler.next())
    slope, intercept, r_value, p_value, std_err =sp.linregress(quasiDen[i],aveDegFilt[i])     
    print(labs[i]+',%.3f,%.3f,%.3f,%.3f' \
    %(slope,slope-1.96*std_err,slope+1.96*std_err,r_value))
    #print ('y=%.3fx+%.3f'%(z[0],z[1]))
