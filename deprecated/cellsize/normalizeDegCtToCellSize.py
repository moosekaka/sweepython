#script to plot the number of brancpoints normalized to cell VOl, \
#number of cell Vols and mito graphs are different therefore have to 
#use a dictionary to find common cells in both sets
 
import sys,math,glob,string,os
import numpy as np
import scipy.stats as sp
import matplotlib.pyplot as plt
import cPickle as pickle
import scipy.stats as sp

with open('CellVolFileNames.pkl') as inpt:
   VolName=pickle.load(inpt)

with open('MitoGraphFileNames.pkl') as inpt:
   GraphName=pickle.load(inpt)

with open('CellVol.pkl') as input:
    cellVol=pickle.load(input)
    Vols=[i[1] for i in cellVol]
    
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

#normalized the # branchpoints to Cell Vols
normDegCt=[];
for i in range(5):
    temp=[]
    for key,value in GrphFilt[i].iteritems():
        temp.append(degMCount[i][value]/Vols[i][VolFilt[i][key]])
    normDegCt.append(temp)

plt.figure();labs=['YPD','YPE','YPR','YPL','YPLatg32']
plt.boxplot([i for i in normDegCt],labels=labs,notch='true');plt.show()

#normalized the # branchpoints to Mito Vols
sumED3=[[i*math.pi*.15*.15 for i in j] for j in sumED2]
normDegMito=[[degMCount[a][i]/sumED3[a][i] for i in range(len(sumED3[a]))] for a in range(5)]
plt.figure();
plt.boxplot([i for i in normDegMito],labels=labs,notch='true');plt.show()


