import pandas as pd
import matplotlib.pyplot as plt
import math,vtk,glob,string,os
import cPickle as pickle

b=glob.glob(os.getcwd()+'\[0-9][0-9]Y*')
cellSurfArea=[]
budUnbud=[]

for k in range(len(b)):
    media=b[k].rsplit('\\',2)
    mediaCur=media[-1].split('_')[0]
    Data=pd.read_table(b[k]+'\BF\Results.txt')

    A=Data['Label']
    Cells=[i.split(':')[2] for i in A]
    Ma=Data['Major'];Mi=Data['Minor'];
    im=0;l=0;L=[[] for i in range(119)]

    D={};
    for i,j in enumerate(Cells):
        D.setdefault((j.split('_')[-3]+j.split('_')[-1]),[]).append((j,Mi[i],Ma[i]))
    momBud=D.values()
    momBud.sort()
    budUnbud.append(momBud)
    CellArea=[]
    for i in momBud:
        temp=0
        for j in i:
            temp+=math.pi*(j[2]*.055/2)*(j[1]*.055/2)
        CellArea.append(temp)

    cellSurfArea.append((mediaCur[2:],CellArea))

Out=cellSurfArea
output=open('CellSurfArea.pkl','wb')
pickle.dump(Out,output)
output.close()
