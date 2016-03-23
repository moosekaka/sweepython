import math,vtk,glob,string,os
import scipy.stats as sp
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
from itertools import chain
import pandas as pd
import matplotlib.gridspec as gridspec
from tvtk.api import tvtk

fileLoc=os.getcwd()+'\N*vtk'
files = glob.glob(fileLoc)
parDir=os.path.dirname(os.getcwd())

r=pd.read_csv(parDir+'\\'+'BackgroundRFPGFP.csv')
g=pd.read_csv(parDir+'\\'+'BackgroundRFPGFP.csv')
backgroundGFP={i[0]:i[1] for i in zip(g['Label'],g['MeanGFP'])}
backgroundRFP={i[0]:i[1] for i in zip(r['Label'],r['MeanRFP'])}
plt.close('all')

with open('05a1ra1g_cellData.pkl','rb') as inpt:
    EnrichBP,DiffBP,arrBrc,BrcMeans,cellAPmeans=pickle.load(inpt)
with open('05a1ra1g_grph.pkl','rb') as inpt:
    G=pickle.load(inpt)
G=G[2]

for a in DiffBP:
    curGrph=G[a]
    reader=tvtk.PolyDataReader()
    reader.set(file_name=files[a])
    reader.update()
    data=reader.output
    scalarsNorm=data.point_data.scalars
    temp=data.point_data
    rawGFP=np.ravel(temp.get_array('rGFP'))
    rawRFP=np.ravel(temp.get_array('rRFP'))
    WidthEq=np.ravel(temp.get_array('WidthEq'))
    WidthEq2=np.ravel(temp.get_array('WidthEq2'))
    scalarsNorm2=np.ravel(temp.get_array('DY2_minmax'))
    fileKey=files[a].rsplit('\\',1)[1][5:][:-13]
    minA=backgroundRFP[fileKey]
    minB=min(backgroundGFP[fileKey],min(rawGFP))

    Norm=[];GFP=[];RFP=[];arrPts=[];W=[];W2=[];Norm2=[]
    rGFP=[];lineId={}

############## get the line of interest ###########################\
    bpts=[curGrph.node[i]['coord'] for i in curGrph.nodes()
    if curGrph.node[i]['degree']>1]
    bptsId = [data.find_point(el) for el in bpts]
    
    epts=[curGrph.node[i]['coord'] for i in curGrph.nodes()
    if curGrph.node[i]['degree']==1]
    eptsId = [data.find_point(el) for el in epts]
       
    firstId=0
    for j in range(data.number_of_lines):
        cellIds=list(data.get_cell(j).point_ids)
        #first and last point on cell        
        temp=[cellIds[0],cellIds[-1]]
    
        Norm.append([scalarsNorm[i] for i in cellIds])    
        GFP.append([rawGFP[i] for i in cellIds])
        RFP.append([rawRFP[i] for i in cellIds])
        W.append([WidthEq[i] for i in cellIds])
        W2.append([WidthEq2[i] for i in cellIds])
        Norm2.append([scalarsNorm2[i] for i in cellIds])   
        
        lineId.setdefault(firstId,[])
        for i in temp:
            if i in bptsId:        
                lineId[firstId].append('b')       
            elif i in eptsId:                
                lineId[firstId].append('e')
        firstId+=len(cellIds)

############################    PLOT  ########################################

    fig,axarr = plt.subplots(6, sharex=True,figsize=(11,8.25))
    axarr[0].set_xlim(-1,firstId+1)
    axarr[0].set_title('RFP Width$_{eq}$',fontsize=11)
    axarr[0].plot([el for lis in W for el in lis],color='red')
    
    axarr[2].set_title('GFP Width$_{eq2}$',fontsize=11)
    axarr[2].plot([el for lis in W2 for el in lis],color='green')
    
    axarr[1].set_title('RFP$_{Background\ Substracted}$',fontsize=11)
    axarr[1].plot([el-minA for lis in RFP for el in lis],color='red',alpha=0.8)         
    
    axarr[3].set_title('GFP$_{Background\ Substracted}$',fontsize=11)
    axarr[3].plot([el-minB for lis in GFP for el in lis],color='green',alpha=0.8)   
    
    axarr[4].set_title('$\Delta \Psi$ GFP$_{Background\ Substracted}$ /\
    RFP$_{widthEq}$ ',fontsize=11)
    axarr[4].plot([el for lis in Norm for el in lis],color='purple')
    
    axarr[5].set_title('$\Delta \Psi$ RFP$_{Background\ Substracted}$ /\
    GFP$_{widthEq2}$',fontsize=11)
    axarr[5].set_xlabel(fileKey)
    axarr[5].plot([el for lis in Norm2 for el in lis],color='purple',alpha=0.7)

    plt.show()
    lineIdl=chain(sorted(lineId.keys())) #sort of like a que object
    thisel=lineIdl.next();# first elem in pump
    
    for indx,elem in enumerate(lineId):         
        axarr[4].axvline(thisel,linestyle='-',color='k',linewidth='.75',alpha=0.4)
        axarr[4].axhline(cellAPmeans[a],linestyle='-',color='k',linewidth='.5',alpha=0.8)
        axarr[5].axvline(thisel,linestyle='-',color='k',linewidth='.75',alpha=0.4)
        axarr[5].axhline(cellAPmeans[a],linestyle='-',color='k',linewidth='.5',alpha=0.8)
        
        if(indx%2):
            height=0.96 
        else:
            height=0.92
            
        if lineId[thisel][0]=='b':
            axarr[4].plot(thisel,height,marker='<',color='m',linewidth='.75',alpha=1,ms=4,mec='m')
        else:
            axarr[4].plot(thisel,height,marker='<',color='c',linewidth='.75',alpha=1,ms=4,mec='c')
        try:            
            nextel=lineIdl.next()                       
            if lineId[thisel][1]=='b':
                axarr[4].plot(nextel,height,marker='>',color='m',linewidth='.75',alpha=1,ms=4,mec='m')
            else:
                axarr[4].plot(nextel,height,marker='>',color='c',linewidth='.75',alpha=1,ms=4,mec='c')
            
            thisel=nextel

        except StopIteration: # no last elt to do stuff to
            pass
    if lineId[thisel][1]=='b':
        axarr[4].plot(firstId,height,marker='>',color='m',linewidth='.75',alpha=1,ms=4,mec='m')
    else:
        axarr[4].plot(firstId,height,marker='>',color='c',linewidth='.75',alpha=1,ms=4,mec='c')
        
    fig.savefig(fileKey+'.png')
    plt.close()