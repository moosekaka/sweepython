#####Build Vtk skel files with normalized values for DY using TVTK####
import networkx as nx
import math,glob,string,os
import scipy.stats as sp
import numpy as np
import cPickle as pickle
import pandas as pd
from tvtk.api import tvtk

backgroundGFP={};backgroundRFP={};Type={};Dates={}
parDir=os.path.dirname(os.getcwd())
parparDir=os.path.dirname(parDir)
b=os.getcwd()
#files=glob.glob(os.getcwd()+'\\*raw*vtk')
files=glob.glob(os.getcwd()+'\\*vtk')

with open(parparDir+'\\'+'fileMetas.pkl','rb') as inpt:
    s=pickle.load(inpt)

for i in s:    
    backgroundRFP[i]=s[i][0]
    backgroundGFP[i]=s[i][1]
    Type[i]=i[:4]
    Dates[i]=i.split('_')[1][:]
#r=pd.read_csv(b+'\\'+'BackgroundRFPGFP.csv')
#g=pd.read_csv(b+'\\'+'BackgroundRFPGFP.csv')
#backgroundGFP={i[0]:i[1] for i in zip(g['Label'],g['meanGFP'])}
#backgroundRFP={i[0]:i[1] for i in zip(r['Label'],r['meanRFP'])}
#with open(parDir+'\\'+'fileDates.pkl','rb') as inpt:
#    s=pickle.load(inpt)

minmaxRFP=[];minmaxGFP=[]
dates=[];bckgrndGFP=[];bckgrndRFP=[]; 
################# Read raw vtk files as input ############################################
for a in [21]: 
    reader=tvtk.PolyDataReader()
    reader.set(file_name=files[a])
    reader.update()
    data=reader.output
    temp=data.point_data
    rawGFP=np.ravel(temp.get_array('rGFP'))
    rawRFP=np.ravel(temp.get_array('rRFP'))

    fileKey=files[a].rsplit('\\',1)[1][5:-13]   
    minmaxRFP.append((min(rawRFP),max(rawRFP)))
    minmaxGFP.append((min(rawGFP),max(rawGFP)))   
    if(backgroundRFP[fileKey]>minmaxRFP[0][0]):
        #ensures minimum values of 1 to avoid big divisions
        minA=minmaxRFP[0][0]-1
    else:
        #ensures minimum values of 1 to avoid big divisions
        minA=backgroundRFP[fileKey]-1

    minB=min(backgroundGFP[fileKey],minmaxGFP[0][0])
    bckgrndGFP.append(minB);
    bckgrndRFP.append(minA);
    #dates.append(s[fileKey])
    
#############delete cells before normalize ############################################
    pts=data.points
    lnsb4dlt=data.lines
    
    polyData=tvtk.PolyData()
    polyData.points=pts
    polyData.lines=lnsb4dlt
    polyData.build_links()
    for n in [35]:
        polyData.delete_cell(n)
    polyData.remove_deleted_cells()
    
    #new lines from polyData object with deleted cells
    pointIds=[];filteredLines=[];
    ######## method for getting pointIds as unique list ##########
    for el in range(polyData.number_of_lines):
        temp=np.ravel(polyData.get_cell(el).point_ids)
        filteredLines.append(temp)
    ##############################################################
    pointIds=np.unique([el for ln in filteredLines for el in ln]  )
    #background Substract
    A=rawRFP-minA
    B=rawGFP-minB   
    minAb=np.min(A)
    #width equivalent
    W=A/minAb              
    #raw DY/W normalized values
    DY=B/W      
    #rescale DY to minmax
    minDY=min([DY[i] for i in pointIds]);
    maxDY=max([DY[i] for i in pointIds]);      
    normDY=((DY-minDY)/(maxDY-minDY))                       
        
##################################Make VTK file ################################################
    
    
    
    

    polyData.point_data.scalars=normDY
    polyData.point_data.scalars.name='DY_minmax'
    polyData.point_data.add_array(W)
    polyData.point_data.get_array(1).name='WidthEq'
    polyData.point_data.add_array(DY)
    polyData.point_data.get_array(2).name='DY_raw'
    polyData.point_data.add_array(rawRFP)
    polyData.point_data.get_array(3).name='rRFP'
    polyData.point_data.add_array(rawGFP)
    polyData.point_data.get_array(4).name='rGFP'
    polyData.point_data.add_array(A)
    polyData.point_data.get_array(5).name='bkstRFP'
    polyData.point_data.add_array(B)
    polyData.point_data.get_array(6).name='bkstGFP'
    polyData.point_data.update()  
    
    writer=tvtk.PolyDataWriter()
    fileString=os.getcwd()+'\\Norm_'+fileKey+'_skeleton.vtk'    
    writer=tvtk.PolyDataWriter()
    writer.set(file_name=fileString)
    writer.set_input(polyData)
    writer.update()
        
    
        


