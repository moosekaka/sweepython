import networkx as nx
import math,vtk,glob,string,os
import scipy.stats as sp
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
from itertools import cycle
import pandas as pd
from tvtk.api import tvtk

    
r=pd.read_csv('BackgroundRFPGFP.csv')
g=pd.read_csv('BackgroundRFPGFP.csv')
backgroundGFP={i[0]:i[1] for i in zip(g['Label'],g['MeanGFP'])}
backgroundRFP={i[0]:i[1] for i in zip(r['Label'],r['MeanRFP'])}         

b=os.getcwd()
plt.close('all')

minmaxRFP=[];minmaxGFP=[];
meanRFP=[];meanGFP=[];
bckgrndGFP=[];bckgrndRFP=[];


for k in [0]:    
    fileLoc=b+'\\2.5_raw*vtk'
    files = glob.glob(fileLoc)      
    
    for a in range(len(files)):
        
        reader=tvtk.PolyDataReader()
        reader.set(file_name=files[a])
        reader.update()
        data=reader.output
        temp=data.point_data
        
        rawRFP=np.ravel(temp.get_array("rRFP"))
        rawGFP=np.ravel(temp.get_array("rGFP"))

############ the min, max and mean values of  the skeletons are here ###########################
        meanRFP.append(np.mean(rawRFP))
        meanGFP.append(np.mean(rawGFP))

        fileKey=files[a].rsplit('\\',1)[1][8:-13]
        minmaxRFP.append((min(rawRFP),max(rawRFP)))
        minmaxGFP.append((min(rawGFP),max(rawGFP))) 
        #background labels are switched becuase Micromanager mislabelled RFP as GFP
        minA=backgroundGFP[fileKey]
        minB=backgroundRFP[fileKey]
        bckgrndGFP.append(minB);
        bckgrndRFP.append(minA);

################################Normalize #####################################################
        pts=data.points
        lns=data.lines
        pointIds=[];lines=[];
        ######## method for getting pointIds as unique list ##########
        for el in range(data.number_of_lines):
            temp=np.ravel(data.get_cell(el).point_ids)
            lines.append(temp)
        ##############################################################
        pointIds=np.unique([el for ln in lines for el in ln]  )
           
        A=rawRFP-minA
        B=rawGFP-minB
        
        minAb=min(A)
        minBa=min(B)
        ## W is the fold change for the RFP channel
        W=A/minAb 
        W2=B/minBa            
    
        DY=B/W 
        DY2=A/W2     
        #rescale DY to minmax      
        #rescale DY to minmax
        minDY=min([DY[i] for i in pointIds]);
        maxDY=max([DY[i] for i in pointIds]);      
        normDY=((DY-minDY)/(maxDY-minDY))  
        
        minDY2=min([DY2[i] for i in pointIds]);
        maxDY2=max([DY2[i] for i in pointIds]);      
        normDY2=((DY2-minDY2)/(maxDY2-minDY2))           
##################################Make VTK file #####################################################            
        polyData=tvtk.PolyData()
        polyData.points=pts
        polyData.lines=lns
        polyData.point_data.scalars=normDY
        polyData.point_data.scalars.name='DY_minmax'
        polyData.point_data.add_array(W)
        polyData.point_data.get_array(1).name='WidthEq'
        polyData.point_data.add_array(DY)
        polyData.point_data.get_array(2).name='DY_raw'
        polyData.point_data.add_array(W2)
        polyData.point_data.get_array(3).name='WidthEq2'
        polyData.point_data.add_array(rawRFP)
        polyData.point_data.get_array(4).name='rRFP'
        polyData.point_data.add_array(rawGFP)
        polyData.point_data.get_array(5).name='rGFP'
        polyData.point_data.add_array(normDY2)
        polyData.point_data.get_array(6).name='DY2_minmax'
        polyData.point_data.add_array(DY2)
        polyData.point_data.get_array(7).name='DY2_raw'
        
        polyData.point_data.update()      
        
        writer=tvtk.PolyDataWriter()
        fileString=os.getcwd()+'\\Norm_'+fileKey+'_skeleton.vtk'    
        writer=tvtk.PolyDataWriter()
        writer.set(file_name=fileString)
        writer.set_input(polyData)
        writer.update()
            
        
            


