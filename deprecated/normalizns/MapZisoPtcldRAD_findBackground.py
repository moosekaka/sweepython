
import networkx as nx
import math,vtk,glob,string,os
#import createEdgeNodeList
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
##########################################################################################
####################### make sure YPx_meanZRFP.pkl exist, if not run
#######  MeanZ_GRFP.py first!!!!!!#########################################################
##########################################################################################
def disteuc(p0,p1):
    return math.sqrt(vtk.vtkMath.Distance2BetweenPoints(p1,p0))
def GetZ(Id,Dim):
    return (Id/(Dim[0]*Dim[1]))

s=dict() #dict for RFP background
t=dict() #dict for GFP background
b=glob.glob(os.getcwd()+'\[0-9][0-9]Y*')
#b=[[]];b[0]=os.getcwd()
#for k in range(len(b)):
for k in [0]:
    print('now on '+ b[k]);
    fskel = glob.glob(b[k]+'\\'+'cellsRFP\*kel*')
    fVR=glob.glob(b[k]+'\\'+'cellsRFP\*resample*')
    fVG = glob.glob(b[k]+'\\'+'cellsGFP\*resample*') #GFP voxels #GFP voxels
    R=[];G=[];

    for x in range(len(fskel)):
        print(x)
        reader=vtk.vtkPolyDataReader()
        reader.SetFileName(fskel[x]);reader.Update()
        reader2=vtk.vtkStructuredPointsReader()
        reader2.SetFileName(fVG[x]);reader2.Update()
        reader3=vtk.vtkStructuredPointsReader()
        reader3.SetFileName(fVR[x]);reader3.Update()
            
        dataSkel=reader.GetOutput() # Skel coords value from skel file
        dataGFPV=reader2.GetOutput()#voxels GFP
        dataRFPV=reader3.GetOutput()
    
        maxR=dataRFPV.GetPointData().GetScalars().GetRange()[1]
        minR=dataRFPV.GetPointData().GetScalars().GetRange()[0]    
        maxG=dataGFPV.GetPointData().GetScalars().GetRange()[1]
        minG=dataGFPV.GetPointData().GetScalars().GetRange()[0]   
        intenGFP=dataGFPV.GetPointData().GetScalars().GetTuple1
        intenRFP=dataRFPV.GetPointData().GetScalars().GetTuple1
        
        tempR=[intenRFP(i) for i in range(dataRFPV.GetNumberOfPoints())]
        tempG=[intenGFP(i) for i in range(dataGFPV.GetNumberOfPoints())]
        R.append(np.percentile(tempR,50))
        G.append(np.percentile(tempG,50))
    s.update({fskel[i].rsplit('\\',1)[1]:R[i] for i in range(len(fskel))})
    t.update({fskel[i].rsplit('\\',1)[1]:G[i] for i in range(len(fskel))})

output=open('Background50.pkl','wb')
pickle.dump((s,t),output)
output.close()
