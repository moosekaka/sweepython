
import networkx as nx
import math,vtk,glob,string,os
#import createEdgeNodeList
import numpy as np
import matplotlib.pyplot as plt
import progressbar
import cPickle as pickle
##########################################################################################
####################### make sure YPx_meanZRFP.pkl exist, if not run
#######  MeanZ_GRFP.py first!!!!!!#########################################################
##########################################################################################
def disteuc(p0,p1):
    return math.sqrt(vtk.vtkMath.Distance2BetweenPoints(p1,p0))
def GetZ(Id,Dim):
    return (Id/(Dim[0]*Dim[1]))


b=glob.glob(os.getcwd())
for k in range(1,2):
    fskel = glob.glob('05*R*skel*')
    #fVR=glob.glob(b[k]+'\\'+'RFP\*resample*')
    #fVG = glob.glob(b[k]+'\\'+'GFP\*resample*') #GFP voxels #GFP voxels
    #media=b[k].rsplit('\\',1)
    R=[];G=[];gN=[]


    for x in range(len(fskel)):
        print(x)
  
        reader=vtk.vtkPolyDataReader()
        reader.SetFileName(fskel[x]);reader.Update()
        #reader2=vtk.vtkStructuredPointsReader()
        #reader2.SetFileName(fVG[x]);reader2.Update()
        #reader3=vtk.vtkStructuredPointsReader()
        #reader3.SetFileName(fVR[x]);reader3.Update()
            
        dataSkel=reader.GetOutput() # Skel coords value from skel file
        #dataGFPV=reader2.GetOutput()#voxels GFP
        #dataRFPV=reader3.GetOutput()
       
################################################################################################
        ptsOld=dataSkel.GetPoints()
        Cell=vtk.vtkCellArray()
        #dim=dataRFPV.GetDimensions()
 
        skelInt=dataSkel.GetPointData().GetScalars('Intensity')   
        
        for i in range(dataSkel.GetNumberOfLines()):
            #print(i)
            oldL=dataSkel.GetCell(i).GetPoints()
            ptID=dataSkel.GetCell(i).GetPointIds()
            oldPt=[oldL.GetPoint(idx) for idx in range(oldL.GetNumberOfPoints())]
            oldId=[ptID.GetId(id) for id in range(ptID.GetNumberOfIds())]
            Cell.InsertNextCell(oldL.GetNumberOfPoints())
            for j in oldId:
                Cell.InsertCellPoint(j)
     
        polyData=vtk.vtkPolyData()
        polyData.SetPoints(ptsOld)
        polyData.SetLines(Cell)
        polyData.GetPointData().SetScalars(skelInt)
        
        writer=vtk.vtkPolyDataWriter()
        fileName=fskel[x].rsplit('.',1)[0]
        #fileName=dirName[1]
        fileString=string.join(('raw',fileName),sep='_')+'.vtk'
        

    
        writer=vtk.vtkPolyDataWriter()
        writer.SetFileName(fileString)
        writer.SetInput(polyData)
        writer.Update()
