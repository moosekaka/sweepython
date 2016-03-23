
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

radius=2.5#radius of averaging

b=glob.glob(os.getcwd()+'\[0-9][0-9]Y*')
#b=(os.getcwd()+'\\102214_a1ra1g_Results')

#b=os.getcwd()
for k in range(len(b)):
    fskel = glob.glob(b[k]+'\\'+'RFP\*kel*')
    fVR=glob.glob(b[k]+'\\'+'RFP\*resample*')
    fVG = glob.glob(b[k]+'\\'+'GFP\*resample*') #GFP voxels #GFP voxels
    media=b[k].rsplit('\\',1)
    R=[];G=[];gN=[]

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
    
################################################################################################
################################################################################################
        ptsOld=dataSkel.GetPoints()
        Cell=vtk.vtkCellArray()
        dim=dataRFPV.GetDimensions()
        intenGFP=dataGFPV.GetPointData().GetScalars().GetTuple1
        intenRFP=dataRFPV.GetPointData().GetScalars().GetTuple1
        skelWidth=dataSkel.GetPointData().GetScalars('Width')   
      
        loc=vtk.vtkPointLocator()
        loc.SetDataSet(dataGFPV)
        loc.BuildLocator()
        result=vtk.vtkIdList()
    
        rawRFP=vtk.vtkDoubleArray()
        rawRFP.SetName("rRFP")
        rawGFP=vtk.vtkDoubleArray()
        rawGFP.SetName("rGFP")

############## add the lines/ cells  for connectivity of skel info #############################    
        for i in range(dataSkel.GetNumberOfLines()):
            #print(i)
            oldL=dataSkel.GetCell(i).GetPoints()
            ptID=dataSkel.GetCell(i).GetPointIds()
            oldPt=[oldL.GetPoint(idx) for idx in range(oldL.GetNumberOfPoints())]
            oldId=[ptID.GetId(id) for id in range(ptID.GetNumberOfIds())]
            Cell.InsertNextCell(len(oldId))
            for j in oldId:
                Cell.InsertCellPoint(j)

############### averaging of pts intensity value surrounding each point in skel #####################

        for n in range(dataSkel.GetNumberOfPoints()):        
            ptOI=tuple(np.ceil(i/0.055) for i in dataSkel.GetPoint(n))
            loc.FindPointsWithinRadius(radius,ptOI,result)
            voxID=[result.GetId(i) for i in range(result.GetNumberOfIds())]
            g=np.mean([intenGFP(m) for m in voxID])
            rawGFP.InsertNextValue(g)
            r=np.mean([intenRFP(m) for m in voxID])
            rawRFP.InsertNextValue(r)
    
        polyData=vtk.vtkPolyData()
        polyData.SetPoints(ptsOld)
        polyData.SetLines(Cell)
        polyData.GetPointData().AddArray(rawRFP)
        polyData.GetPointData().AddArray(rawGFP)
        polyData.GetPointData().AddArray(skelWidth)

        list1=[]
        writer=vtk.vtkPolyDataWriter()
        dirName=fskel[x].rsplit('\\',1)
        fileName=dirName[1]
        fileString=os.getcwd()+'\\'+str(radius)+'_raw_'+fileName
    
        writer=vtk.vtkPolyDataWriter()
        writer.SetFileName(fileString)
        writer.SetInput(polyData)
        writer.Update()
