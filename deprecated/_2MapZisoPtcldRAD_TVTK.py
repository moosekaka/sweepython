"""
       Build vtk skel files with normalized values for DY using TVTK
"""
import glob
import os
import numpy as np
import cPickle as pickle
from tvtk.api import tvtk

# pylint: disable=C0103
backgroundGFP = {}
backgroundRFP = {}
Type = {}
Dates = {}
minmaxRFP = []
minmaxGFP = []
dates = []
bckgrndGFP = []
bckgrndRFP = []
parDir = os.path.dirname(os.getcwd())
parparDir = os.path.dirname(parDir)
b = os.getcwd()
files = glob.glob(os.getcwd()+'\\*vtk')
with open(parparDir+'\\'+'fileMetas.pkl', 'rb') as inpt:
    s = pickle.load(inpt)
# pylint: enable=C0103

for i in s:
    backgroundRFP[i] = s[i][0]
    backgroundGFP[i] = s[i][1]

# =============================================================================
#       Read raw vtk files as input
# =============================================================================
for a in range(len(files)):
    print"working..."
    reader = tvtk.PolyDataReader()
    reader.set(file_name=files[a])
    reader.update()
    data = reader.output
    temp = data.point_data
    rawGFP = np.ravel(temp.get_array('rGFP'))
    rawRFP = np.ravel(temp.get_array('rRFP'))
    tubeWidth = np.ravel(temp.get_array('Width'))
    fileKey = files[a].rsplit('\\', 1)[1][:-13]
    minmaxRFP.append((min(rawRFP), max(rawRFP)))
    minmaxGFP.append((min(rawGFP), max(rawGFP)))
    if backgroundRFP[fileKey] > minmaxRFP[a][0]:

        #   ensures minimum values of 1 to avoid big divisions
        minA = minmaxRFP[a][0]-1
    else:
        #   ensures minimum values of 1 to avoid big divisions
        minA = backgroundRFP[fileKey]-1

    minB = min(backgroundGFP[fileKey], minmaxGFP[a][0])
    bckgrndGFP.append(minB)
    bckgrndRFP.append(minA)

#   Normalize
    pts = data.points
    lns = data.lines
    pointIds = []
    lines = []
#   method for getting pointIds as unique list
    for el in range(data.number_of_lines):
        temp = np.ravel(data.get_cell(el).point_ids)
        lines.append(temp)
    pointIds = np.unique([el for ln in lines for el in ln])

#   background Substract
    A = rawRFP-minA
    B = rawGFP-minB
    minAb = np.min(A)
#    width equivalent
    W = A/minAb
#   raw DY/W normalized values
    DY = B/W
#   rescale DY to minmax
    minDY = min([DY[i] for i in pointIds])
    maxDY = max([DY[i] for i in pointIds])
    normDY = ((DY-minDY)/(maxDY-minDY))

# =============================================================================
#           Make VTK file
# =============================================================================
    polyData = tvtk.PolyData()
    polyData.points = pts
    polyData.lines = lns
    polyData.point_data.scalars = normDY
    polyData.point_data.scalars.name = 'DY_minmax'
    polyData.point_data.add_array(W)
    polyData.point_data.get_array(1).name = 'WidthEq'
    polyData.point_data.add_array(DY)
    polyData.point_data.get_array(2).name = 'DY_raw'
    polyData.point_data.add_array(rawRFP)
    polyData.point_data.get_array(3).name = 'rRFP'
    polyData.point_data.add_array(rawGFP)
    polyData.point_data.get_array(4).name = 'rGFP'
    polyData.point_data.add_array(A)
    polyData.point_data.get_array(5).name = 'bkstRFP'
    polyData.point_data.add_array(B)
    polyData.point_data.get_array(6).name = 'bkstGFP'
    polyData.point_data.update()
    polyData.point_data.add_array(tubeWidth)
    polyData.point_data.get_array(7).name = 'tubeWidth'
    polyData.point_data.update()
# =============================================================================
#           Output
# =============================================================================
    writer = tvtk.PolyDataWriter()
    fileString = os.getcwd()+'\\Norm_'+fileKey+'_skeleton.vtk'
    writer = tvtk.PolyDataWriter()
    writer.set(file_name=fileString)
    writer.set_input(polyData)
    writer.update()
    print"done!"
