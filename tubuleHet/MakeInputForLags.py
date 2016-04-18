# -*- coding: utf-8 -*-
"""
Script to create a distribution of lags for Intensities along edges of cells
@author: sweel
"""
import os
import os.path as op
import cPickle as pickle
from tubuleHet.autoCor.fitDistr import fitDist
from pipeline.make_networkx import makegraph
import wrappers as wr
import vtk
from tvtk.api import tvtk
# =============================================================================
#           Data initialization
# =============================================================================
datadir = op.join(os.getcwd(), 'data')
rawdir = op.join(os.getcwd(), 'output')
vtkF = wr.ddwalk(op.join(rawdir, 'normalizedVTK'),
                 '*skeleton.vtk', start=5, stop=-13)

for mediatype in sorted(vtkF.keys())[:]:
    for filekey in sorted(vtkF[mediatype].keys())[:]:
        lNorm = []
        lNormP = []
        llineId = []
        randNDY = []
        randUDY = []

        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(vtkF[mediatype][filekey])
        reader.Update()
        vtkdata = reader.GetOutput()
        tvtkdata = tvtk.to_tvtk(vtkdata)

        _, _, nxgrph = makegraph(vtkdata, filekey)
        output = fitDist(tvtkdata, nxgrph)

        lineId = output[-1:]
        #    0:4 for scaled, 5:8 for unscaled
        sampNRaw, sampURaw, unScaled, rPermute = output[4:8]
        lineId = output[-1:]
        lNorm.append(unScaled)
        lNormP.append(rPermute)
        llineId.append(lineId)
        randNDY.append(sampNRaw)
        randUDY.append(sampURaw)
        print "append norm dist %s" % filekey

        with open(
         op.join(rawdir, 'fitted_data', "%s.pkl" % filekey), 'wb') as OUT:
            pickle.dump((lNorm, lNormP, randNDY, randUDY, llineId),
                        OUT, protocol=2)
