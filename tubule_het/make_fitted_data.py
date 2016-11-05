# -*- coding: utf-8 -*-
"""
Script to create a distribution of lags for Intensities along edges of cells
@author: sweel
"""
import os
import os.path as op
import cPickle as pickle
from tubule_het.autoCor.fitDistr import fitDist
from pipeline.make_networkx import makegraph
import wrappers as wr
import vtk
from tvtk.api import tvtk
# =============================================================================
#           Data initialization
# =============================================================================
rawdir = op.join(os.getcwd(), 'old_w_new')
vtkF = wr.ddwalk(op.join(rawdir, 'normalizedVTK'),
                 '*skeleton.vtk', start=5, stop=-13)

for mediatype in sorted(vtkF.keys())[:]:
    for filekey in sorted(vtkF[mediatype].keys())[:]:
        lNorm = []
        lNormP = []
        llineId = []
        randNDY = []
        randUDY = []

        #scaled versions
        lNorm_s = []
        lNormP_s = []
        randNDY_s = []
        randUDY_s = []

        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(vtkF[mediatype][filekey])
        reader.Update()
        vtkdata = reader.GetOutput()
        tvtkdata = tvtk.to_tvtk(vtkdata)

        _, _, nxgrph = makegraph(vtkdata, filekey)
        output = fitDist(tvtkdata, nxgrph)

        lineId = output[-1:]
        #    0:4 for scaled, 5:8 for unscaled
        sampN, sampU, Scaled, sPermute = output[0:4]
        sampNRaw, sampURaw, unScaled, rPermute = output[4:8]
        lineId = output[-1:]
        llineId.append(lineId)

        lNorm.append(unScaled)
        lNormP.append(rPermute)
        randNDY.append(sampNRaw)
        randUDY.append(sampURaw)

        lNorm_s.append(Scaled)
        lNormP_s.append(sPermute)
        randNDY_s.append(sampN)
        randUDY_s.append(sampU)

        print "append norm dist %s" % filekey

        with open(
         op.join(rawdir, 'fitted_data', "%s.pkl" % filekey), 'wb') as OUT:
            pickle.dump((lNorm, lNormP, randNDY, randUDY, llineId),
                        OUT, protocol=2)
        with open(
         op.join(rawdir,
                 'fitted_data_scaled', "%s.pkl" % filekey), 'wb') as OT:
            pickle.dump((lNorm_s, lNormP_s, randNDY_s, randUDY_s, llineId),
                        OT, protocol=2)
