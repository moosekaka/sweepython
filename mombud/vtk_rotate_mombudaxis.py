# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 00:51:20 2015
	module to transform raw vtk to rotated vtk along mom bud axis

@author: sweel
"""

import os
import os.path as op
import numpy as np
from mayavi import mlab
import fnmatch
import pandas as pd
from tvtk.api import tvtk
from collections import defaultdict
from mombud.vtk_viz import vtkvizfuncs as vz
import cPickle as pickle
# pylint: disable=C0103
# pylint:disable=E1101
vtkF = defaultdict(dict)
mombud = defaultdict(dict)
datadir = op.join(os.getcwd(), 'data')

# =============================================================================
# filelist and graph list
# =============================================================================
for root, dirs, files in os.walk(op.join(datadir, 'normSkel')):
    for i in files:
        if fnmatch.fnmatch(i, '*skeleton.vtk'):
            media = root.rsplit(os.sep, 1)[1]
            vtkF[media][i[5:-13]] = op.join(root, i)
for root, dirs, files in os.walk(op.join(datadir, 'csv')):
    for i in files:
        if fnmatch.fnmatch(i, 'YP*csv'):
            mombud[i[:-4]] = op.join(root, i)

filekeys = {item: vtkF[media][item] for media
            in sorted(vtkF.keys()) for item
            in sorted(vtkF[media].keys())}

DataSize = pd.read_table(op.join(datadir, 'csv', 'Results.txt'))
df = DataSize.ix[:, 1:]
df['cell'] = df.ix[:, 'Label'].apply(lambda x: x.partition(':')[2])
df['vol'] = 4 / 3 * np.pi * (df.Major * .055 / 2) * (df.Minor * .055 / 2) ** 2

# =============================================================================
# Draw cell using cellplot and edgeplot
# =============================================================================

if __name__ == "__main__":
    dfmb = pd.DataFrame(columns=['base', 'neck', 'tip', 'media'])
    mlab.close(all=True)
    for _, key in enumerate(sorted(mombud.keys())[-5:-3]):
        df1 = pd.read_csv(op.join(datadir, 'csv', '%s.csv' % key),
                          header=0,
                          names=['x', 'y', 'z'],
                          index_col=0)
        tip = np.array(df1.ix['tip'])
        base = np.array(df1.ix['base'])
        neck = np.array(df1.ix['neck'])

        filekey = key
        df2 = vz.getelipspar(filekey, df)
        df2 = df2.sort_values('vol')
        df2.reset_index(drop=True, inplace=True)
        df2.index = ['bud', 'mom']
        df2['center'] = zip((df2.X - 25) * .055, (225 - df2.Y) * .055)
        figone = mlab.figure(figure=filekey,
                             size=(800, 600),
                             bgcolor=(0., 0., 0.))

        figone.scene.disable_render = True
        vtkobj, tubeout = vz.cellplot(figone, filekeys[filekey])
        xmin, xmax, ymin, ymax, zmin, zmax = vtkobj.outputs[0].bounds
# ============================================================================
# zposition of center slice
# =============================================================================
        try:
            zp = df1.ix['centerpt'][0]
        except KeyError:
            zp = (zmax - zmin) / 2

        vz.adjustlut(tubeout)
        vz.drawelips('mom', df2, zpos=zp)
        vz.drawelips('bud', df2, zpos=zp)

# =============================================================================
#       get orientation vector defining mom bud axis
# ======================================================================.=======
        tr, rot, scale1 = vz.arrowvect(base, tip, neck)
        arrsource = tvtk.ArrowSource(shaft_radius=.01,
                                     shaft_resolution=18,
                                     tip_length=.15,
                                     tip_radius=.05,
                                     tip_resolution=18)
        transformPD = tvtk.TransformPolyDataFilter()
        transformPD = tvtk.TransformPolyDataFilter(input=arrsource.output,
                                                   transform=tr)
# =============================================================================
# All the transformations objects
# =============================================================================
        # ccw 90 rotation and TR to mother bud coord system (for 2nd arrow)
        ccw90 = np.eye(4)
        ccw90[0, 0] = 0
        ccw90[1, 1] = 0
        ccw90[0, 1] = -1
        ccw90[1, 0] = 1
        trans1 = tvtk.Transform()
        trans1.set_matrix(ccw90.flatten())
        trans1.scale(1 / 3., 1 / 3., 1 / 3.)
        trans1.post_multiply()
        trans1.concatenate(tr)

        # inverse transfrom from mother bud coords to cartesian coord
        trans2 = tvtk.Transform()
        rot.transpose()
        trans2.translate(-base)
        trans2.post_multiply()  # translate, THEN rotate
        trans2.concatenate(rot)
        trans2.translate([-1, 0, 0])

        # transform to scale and translate default arrowsource
        trans3 = tvtk.Transform()
        trans3.scale(scale1, scale1, scale1)
        trans3.post_multiply()
        trans3.translate([-1, 0, 0])

        # transform for second arrow (rotates 90ccw) at origin
        trans4 = tvtk.Transform()
        trans4.scale(scale1 / 3, scale1 / 3, scale1 / 3)
        trans4.post_multiply()
        trans4.concatenate(ccw90.flatten())
        trans4.translate([-1, 0, 0])

# =============================================================================
#       Draw all the transformed data
# =============================================================================
        # mother bud axis arrow in mother bud coord system
        arr_mombud = mlab.pipeline.surface(transformPD.output,
                                           figure=figone,
                                           opacity=.33)
        # second arrow, perpendicular to arr_mombud
        a2act = mlab.pipeline.surface(arrsource.output,
                                      figure=figone,
                                      opacity=.33)
        a2act.actor.actor.user_transform = trans1

        tippt = tvtk.SphereSource(center=tip, radius=.15)
        mlab.pipeline.surface(tippt.output,
                              figure=figone,
                              color=(.3, 1., .3),
                              opacity=.33)
        basept = tvtk.SphereSource(center=base, radius=.15)
        mlab.pipeline.surface(basept.output,
                              figure=figone,
                              color=(.1, .3, 1),
                              opacity=.33)
        neckpt = tvtk.SphereSource(center=neck, radius=.15)
        mlab.pipeline.surface(neckpt.output,
                              figure=figone,
                              color=(1, .1, .1),
                              opacity=.33)

        cell_t = tvtk.TransformPolyDataFilter(input=vtkobj.outputs[0],
                                              transform=trans2).output
        mom_t, _ = vz.drawelips('mom', df2, zpos=zp)
        bud_t, _ = vz.drawelips('bud', df2, zpos=zp)
        mom_t.actor.actor.user_transform = trans2
        bud_t.actor.actor.user_transform = trans2

        # transform the arrows and spheres in mombud axis coords back to origin
        arr_mombud_t = mlab.pipeline.surface(arrsource.output,
                                             figure=figone,
                                             opacity=0.33)
        arr_mombud_t.actor.actor.user_transform = trans3
        a2act_t = mlab.pipeline.surface(arrsource.output,
                                        figure=figone,
                                        opacity=0.33)
        a2act_t.actor.actor.user_transform = trans4
        base_t = mlab.pipeline.surface(basept.output,
                                       figure=figone,
                                       color=(.1, .3, 1),
                                       opacity=0.33)
        tip_t = mlab.pipeline.surface(tippt.output,
                                      figure=figone,
                                      opacity=0.33,
                                      color=(.3, 1., .3))
        neck_t = mlab.pipeline.surface(neckpt.output,
                                       figure=figone,
                                       color=(1, .1, .1),
                                       opacity=.33)
        neck_t.actor.actor.user_transform = trans2
        base_t.actor.actor.user_transform = trans2
        tip_t.actor.actor.user_transform = trans2
        dftemp = pd.Series({'base': base_t.actor.actor.center,
                            'neck': neck_t.actor.actor.center,
                            'tip': tip_t.actor.actor.center,
                            'media': key[:3],
                            'bud': df2.ix['bud', 'vol'],
                            'mom': df2.ix['mom', 'vol']},
                           name=key)
#        mlab.close(all=True)
        dfmb = dfmb.append(dftemp)

        # THIS IS THE TRANSFORMED CELL VTK POLYDATA THAT WE WANT!!
        cell_t2 = mlab.pipeline.surface(cell_t, figure=figone)
        cell_t2.actor.mapper.scalar_visibility = True
        cell_t2.module_manager.lut_data_mode = 'point data'
        vz.adjustlut(cell_t2)

        t2tube = mlab.pipeline.tube(cell_t2, figure=figone)
        t2tube.filter.radius = .07
        t2surfTube = mlab.pipeline.surface(t2tube)
        t2surfTube.actor.mapper.scalar_visibility = True
        vz.adjustlut(t2surfTube)

        figone.scene.disable_render = False
        mlab.view(0, 0, 180)
        mlab.view(distance='auto')
        # rotated vtk coordinate files
#        w = tvtk.PolyDataWriter(input=cell_t, file_name='%s.vtk' % key)
#        w.write()
#    with open(op.join(datadir,
#                      'transformedData',
#                      'mombudtrans.pkl'), 'wb') as output:
#        pickle.dump(dfmb, output)
