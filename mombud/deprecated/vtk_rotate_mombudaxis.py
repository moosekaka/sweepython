# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 00:51:20 2015
module to transform raw vtk to rotated vtk along mom bud axis
@author: sweel
"""
import os
import os.path as op
import cPickle as pickle
import numpy as np
from mayavi import mlab
import pandas as pd
from tvtk.api import tvtk
from mombud.functions import vtkvizfuncs as vz
from mombud.classes.vtk_pick_mombud_class import MitoSkel, CellEllipse
import wrappers as wr
# pylint: disable=C0103

datadir = op.join(os.getcwd(), 'data')
rawdir = op.join(os.getcwd(), 'output')

# xkcd palette colors for labels
def_cols = dict(colors=['medium blue', 'bright green', 'red'],
                labels=['base', 'tip', 'neck'])
cur_col, palette = vz.generate_color_labels(**def_cols)

# filelist and graph list
vtkF = wr.ddwalk(op.join(rawdir, 'normalizedVTK'),
                 '*skeleton.vtk', start=5, stop=-13)

mombud = wr.swalk(op.join(datadir, 'csv'),
                  '*csv', stop=-4)

filekeys = {item: vtkF[media][item] for media
            in sorted(vtkF.keys()) for item
            in sorted(vtkF[media].keys())}

D = {}  # holder for original bud,neck, tip points
dfmb = pd.DataFrame(columns=['base', 'neck', 'tip', 'media'])

# cell tracing info
DataSize = pd.read_table(op.join(datadir, 'csv', 'Results.txt'))
df_celltracing = DataSize.ix[:, 1:]
df_celltracing['cell'] = \
    df_celltracing.ix[:, 'Label'].apply(lambda x: x.partition(':')[2])
counter = df_celltracing.groupby('cell').Label.count()
hasbuds = \
    df_celltracing[df_celltracing.cell.isin(counter[counter > 1].index.values)]

mlab.close(all=True)

##############################################################################
if __name__ == "__main__":
    # Draw cell using cellplot and edgeplot
    figone = mlab.figure(figure='test',
                         size=(800, 600),
                         bgcolor=(0., 0., 0.))
    figone.scene.off_screen_rendering=False

    for key in sorted(mombud.keys())[0::150]:
        mlab.clf(figure=figone)
        print "now on %s" % key

        # get original cursor points
        df_cursorpts = pd.read_csv(mombud[key],
                                   header=0,
                                   names=['x', 'y', 'z'],
                                   index_col=0)
        D['tip'] = np.array(df_cursorpts.ix['tip'])
        D['base'] = np.array(df_cursorpts.ix['base'])
        D['neck'] = np.array(df_cursorpts.ix['neck'])
        # get rotation matrix transform
        t, rot, scale1 = vz.arrowvect(D['base'], D['tip'], D['neck'])
        tr_filt = vz.inverse_tr(rot, D['base'])

        # original vtk skel
        vtkob = vz.setup_vtk_source(op.join(filekeys[key]))
        mitoskel = MitoSkel(data_src=vtkob)
        trans_obj = tvtk.TransformPolyDataFilter(
                input=mitoskel.data_src.data,
                transform=tr_filt).output

        # this is just to visualize the VTK actor
        mitoskel.viz_skel(figure=figone)
        mitoskel.transform(tr_filt)


        df_ellipse = vz.getelipspar(key, df_celltracing, useold=True)
        for mb in ['mom', 'bud']:
            mb_glyph = CellEllipse(name='%s' % mb, dataframe=df_ellipse)
            mb_glyph.make_surf(figure=figone)
            mb_glyph.adjust_ellipse()
            x, y, _ = mb_glyph.surf.actor.actor.position
            mb_glyph.surf.actor.actor.set(
                position=[x, y, df_cursorpts.ix['centerpt', 'x']])
            mb_glyph.surf.actor.actor.user_transform = tr_filt

        # Dataframe to save parameters of transformed object
        df = pd.Series({}, name=key)

        # transform original bud, neck and tip points and writeout
        for part in D:
            center = D[part]
            src = tvtk.SphereSource(center=D[part], radius=.15)
            label = cur_col[part]
            pt_glyph = mlab.pipeline.surface(src.output,
                                             color=palette[label],
                                             name='%s_trnf' % part,
                                             figure=figone)
            pt_glyph.actor.actor.user_transform = tr_filt
            df[part] = pt_glyph.actor.actor.center

        mlab.view(0,0,distance='auto')

## All the transformations objects
#        # ccw 90 rotation and TR to mother bud coord system (for 2nd arrow)
#        ccw90 = np.eye(4)
#        ccw90[0, 0] = 0
#        ccw90[1, 1] = 0
#        ccw90[0, 1] = -1
#        ccw90[1, 0] = 1
#        trans1 = tvtk.Transform()
#        trans1.set_matrix(ccw90.flatten())
#        trans1.scale(1 / 3., 1 / 3., 1 / 3.)
#        trans1.post_multiply()
#        trans1.concatenate(tr)
#
#        # inverse transfrom from mother bud coords to cartesian coord
#        trans2 = tvtk.Transform()
#        rot.transpose()
#        trans2.translate(-base)
#        trans2.post_multiply()  # translate, THEN rotate
#        trans2.concatenate(rot)
#        trans2.translate([-1, 0, 0])
#
#        # transform to scale and translate default arrowsource
#        trans3 = tvtk.Transform()
#        trans3.scale(scale1, scale1, scale1)
#        trans3.post_multiply()
#        trans3.translate([-1, 0, 0])
#
#        # transform for second arrow (rotates 90ccw) at origin
#        trans4 = tvtk.Transform()
#        trans4.scale(scale1 / 3, scale1 / 3, scale1 / 3)
#        trans4.post_multiply()
#        trans4.concatenate(ccw90.flatten())
#        trans4.translate([-1, 0, 0])
#
##       Draw all the transformed data
#        # mother bud axis arrow in mother bud coord system
#        arr_mombud = mlab.pipeline.surface(transformPD.output,
#                                           figure=figone,
#                                           opacity=.33)
#        # second arrow, perpendicular to arr_mombud
#        a2act = mlab.pipeline.surface(arrsource.output,
#                                      figure=figone,
#                                      opacity=.33)
#        a2act.actor.actor.user_transform = trans1
#
#        tippt = tvtk.SphereSource(center=tip, radius=.15)
#        mlab.pipeline.surface(tippt.output,
#                              figure=figone,
#                              color=(.3, 1., .3),
#                              opacity=.33)
#        basept = tvtk.SphereSource(center=base, radius=.15)
#        mlab.pipeline.surface(basept.output,
#                              figure=figone,
#                              color=(.1, .3, 1),
#                              opacity=.33)
#        neckpt = tvtk.SphereSource(center=neck, radius=.15)
#        mlab.pipeline.surface(neckpt.output,
#                              figure=figone,
#                              color=(1, .1, .1),
#                              opacity=.33)
#
#        cell_t = tvtk.TransformPolyDataFilter(input=vtkobj.outputs[0],
#                                              transform=trans2).output
#        mom_t, _ = vz.drawelips('mom', df2, zpos=zp)
#        bud_t, _ = vz.drawelips('bud', df2, zpos=zp)
#        mom_t.actor.actor.user_transform = trans2
#        bud_t.actor.actor.user_transform = trans2
#
#        # transform the arrows and spheres in mombud axis coords back to origin
#        arr_mombud_t = mlab.pipeline.surface(arrsource.output,
#                                             figure=figone,
#                                             opacity=0.33)
#        arr_mombud_t.actor.actor.user_transform = trans3
#        a2act_t = mlab.pipeline.surface(arrsource.output,
#                                        figure=figone,
#                                        opacity=0.33)
#        a2act_t.actor.actor.user_transform = trans4
#        base_t = mlab.pipeline.surface(basept.output,
#                                       figure=figone,
#                                       color=(.1, .3, 1),
#                                       opacity=0.33)
#        tip_t = mlab.pipeline.surface(tippt.output,
#                                      figure=figone,
#                                      opacity=0.33,
#                                      color=(.3, 1., .3))
#        neck_t = mlab.pipeline.surface(neckpt.output,
#                                       figure=figone,
#                                       color=(1, .1, .1),
#                                       opacity=.33)
#        neck_t.actor.actor.user_transform = trans2
#        base_t.actor.actor.user_transform = trans2
#        tip_t.actor.actor.user_transform = trans2
#        dftemp = pd.Series({'base': base_t.actor.actor.center,
#                            'neck': neck_t.actor.actor.center,
#                            'tip': tip_t.actor.actor.center,
#                            'media': key[:3],
#                            'bud': df2.ix['bud', 'vol'],
#                            'mom': df2.ix['mom', 'vol']},
#                           name=key)
##        mlab.close(all=True)
#        dfmb = dfmb.append(dftemp)
#
#        # THIS IS THE TRANSFORMED CELL VTK POLYDATA THAT WE WANT!!
#        cell_t2 = mlab.pipeline.surface(cell_t, figure=figone)
#        cell_t2.actor.mapper.scalar_visibility = True
#        cell_t2.module_manager.lut_data_mode = 'point data'
#        vz.adjustlut(cell_t2)
#
#        t2tube = mlab.pipeline.tube(cell_t2, figure=figone)
#        t2tube.filter.radius = .07
#        t2surfTube = mlab.pipeline.surface(t2tube)
#        t2surfTube.actor.mapper.scalar_visibility = True
#        vz.adjustlut(t2surfTube)
##
##        figone.scene.disable_render = False
##        mlab.show()
##        mlab.view(0, 0, 180)
##        mlab.view(distance='auto')
##         rotated vtk coordinate files
#        w = tvtk.PolyDataWriter(input=cell_t, file_name='%s.vtk' %
#                                op.join(datadir, 'transformedData', key))
#        w.write()
##    with open(op.join(datadir,
##                      'mombudtrans.pkl'), 'wb') as output:
##        pickle.dump(dfmb, output)
