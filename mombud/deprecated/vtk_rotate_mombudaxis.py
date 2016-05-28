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

# vtk data and cell picking points data
vtkF = wr.ddwalk(op.join(rawdir, 'normalizedVTK'),
                 '*skeleton.vtk', start=5, stop=-13)

mombud = wr.swalk(op.join(datadir, 'csv'),
                  '*csv', stop=-4)

filekeys = {item: vtkF[media][item] for media
            in sorted(vtkF.keys()) for item
            in sorted(vtkF[media].keys())}

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
    D = {}  # holder for original bud,neck, tip points
    dfmb = pd.DataFrame(columns=['base', 'neck', 'tip', 'media'])
    # Draw cell using cellplot and edgeplot
    figone = mlab.figure(figure='test',
                         size=(800, 600),
                         bgcolor=(0., 0., 0.))
    figone.scene.off_screen_rendering = True

    for key in sorted(mombud.keys())[:]:
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
        if 'centerpt'  in df_cursorpts:
            zpos = df_cursorpts.ix['centerpt', 'x']
        else:
            zpos = np.mean(vtkob.data.bounds[4:])

        for mb in ['mom', 'bud']:
            mb_glyph = CellEllipse(name='%s' % mb, dataframe=df_ellipse)
            mb_glyph.make_surf(figure=figone)
            mb_glyph.adjust_ellipse()
            x, y, _ = mb_glyph.surf.actor.actor.position
            mb_glyph.surf.actor.actor.set(
                position=[x, y, zpos])
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
        mlab.view(0, 0)

        df['mom'] = df_ellipse.ix['mom', 'vol']
        df['bud'] = df_ellipse.ix['bud', 'vol']
        df['media'] = key.partition("_")[0]
        dfmb = dfmb.append(df)

        w = tvtk.PolyDataWriter(input=trans_obj, file_name='%s.vtk' %
                                op.join(datadir, 'transformedData', key))

        w.write()
        mlab.savefig(op.join(datadir, 'transformedData', '%s.png' % key))

    with open(op.join(datadir, 'transformedData',
                      'mombudtrans_new.pkl'), 'wb') as output:
        pickle.dump(dfmb, output)
