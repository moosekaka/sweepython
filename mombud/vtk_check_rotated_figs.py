# -*- coding: utf-8 -*-
"""
Batch visualize skel and surface of rotated cell/ mitoskels
Also writes out to a pickle file the transformed mom, bud and neck coordinates
which is used for mombud analysis
"""
import os
import os.path as op
import cPickle as pickle
from mayavi import mlab
import pandas as pd
import numpy as np
from tvtk.api import tvtk
from wrappers import swalk, UsageError
from mombud.functions import vtkvizfuncs as vz
from mombud.classes.vtk_pick_mombud_class import CellEllipse
# pylint: disable=C0103

datadir = op.join(os.getcwd(), 'mutants', 'transformedData2')
rawdir = op.join(os.getcwd(), 'mutants', 'transformedData2')
surfdir = op.join(os.getcwd(), 'mutants', 'surfaceFiles')
# xkcd palette colors for labels
def_cols = dict(colors=['medium blue', 'bright green', 'red'],
                labels=['base', 'tip', 'neck'])
cur_col, palette = vz.generate_color_labels(**def_cols)

# cell tracing info
DataSize = pd.read_table(op.join(datadir, op.pardir, 'csv', 'Results.txt'))
df_celltracing = DataSize.ix[:, 1:]
df_celltracing['cell'] = \
    df_celltracing.ix[:, 'Label'].apply(lambda x: x.partition(':')[2])
counter = df_celltracing.groupby('cell').Label.count()
hasbuds = \
    df_celltracing[df_celltracing.cell.isin(counter[counter > 1].index.values)]

# vtk data and picked bud, neck, tip inputs
try:
    vtkF = swalk(op.join(rawdir),
                 '*.vtk', start=0, stop=-4)
    vtkS = swalk(surfdir,
                 '*.vtk', start=0, stop=-12)
    mombud = swalk(op.join(datadir), '*.csv', stop=-4)

except UsageError:
    raise

##############################################################################
if __name__ == "__main__":
    WRITE_PICKLE = False  # don't overwrite old pickle file by default
    WRITE_PNG = False

    # Figure to render on
    figone = mlab.figure(size=(800, 600), bgcolor=(.1, .1, .1))
    figone.scene.off_screen_rendering = True

    D = {}  # holder for original bud,neck, tip points
    dfmb = pd.DataFrame(columns=['base', 'neck', 'tip', 'media'])
    for key in vtkF.keys()[:]:
        mlab.clf(figure=figone)  # clear current figure

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

        # setup mom bud shell ellipse
        df_ellipse = vz.getelipspar(key, df_celltracing, useold=False)
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

        # plot the vtk skeleton, surface and ellipse objects
        vz.callreader(vtkF[key])
        vtkobj, _ = vz.cellplot(figone, vtkF[key])
        vact = vz.rendsurf(vtkS[key],
                           color="robin's egg blue")
        vact.actor.actor.user_transform = tr_filt

        # transform original bud, neck and tip points of arrow to be
        # parallel to x-axis unit vector
        for part in D:
            center = D[part]
            src = tvtk.SphereSource(center=D[part], radius=.15,
                                    theta_resolution=32,
                                    phi_resolution=32)
            pt_glyph = mlab.pipeline.surface(src.output,
                                             color=palette[cur_col[part]],
                                             name='%s_trnf' % part,
                                             figure=figone)
            pt_glyph.actor.actor.user_transform = tr_filt
            df[part] = pt_glyph.actor.actor.center
        df['mom'] = df_ellipse.ix['mom', 'vol']
        df['bud'] = df_ellipse.ix['bud', 'vol']
        df['media'] = key.partition("_")[0]
        dfmb = dfmb.append(df)
        mlab.view(0, 0)

        if WRITE_PNG:
            mlab.savefig(op.join(rawdir, '%s.png' % key))

        # dump to pickle file for mom bud analysis
        if WRITE_PICKLE:
            with open(op.join(datadir,
                              'mombudtrans.pkl'), 'wb') as output:
                pickle.dump(dfmb, output)
