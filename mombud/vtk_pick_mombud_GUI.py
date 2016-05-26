# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 17:54:10 2015
       module to pick mom, bud ,neck positions
"""
import os
import os.path as op
from traits.api import HasTraits, Instance,\
    on_trait_change, Dict, Range, Button, Str, Array
from traitsui.api import View, Item, Group, HSplit
import pandas as pd
import numpy as np
from tvtk.api import tvtk
from mayavi.sources.vtk_data_source import VTKDataSource
from mayavi.sources.api import ParametricSurface
from mayavi import mlab
from mayavi.core.api import PipelineBase, Source, Engine
from mayavi.core.ui.api import SceneEditor, MlabSceneModel, EngineView, \
                               MayaviScene

from seaborn import xkcd_palette as scolor
import wrappers as wr
from pandas import DataFrame
# pylint: disable=C0103
datadir = op.join(os.getcwd(), 'mutants')
#xkcd palette colors
colors = ["medium green",
          "bright blue",
          "red"]
palette = {col:rgb for col, rgb in zip(colors, scolor(colors))}


def getelipspar(fname, dfin):
    """
    Return parameters for ellipse from cell tracing dataframe

    Parameters:
    -----------
    fname: Str
        cell name id, e.g. `MFB1_032016_002_RFPstack_000`
    dfin: DataFrame
        Dataframe file from pandas.csv read of a cell tracing data

    Returns:
    --------
    dftemp
    """
    # selection returns a view by default, we want a copy!
    dfout = dfin[dfin.cell == fname].copy()
    dfout['Major'] = dfout.Major * .055/2
    dfout['Minor'] = dfout.Minor * .055/2
    dfout['vol'] =  \
        4 / 3 * np.pi * dfout.Major * dfout.Minor
    dfout = dfout.sort_values(by='vol')
    dfout.index = ['bud', 'mom']
    # reverse the y-coordinate system (VTK vs ImageJ)
    dfout['center'] = zip((dfout.X) * .055, (250 - dfout.Y) * .055)
    return dfout


def setup_data(fname):
    """Given a VTK file name `fname`, this creates a mayavi2 reader
    for it and adds it to the pipeline.  It returns the reader
    created.
    """
    dat = tvtk.PolyDataReader(file_name=fname)
    dat.update()   # very IMPORTANT!!!
    src = VTKDataSource(data=dat.output)
    return src


def setup_ellipsedata(mom_bud, dF):
    """
    Returns a dict, `D` of ellipse parameters of `cell_ID` from a
    dataframe `dF`.
    """
    D = {}
    D['major'] = dF.ix[mom_bud, 'Major']
    D['minor'] = dF.ix[mom_bud, 'Minor']
    D['angle'] = dF.ix[mom_bud, 'Angle']
    D['xc'] = dF.ix[mom_bud, 'center'][0]
    D['yc'] = dF.ix[mom_bud, 'center'][1]
    D['zpos'] = 0
    return D


def getellipsesource(major, minor):
    """
    Convenience wrapper to generate a Mayavi Source object based on
    `major` and `minor` radius of ellipse parameters
    """
    source = ParametricSurface()
    source.function = 'ellipsoid'
    source.parametric_function.set(x_radius=major,
                                   y_radius=minor,
                                   z_radius=minor)
    return source


def arrowvect(base, tip, neck):
    """
    Draws a vector based on base (mom end) and tip (bud end) of cell.
    Calculates the transformation matrix trans and returns it along with the
    rotation matrix.

    Parameters
    ----------
    base, tip, neck : float array with shape (3L,)
        coordinates for the apical ends of mom and bud cells

    Returns
    -------
    trans : vtkTransform
        Transform filter

    matrix : vtkMatrix4x4
        Rotation matrix

    length : float
      scale factor for arrow length
    """
    normalizedX = np.zeros(3)
    normalizedY = np.zeros(3)
    normalizedZ = np.zeros(3)
    AP = np.zeros(3)
    math = tvtk.Math()
    math.subtract(tip, base, normalizedX)  # normalizedX = arrow unit vector
    math.subtract(neck, base, AP)
    length = math.norm(normalizedX)
    math.normalize(normalizedX)
    math.normalize(AP)  # another unit vector used to fix the local x-y plane

    x1, x2, x3 = normalizedX
    t1, t2, t3 = AP
    l3 = -t3 / (t1 + t2)
    m3 = (t3 * x1 - x3 * t1 - x3 * t2) / (x2 * t1 + t2 * x2)
    D = np.sqrt((t3 / (t1 + t2))**2 +
                ((t3 * x1 - x3 * t1 - x3 * t2) / (x2 * t1 + t2 * x2))**2 + 1)
    z1 = l3 / D
    z2 = m3 / D
    z3 = 1 / D
    normalizedZ = np.array([z1, z2, z3])
    math.cross(normalizedZ, normalizedX, normalizedY)
    matrix = tvtk.Matrix4x4()
    matrix.identity()
    for el in range(3):  # rotation matrix to x-axis
        matrix.set_element(el, 0, normalizedX[el])
        matrix.set_element(el, 1, normalizedY[el])
        matrix.set_element(el, 2, normalizedZ[el])
    trans = tvtk.Transform()
    trans.translate(base)  # translate origin to base of arrow
    trans.concatenate(matrix)  # rotate around the base of arrow
    trans.scale(length, length, length)
    return trans, matrix, length


def adjustellipse(surf, data):
    """
    Adjust ellipse position
    """
    actor = surf.actor
    actor.property.opacity = .35
    actor.property.color = (.9, .9, .0)
    actor.mapper.scalar_visibility = False
    actor.property.backface_culling = True
    actor.property.specular = 0.1
    actor.property.frontface_culling = True
    actor.actor.position = np.array([data['xc'],
                                     data['yc'],
                                     data['zpos']])
    actor.actor.orientation = np.array([0, 0, data['angle']])


def adjustlut(surf):
    """
    Adjust lut colormap of skeleton
    """
    mmgr = surf.module_manager.scalar_lut_manager
    mmgr.show_legend = False
    mmgr.reverse_lut = True
    mmgr.lut_mode = 'RdBu'
    mmgr.number_of_labels = 5
    mmgr.scalar_bar.label_format = '%4.f'
    mmgr.label_text_property.font_size = 12
    mmgr.scalar_bar_representation.position = [.85, .25]
    mmgr.scalar_bar_representation.position2 = [.1, .4]


class CellEllipse(HasTraits):
    """
    Ellipse class container for mom bud cells
    """
    name = Str(desc='either mom or bud')
    data = Dict(desc='dictionary of ellipse parameters')
    dataframe = Instance(DataFrame)

    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        self.data = setup_ellipsedata(self.name, self.dataframe)
        self.src = getellipsesource(self.data['major'],
                                    self.data['minor'])


class MombudPicker(HasTraits):
    """
    A Traits and Mayavi based object for picking moms and bud and setting
    orientation to the cell axis
    """
    name = Str()
    data_src3d = Instance(Source)
    arrow_src = Instance(Source)
    momellipse = Instance(CellEllipse)
    budellipse = Instance(CellEllipse)

    scene1 = Instance(MlabSceneModel, args=())
    scene2 = Instance(MlabSceneModel, args=())
    arrow_actor = Instance(PipelineBase)
    emom = Instance(PipelineBase)
    ebud = Instance(PipelineBase)
    neck = None
    base = None
    tip = None

    button1 = Button('Mom')
    button2 = Button('Bud')
    button3 = Button('Neck')
    button4 = Button('SaveOutput')
    button5 = Button('Arrow')
    transform = Button('Transform')
    cursors = Dict({'base': None, 'tip': None, 'neck': None})
    # colors defined from xkcd pallete
    cur_col = Dict(
        {part: col for col, part in zip(colors, ['tip', 'base', 'neck'])})

    engine_view = Instance(EngineView)

    engine1 = Instance(Engine, args=())

    def __init__(self, label='z_position', **traits):
        # init the parent class HasTraits
        HasTraits.__init__(self, **traits)
        self.scene1.mayavi_scene.name = 'scene1'
        self.engine_view = EngineView(engine=self.scene1.engine)
        self.scene2.mayavi_scene.name = 'scene2'
        zmin, zmax = self.data_src3d.outputs[0].bounds[4:]
        self.momellipse.data['zpos'] = zinit = np.mean((zmin, zmax))
        self.budellipse.data['zpos'] = zinit = np.mean((zmin, zmax))
        trait = Range(zmin, zmax, zinit)
        self.add_trait(label, trait)

    def _scene1_default(self):
        " The default initializer for 'scene1' "
        self.engine1.start()
        scene1 = MlabSceneModel(engine=self.engine1)
        return scene1

    def _scene2_default(self):
        " The default initializer for 'scene1' "
        self.engine1.start()
        scene2 = MlabSceneModel(engine=self.engine1)
        return scene2

    @on_trait_change('scene1.activated')
    def _display_scene1(self):
        self.scene1.picker.show_gui = False

        # cursor to mark mom/neck/bud locations
        for key in self.cursors:
            self.cursors[key] = mlab.points3d(0, 0, 0,
                                              mode='2dcross',
                                              scale_factor=.25,
                                              color=palette[self.cur_col[key]],
                                              name=key)

        # select which scalartype to show on the skeleton
        self.data_src3d.point_scalars_name = 'DY_raw'

        self._tubify('data_src3d', 'scene1')
        self.arrow_actor = mlab.pipeline.surface(self.arrow_src,
                                                 name='arrow',
                                                 opacity=0.5,
                                                 figure=self.
                                                 scene1.mayavi_scene)

        # draw and mom/bud ellipse surface and adjust the positions
        for key in ['mom', 'bud']:
            ellipse = getattr(self, '%sellipse' % key)
            obj = mlab.pipeline.surface(ellipse.src,
                                        name='%sSurf' % key,
                                        figure=self.scene1.mayavi_scene)
            setattr(self, 'e%s' % key, obj)
            adjustellipse(obj, ellipse.data)
            self._update_zpos(obj)

        self._labelscene(self.name, 'scene1')
        self.scene1.mlab.view(0, 0, 180)
        self.scene1.scene.background = (0, 0, 0)

    @on_trait_change('scene2.activated')
    def _display_scene2(self):
        mlab.clf(figure=self.scene2.mayavi_scene)
        self._labelscene('Transformed View', 'scene2')
        self.scene2.scene.background = (0, 0, 0)
        # Keep the view always pointing up
        self.scene2.scene.interactor.interactor_style = \
            tvtk.InteractorStyleTerrain()

    def _tubify(self, data, scene_name):
        datasrc = getattr(self, data)
        scene = getattr(self, scene_name)
        tube = mlab.pipeline.tube(datasrc, tube_sides=32,
                                  figure=scene.mayavi_scene)
        surfTube = mlab.pipeline.surface(tube)
        adjustlut(surfTube)

    def _labelscene(self, label, scene_name):
        scene = getattr(self, scene_name)
        # label text and adjust view
        vtext = tvtk.TextActor()
        vtext.set(input=label, text_scale_mode='viewport')
        vtext.text_property.set(font_size=14)
        scene.mayavi_scene.scene.add_actor(vtext)

    def _update_curs(self, part):
        """
        Logic to update the cursor points `part` based on mayavi point picker
        """
        if hasattr(self, 'part') is not None:
            setattr(self, '%s' % part,
                    Array(value=(0, 0, 0),
                          shape=(3,)))
        setattr(self, part, self.scene1.picker.pointpicker.pick_position)
        array = getattr(self, '%s' % part)
        self.cursors[part].actor.actor.set(position=array)

    def _drawarrow(self):
        tr, _, _ = arrowvect(self.base, self.tip, self.neck)
        self.arrow_actor.actor.actor.user_transform = tr

    def _savecsv(self):
        output = op.join(datadir, '%s.csv' % self.name)
        f = open(output, 'w')
        f.write('%s\n' % self.name)
        for part in ['neck', 'base', 'tip']:
            out = getattr(self, part)
            f.write('{},{},{},{}\n'.format(part, *tuple(out)))
        f.write('centerpt,{}\n'.format(self.z_position))
        f.close()
        print('coordinates for base, tip and neck recorded for {}!'
              .format(self.name))

    def _save_transform_vtk(self, vtkobj):
        data = getattr(self, vtkobj)
        w = tvtk.PolyDataWriter(
            input=data,
            file_name=op.join(datadir, '%s.vtk' % self.name))
        w.write()
        print 'transformed vtk saved as {}.vtk!'.format(self.name)

    def _update_zpos(self, obj):
        if isinstance(obj, basestring):
            objactor = getattr(self, '%s' % obj)
        else:
            objactor = obj
        x, y, _ = objactor.actor.actor.position
        objactor.actor.actor.set(position=[x, y, self.z_position])

    @on_trait_change('z_position')
    def _update_range(self):
        self._update_zpos('emom')
        self._update_zpos('ebud')

    @on_trait_change('button1')
    def _updatemom(self):
        self._update_curs('base')
        if self.tip and self.neck:
            self._drawarrow()

    @on_trait_change('button2')
    def _updatebud(self):
        self._update_curs('tip')
        if self.base and self.neck:
            self._drawarrow()

    @on_trait_change('button3')
    def _updateneck(self):
        self._update_curs('neck')

    @on_trait_change('button4')
    def _savecoords(self):
        if self.base and self.tip and self.neck:
            self._savecsv()
            self._save_transform_vtk('trans_obj')
        else:
            print "please finish selecting all three points!"

    @on_trait_change('button5')
    def _redraw_arrow(self):
        if self.base and self.tip and self.neck:
            self._drawarrow()
        else:
            print "please finish selecting all three points!"

    @on_trait_change('transform')
    def _draw_transformed(self):
        mlab.view(0, 0, figure=self.scene1.mayavi_scene)
        self.scene1.mayavi_scene.scene.reset_zoom()
        if self.base and self.tip and self.neck:
            _, rot, scale1 = arrowvect(self.base, self.tip, self.neck)
            tr_filt = tvtk.Transform()
            rot.transpose()
            tr_filt.translate(np.negative(self.base))
            tr_filt.post_multiply()  # translate, THEN rotate
            tr_filt.concatenate(rot)
            tr_filt.translate([-1, 0, 0])
            skel_polydata = self.data_src3d.outputs[0]
            # this is the transformed VTK object of interest
            self.trans_obj = tvtk.TransformPolyDataFilter(
                input=skel_polydata, transform=tr_filt).output
            mlab.clf(figure=self.scene2.mayavi_scene)
            self._tubify('trans_obj', 'scene2')

            # transf mom bud shells
            for key in ['mom', 'bud']:
                obj = getattr(self, 'e%s' % key)
                ellipse = getattr(self, '%sellipse' % key)
                temp = VTKDataSource(data=mlab.pipeline.get_vtk_src(obj)[0])
                surf = mlab.pipeline.surface(temp,
                                             name='%s_trnf' % key,
                                             figure=self.scene2.mayavi_scene)
                adjustellipse(surf, ellipse.data)
                self._update_zpos(surf)
                surf.actor.actor.user_transform = tr_filt

            # transf cursor pts
            for key in self.cursors:
                center = getattr(self, key)
                src = tvtk.SphereSource(center=center, radius=.15)
                surf = mlab.pipeline.surface(src.output,
                                             color=palette[self.cur_col[key]],
                                             name='%s_trnf' % key)
                surf.actor.actor.user_transform = tr_filt
            self.scene2.mayavi_scene.scene.reset_zoom()
        else:
            print "please finish selecting all three points!"

    # GUI layout
    # pylint: disable=C0330
    view = View(HSplit(
                 Group(
                       Item('engine_view',
                            style='custom',
                            resizable=True),
                       show_labels=False
                  ),
                  Group(
                       Item('scene1',
                            editor=SceneEditor(scene_class=MayaviScene),
                            height=600,
                            width=600),
                       'z_position',
                       'button1',
                       'button2',
                       'button3',
                       show_labels=False,
                  ),
                  Group(
                       Item('scene2',
                            editor=SceneEditor(), height=600,
                            width=600, show_label=False),
                       '_',
                       'button4',
                       'button5',
                       'transform',
                       show_labels=False,
                  ),
                ),
                resizable=True,
               )
    # pylint: enable=C0330
##############################################################################
if __name__ == "__main__":
    DataSize = pd.read_table(op.join(datadir, 'Results.txt'))
    df = DataSize.ix[:, 1:]
    df['cell'] = df.ix[:, 'Label'].apply(lambda x: x.partition(':')[2])
    counter = df.groupby('cell').Label.count()
    hasbuds = df[df.cell.isin(counter[counter > 1].index.values)]

    mlab.close(all=True)
#    vtkF = wr.swalk(datadir, '*csv', start=0, stop=-4)
    Dcells = {key:None for key in hasbuds.cell.values}
    for i in hasbuds.cell.unique()[0:1]:
        filename = i
        # setup VTK input dataset
        vtkob = setup_data(op.join(datadir,
                                   'normalizedVTK/Norm_%s_skeleton.vtk' %
                                   filename))

        # setup cell ellipse objects
        df2 = getelipspar(filename, df)
        mom = CellEllipse(name='mom', dataframe=df2)
        bud = CellEllipse(name='bud', dataframe=df2)
        arrow = VTKDataSource(
            data=tvtk.ArrowSource(shaft_radius=.01,
                                  shaft_resolution=18,
                                  tip_length=.15,
                                  tip_radius=.05,
                                  tip_resolution=18).output
                              )
        m = MombudPicker(name=filename,
                         data_src3d=vtkob,
                         momellipse=mom,
                         budellipse=bud,
                         arrow_src=arrow,
                         )

        m.configure_traits()
        m.scene1.mayavi_scene.scene.reset_zoom()
