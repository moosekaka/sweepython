# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 17:54:10 2015
       module to pick mom, bud ,neck positions
"""
import os
import os.path as op
from traits.api import HasTraits, Instance,\
    on_trait_change, Dict, Range, Button, Str, Array, PrototypedFrom, Float
from traitsui.api import View, Item, HGroup, Group, HSplit, VSplit
import pandas as pd
import numpy as np
from tvtk.api import tvtk
from tvtk.common import configure_input_data
from mayavi.sources.vtk_data_source import tvtk, VTKDataSource
from mayavi.sources.api import ParametricSurface
from mayavi import mlab
from mayavi.core.api import PipelineBase, Source
from mayavi.core.ui.api import SceneEditor, MlabSceneModel, EngineView

from seaborn import xkcd_palette as scolor
import wrappers as wr
# pylint: disable=C0103
datadir = op.join(os.getcwd(), 'mutants')
#xkcd palette colors
colors = ["medium green",
          "light blue",
          "red"]
palette = {col:rgb for col, rgb in zip(colors, scolor(colors))}

def getelipspar(fname, datafile):
    """ parameters for ellipse from cell tracing """
    dftemp = datafile[datafile.cell == fname]
    dftemp = dftemp.sort_values(by='vol')
    dftemp.index = ['bud', 'mom']
    # reverse the y-coordinate system (VTK vs ImageJ)
    dftemp['center'] = zip((dftemp.X) * .055, (250 - dftemp.Y) * .055)
    return dftemp


def setup_data(fname):
    """Given a VTK file name `fname`, this creates a mayavi2 reader
    for it and adds it to the pipeline.  It returns the reader
    created.
    """
    dat = tvtk.PolyDataReader(file_name=fname)
    dat.update()   # very IMPORTANT!!!
    src = VTKDataSource(data=dat.output)
    return src


def setup_ellipsedata(strg, dF):
    """
    setup dictionary based on pandas dataframe to serve as input data for
    CellEllipse class
    """
    D = {}
    D['major'] = dF.ix[strg, 'Major'] * .055 / 2
    D['minor'] = dF.ix[strg, 'Minor'] * .055 / 2
    D['angle'] = dF.ix[strg, 'Angle']
    D['xc'] = dF.ix[strg, 'center'][0]
    D['yc'] = dF.ix[strg, 'center'][1]
    D['zpos'] = 0
    return D


def getellipsesource(datadic):
    """
    Convenience wrapper to generate a Mayavi Source object based on datadic
    """
    source = ParametricSurface()
    source.function = 'ellipsoid'
    source.parametric_function.set(x_radius=datadic['major'],
                                   y_radius=datadic['minor'],
                                   z_radius=datadic['minor'])
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
    name = Str()
    data = Dict()

    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        self.src = getellipsesource(self.data)


class MombudPicker(HasTraits):
    """
    A Traits and Mayavi based object for picking moms and bud and setting
    orientation to the cell axis
    """
    name = Str()
    data_src3d = Instance(Source)
    scene3d = Instance(MlabSceneModel, (), editor=SceneEditor())
    scene2 = Instance(MlabSceneModel, (), editor=SceneEditor())
    arrow_src = None
    arrow_actor = None
    momellipse = Instance(CellEllipse, ())
    budellipse = Instance(CellEllipse, ())
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
    transform = Button('Tranform')
    cursors = Dict({'base': None, 'tip': None, 'neck': None})
    # colors defined from xkcd pallete
    cur_col = Dict(
        {part: col for col, part in zip(colors, ['tip', 'base', 'neck'])})

    engine_view = Instance(EngineView)

    def __init__(self, label='z_position', **traits):
        # init the parent class HasTraits
        HasTraits.__init__(self, **traits)
        self.scene3d.mayavi_scene.name = self.name
        self.engine_view = EngineView(engine=self.scene3d.engine)
        zmin, zmax = self.data_src3d.outputs[0].bounds[4:]
        self.momellipse.data['zpos'] = zinit = np.mean((zmin, zmax))
        trait = Range(zmin, zmax, zinit)
        self.add_trait(label, trait)

    @on_trait_change('scene3d.activated')
    def _display_scene3d(self):
        self.scene3d.picker.show_gui = False

        # cursor to mark mom/neck/bud locations
        for key in self.cursors:
            self.cursors[key] = mlab.points3d(0, 0, 0,
                                              mode='2dcross',
                                              scale_factor=.25,
                                              color=palette[self.cur_col[key]],
                                              name=key)

        # the scalar value data for the skeleton , adjust the LUT
        tube = mlab.pipeline.tube(self.data_src3d, tube_sides=32,
                                  figure=self.scene3d.mayavi_scene)
        self.data_src3d.point_scalars_name = 'DY_raw'
        surfTube = mlab.pipeline.surface(tube)
        adjustlut(surfTube)

        # draw and mom/bud ellipse surface and adjust the positions
        self.emom = mlab.pipeline.surface(self.momellipse.src,
                                          name='momSurf',
                                          figure=self.scene3d.mayavi_scene)
        self.ebud = mlab.pipeline.surface(self.budellipse.src,
                                          name='budSurf',
                                          figure=self.scene3d.mayavi_scene)
        adjustellipse(self.emom, self.momellipse.data)
        adjustellipse(self.ebud, self.budellipse.data)
        self._update_z()  # update z to default range value

        self._labelscene(self.name, 'scene3d')
        self.scene3d.mlab.view(0, 0, 180)
        self.scene3d.scene.background = (0, 0, 0)

    @on_trait_change('scene2.activated')
    def display_scene2(self):
        mlab.clf(figure=self.scene2.mayavi_scene)
        self._labelscene('Transformed View', 'scene2')
        self.scene2.scene.background = (0, 0, 0)
        # Keep the view always pointing up
        self.scene2.scene.interactor.interactor_style = \
            tvtk.InteractorStyleTerrain()

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
        setattr(self, part, self.scene3d.picker.pointpicker.pick_position)
        array = getattr(self, '%s' % part)
        self.cursors[part].actor.actor.set(position=array)

    def _drawarrow(self):
        # initialize arrow object
        if self.arrow_src is None:
            self.arrow_src = tvtk.ArrowSource(shaft_radius=.01,
                                              shaft_resolution=18,
                                              tip_length=.15,
                                              tip_radius=.05,
                                              tip_resolution=18)
        # remove the previous arrow object if it exists
        if self.arrow_actor is not None:
            self.scene3d.mayavi_scene.scene.remove_actor(self.arrow_actor)

        # tr is the transform filter (transformation matrix object)
        tr, _, _ = arrowvect(self.base, self.tip, self.neck)

        tranf_output = tvtk.TransformPolyDataFilter(
            input=self.arrow_src.output, transform=tr).output
        mapper = tvtk.PolyDataMapper()
        self.arrow_actor = tvtk.Actor()
        mapper.set(input=tranf_output)
        self.arrow_actor.set(mapper=mapper)
        self.scene3d.mayavi_scene.scene.add_actor(self.arrow_actor)

    @on_trait_change('z_position')
    def _update_z(self):
        self.emom.actor.actor.set(position=[self.momellipse.data['xc'],
                                            self.momellipse.data['yc'],
                                            self.z_position])
        self.ebud.actor.actor.set(position=[self.budellipse.data['xc'],
                                            self.budellipse.data['yc'],
                                            self.z_position])

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
        output = op.join(datadir, '%s.csv' % self.name)
        f = open(output, 'w')
        f.write('%s\n' % self.name)
        for part in ['neck', 'base', 'tip']:
            out = getattr(self, part)
            f.write('{},{},{},{}\n'.format(part, *tuple(out)))
        f.write('centerpt,{}\n'.format(self.z_position))
        f.close()
        print 'results recorded for {}!'.format(self.name)

    @on_trait_change('button5')
    def redraw_arrow(self):
        self._drawarrow()

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
                       Item('scene3d',
                            editor=SceneEditor(),
                            height=600,
                            width=600),
                       'button1',
                       'button2',
                       'button3',
                       show_labels=False,
                  ),
                  Group(
                       Item('scene2',
                            editor=SceneEditor(), height=600,
                            width=600, show_label=False),
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
    df['vol'] = 4 / 3 * np.pi * (df.Major * .055 / 2) * (df.Minor * .055 / 2)
    counter = df.groupby('cell').Label.count()
    hasbuds = df[df.cell.isin(counter[counter > 1].index.values)]

    mlab.close(all=True)
#    vtkF = wr.ddwalk(datadir, '*csv', start=0, stop=-4)

    for i in hasbuds.cell.unique()[10:11]:
        filename = i
        vtkob = setup_data(op.join(datadir,
                                   'normalizedVTK/Norm_%s_skeleton.vtk' %
                                   filename))
        zpos_init = np.mean(vtkob.outputs[0].bounds[4:])

        # setup input dataset
        df2 = getelipspar(filename, df)
        Dmom = setup_ellipsedata('mom', df2)
        Dbud = setup_ellipsedata('bud', df2)
        mom = CellEllipse(name='mom', data=Dmom)
        bud = CellEllipse(name='bud', data=Dbud)

        m = MombudPicker(name=filename,
                         data_src3d=vtkob,
                         momellipse=mom,
                         budellipse=bud)

        m.configure_traits()
        m.scene3d.mayavi_scene.scene.reset_zoom()
