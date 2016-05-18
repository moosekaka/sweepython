# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 17:54:10 2015
       module to pick mom, bud ,neck positions
"""
import os
import os.path as op
from traits.api import HasTraits, Instance,\
                       on_trait_change, Dict, Range, Button, Str, Array
from traitsui.api import View, Item, HGroup, Group, HSplit, VSplit, HGroup
import pandas as pd
import numpy as np
from tvtk.api import tvtk
from mayavi.sources.vtk_data_source import tvtk, VTKDataSource
from mayavi.sources.api import ParametricSurface
from mayavi import mlab
from mayavi.core.api import PipelineBase, Source
from mayavi.core.ui.api import SceneEditor, MlabSceneModel, EngineView
import wrappers as wr
## pylint: disable=C0103
datadir = op.join(os.getcwd(), 'mutants')


def adjustlut(vtksurface):
    """ adjust lut colormap
    """
    mmgr = vtksurface.module_manager.scalar_lut_manager
    mmgr.show_legend = True
    mmgr.reverse_lut = True
    mmgr.lut_mode = 'RdBu'
    mmgr.number_of_labels = 5
    mmgr.scalar_bar.label_format = '%4.f'
    mmgr.label_text_property.font_size = 12
    mmgr.scalar_bar_representation.position = [.85, .25]
    mmgr.scalar_bar_representation.position2 = [.1, .4]


def getelipspar(filename, df):
    """ parameters for ellipse from cell tracing """
    dftemp = df[df.cell == filename]
    dftemp.sort_values(by='vol', inplace=True)
    dftemp.index = ['bud', 'mom']
    # reverse the y-coordinate system (VTK vs ImageJ)
    dftemp['center'] = zip((dftemp.X)*.055, (250 - dftemp.Y)*.055)
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


def setup_ellipsedata(strg, df):
    D = {}
    D['major'] = df.ix[strg, 'Major']*.055/2
    D['minor'] = df.ix[strg, 'Minor']*.055/2
    D['angle'] = df2.ix[strg, 'Angle']
    D['xc'] = df2.ix[strg, 'center'][0]
    D['yc'] = df2.ix[strg, 'center'][1]
    D['zpos'] = 0
    return D


def getellipsesource(data):
    #    x, y = df2.ix[strg, 'center']
    #    mlab.points3d(np.array(x), np.array(y), np.array(zpos),
    #                  mode='2dcross',
    #                  scale_factor=.3,
    #                  color=(.5, .7, 0.),)
    source = ParametricSurface()
    source.function = 'ellipsoid'
    source.parametric_function.set(x_radius=data['major'],
                                   y_radius=data['minor'],
                                   z_radius=data['minor'])
    return source


def adjustElip(ee, data):
    actor = ee.actor
    actor.property.opacity = .35
    actor.property.color = (.9, .9, .0)
    actor.mapper.scalar_visibility = False
    actor.property.backface_culling = True
    actor.property.specular = 0.1
    actor.property.frontface_culling = True
    actor.actor.position = np.array([data['xc'], data['yc'], data['zpos']])
    actor.actor.orientation = np.array([0, 0, data['angle']])


def arrowvect(B, A, C):
    """draws a vector based on base, B and tip, A.
    calculates the transformation matrix trans and returns it along with the
    rotation matrix
    """
    normalizedX = np.zeros(3)
    normalizedY = np.zeros(3)
    normalizedZ = np.zeros(3)
    AP = np.zeros(3)
    math = tvtk.Math()
    math.subtract(A, B, normalizedX)  # normalizedX is the arrow unit vector
    math.subtract(C, B, AP)
    length = math.norm(normalizedX)
    math.normalize(normalizedX)
    math.normalize(AP)  # another unit vector used to fix the local x-y plane

    x1, x2, x3 = normalizedX
    t1, t2, t3 = AP
    l3 = -t3/(t1+t2)
    m3 = (t3*x1 - x3*t1 - x3*t2) / (x2*t1 + t2*x2)
    D = np.sqrt((t3 / (t1 + t2))**2 +
                ((t3*x1 - x3*t1 - x3*t2) / (x2*t1 + t2*x2))**2 + 1)
    z1 = l3/D
    z2 = m3/D
    z3 = 1/D
    normalizedZ = np.array([z1, z2, z3])
    math.cross(normalizedZ, normalizedX, normalizedY)
    matrix = tvtk.Matrix4x4()
    matrix.identity()
    for el in range(3):  # rotation matrix to x-axis
        matrix.set_element(el, 0, normalizedX[el])
        matrix.set_element(el, 1, normalizedY[el])
        matrix.set_element(el, 2, normalizedZ[el])
    trans = tvtk.Transform()
    trans.translate(B)  # translate origin to base of arrow
    trans.concatenate(matrix)  # rotate around the base of arrow
    trans.scale(length, length, length)
    return trans, matrix, length


class MombudPicker(HasTraits):
    name = Str()
    data_src3d = Instance(Source)
    data_mom = Dict()
    data_bud = Dict()
    scene3d = Instance(MlabSceneModel, (), editor=SceneEditor())
    momelps_src = Instance(Source)
    budelps_src = Instance(Source)
    arrow_src = Instance(Source)
    z_position = Range(-10., 10., 0)
    neck = Array()
    base = Array()
    tip = Array()
    button1 = Button('Mom')
    button2 = Button('Bud')
    button3 = Button('Neck')
    button4 = Button('SaveOutput')
    button5 = Button('Arrow')

    engine_view = Instance(EngineView)

    def __init__(self, **traits):
        # init the parent class HasTraits
        HasTraits.__init__(self, **traits)

        self.scene3d.mayavi_scene.name = self.name
        self.engine_view = EngineView(engine=self.scene3d.engine)

        self.data_mom['zpos'] = np.mean(
            self.data_src3d.outputs[0].bounds[4:])
        print "%s has init. mean z pos = %6.4f" % (self.name,
                                                   self.data_mom['zpos'])

    @on_trait_change('scene3d.activated')
    def display_scene3d(self):
        self.scene3d.picker.show_gui=False
        x, y, z = self.data_src3d.outputs[0].center
        mlab.points3d(0, 0, 0,
                      mode='2dcross',
                      scale_factor=.25,
                      color=(.9, .1, .1),
                      name='neck')
        mlab.points3d(0, 0, 0,
                      mode='2dcross',
                      scale_factor=.25,
                      color=(.0, .25, .9),
                      name='base')
        mlab.points3d(0, 0, 0,
                      mode='2dcross',
                      scale_factor=.25,
                      color=(.2, .7, .2),
                      name='tip')

        tube = mlab.pipeline.tube(self.data_src3d,
                                  figure=self.scene3d.mayavi_scene)
        emom = mlab.pipeline.surface(self.momelps_src,
                                     name='momSurf',
                                     figure=self.scene3d.mayavi_scene)
        ebud = mlab.pipeline.surface(self.budelps_src,
                                     name='budSurf',
                                     figure=self.scene3d.mayavi_scene)
        adjustElip(emom, self.data_mom)
        adjustElip(ebud, self.data_bud)

        self.data_src3d.point_scalars_name = 'DY_raw'
        surfTube = mlab.pipeline.surface(tube)
        mod_mngr = tube.children[0]
        mmgr = mod_mngr.scalar_lut_manager
        mmgr.scalar_bar.title = 'DY_raw'
        mmgr.data_name = 'DY_raw'
        tube.filter.number_of_sides = 32
        adjustlut(surfTube)
        self.scene3d.mlab.view(0, 0, 180)
        self.scene3d.scene.background = (0, 0, 0)

    @on_trait_change('z_position')
    def update_z(self):
        momsurf = self.scene3d.mayavi_scene.children[-2].children[0].children[0]
        budsurf = self.scene3d.mayavi_scene.children[-1].children[0].children[0]
        momsurf.actor.actor.set(position=[self.data_mom['xc'],
                                          self.data_mom['yc'],
                                          self.z_position])
        budsurf.actor.actor.set(position=[self.data_bud['xc'],
                                          self.data_bud['yc'],
                                          self.z_position])

    @on_trait_change('button1')
    def updatemom(self):
        x, y, z = self.base = self.scene3d.picker.pointpicker.pick_position
        mompick = self.scene3d.mayavi_scene.children[1].children[0].children[0]
        mompick.actor.actor.set(position=[x, y, z])

    @on_trait_change('button2')
    def updatebud(self):
        x, y, z = self.tip = self.scene3d.picker.pointpicker.pick_position
        budpick = self.scene3d.mayavi_scene.children[2].children[0].children[0]
        budpick.actor.actor.set(position=[x, y, z])

    @on_trait_change('button3')
    def updateneck(self):
        x, y, z = self.neck = self.scene3d.picker.pointpicker.pick_position
        neckpick = self.scene3d.mayavi_scene.children[0].children[0].children[0]
        neckpick.actor.actor.set(position=[x, y, z])

    @on_trait_change('button4')
    def savecoords(self):
        output = op.join(datadir, '%s.csv' % self.name)
        f = open(output, 'w')
        f.write('%s\n' % self.name)
        f.write('neck,{},{},{}\n'.format(*tuple(self.neck)))
        f.write('base,{},{},{}\n'.format(*tuple(self.base)))
        f.write('tip,{},{},{}\n'.format(*tuple(self.tip)))
        f.write('centerpt,{}\n'.format(self.z_position))
        f.close()
        print 'results recorded for {}!'.format(self.name)

    # The layout of the dialog created
    view = View(HSplit(Item(name='engine_view',
                            style='custom',
                            resizable=True),
                       VSplit(Item('scene3d',
                              editor=SceneEditor(),
                              height=800,
                              width=600),
                              Group('_', 'z_position',
                                    'button1',
                                    'button2',
                                    'button3',
                                    'button4',
                                    show_labels=False))), resizable=True)

if __name__ == "__main__":
    DataSize = pd.read_table(op.join(datadir, 'Results.txt'))
    df = DataSize.ix[:, 1:]
    df['cell'] = df.ix[:, 'Label'].apply(lambda x: x.partition(':')[2])
    df['vol'] = 4 / 3 * np.pi * (df.Major * .055 / 2) * (df.Minor * .055 / 2)
    counter = df.groupby('cell').Label.count()
    hasbuds = df[df.cell.isin(counter[counter>1].index.values)]

    mlab.close(all=True)
#    vtkF = wr.ddwalk(datadir, '*csv', start=0, stop=-4)

    for i in hasbuds.cell.unique()[10:11]:
        filename = i
        vtkob = setup_data(op.join(datadir,
                                   'normalizedVTK/Norm_%s_skeleton.vtk' %
                                   filename))
        zpos_init = np.mean(vtkob.outputs[0].bounds[4:])

        df2 = getelipspar(filename, df)
        Dmom = setup_ellipsedata('mom', df2)
    #    Dmom['zpos'] = zpos_init
        Dbud = setup_ellipsedata('bud', df2)
    #    Dbud['zpos'] = zpos_init
        m = MombudPicker(name=filename,
                         data_src3d=vtkob,
                         data_mom=Dmom,
                         data_bud=Dbud,
                         momelps_src=getellipsesource(Dmom),
                         budelps_src=getellipsesource(Dbud))
        m.configure_traits()
        m.scene3d.mayavi_scene.scene.reset_zoom()
