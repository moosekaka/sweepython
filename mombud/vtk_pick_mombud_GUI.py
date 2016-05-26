# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 17:54:10 2015
       module to pick mom, bud ,neck positions
"""
import os
import os.path as op
from traits.api import HasTraits, Instance,\
    on_trait_change, Dict, Range, Button, Str, Array, Trait
from traitsui.api import View, Item, Group, HSplit
import pandas as pd
from pandas import DataFrame
import numpy as np
from tvtk.api import tvtk
from tvtk.tvtk_classes.arrow_source import ArrowSource

from mayavi.sources.vtk_data_source import VTKDataSource
from mayavi.sources.api import ParametricSurface
from mayavi import mlab
from mayavi.core.api import Source, Engine, PipelineBase
from mayavi.core.ui.api import SceneEditor, MlabSceneModel, EngineView, \
                               MayaviScene
from seaborn import xkcd_palette as scolor
import wrappers as wr
# pylint: disable=C0103, E1136
datadir = op.join(os.getcwd(), 'mutants')
#xkcd palette colors
colors = ["medium green",
          "light blue",
          "red"]
palette = {col: rgb for col, rgb in zip(colors, scolor(colors))}


def setup_data(fname):
    """Given a VTK file name `fname`, this creates a mayavi2 reader
    for it and adds it to the pipeline.  It returns the reader
    created.
    """
    dat = tvtk.PolyDataReader(file_name=fname)
    dat.update()   # very IMPORTANT!!!
    src = VTKDataSource(data=dat.output)
    return src


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


class ArrowClass(HasTraits):
    """
    Container for VTKSource of arrow glpyh
    """
    def_par = dict(shaft_radius=0.01,
                   shaft_resolution=24,
                   tip_length=.1,
                   tip_radius=.035,
                   tip_resolution=18)
    arrow_src = Instance(ArrowSource, kw=def_par)

    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        self.vtk_src = VTKDataSource(data=self.arrow_src.output)

    def viz_arrow(self, **kwargs):
        obj = mlab.pipeline.surface(self.vtk_src, **kwargs)
        setattr(self, 'surf', obj)

    def transform(self, transform_filter):
        actor = self.surf.actor.actor
        actor.user_transform = transform_filter


class MitoSkel(HasTraits):
    """
    Container for VTKSource of mitoskel
    """
    data_src = Instance(Source)

    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)

    def viz_skel(self, **kwargs):
        tube = mlab.pipeline.tube(self.data_src, tube_sides=32, **kwargs)
        skelsurf = mlab.pipeline.surface(tube, **kwargs)
        setattr(self, 'surf', skelsurf)
        mmgr = skelsurf.module_manager.scalar_lut_manager
        mmgr.show_legend = False
        mmgr.reverse_lut = True
        mmgr.lut_mode = 'RdBu'
        mmgr.number_of_labels = 5
        mmgr.scalar_bar.label_format = '%4.f'
        mmgr.label_text_property.font_size = 12
        mmgr.scalar_bar_representation.position = [.85, .25]
        mmgr.scalar_bar_representation.position2 = [.1, .4]

    def transform(self, transform_filter):
        actor = self.surf.actor.actor
        actor.user_transform = transform_filter


class CellEllipse(HasTraits):
    """
    Ellipse class container for mom bud cells
    """
    name = Str(desc='must be of the form mom_xx or bud_xx')
    data = Dict(desc='dictionary of ellipse parameters')
    dataframe = Instance(DataFrame)

    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        # name splitting needs the name format to be as described above in
        # order to select the right data
        self.data = setup_ellipsedata(self.name.split('_')[0], self.dataframe)
        self.src = getellipsesource(self.data['major'],
                                    self.data['minor'])

    def make_surf(self, **kwargs):
        obj = mlab.pipeline.surface(self.src,
                                    name=self.name,
                                    **kwargs)
        setattr(self, 'surf', obj)

    def adjust_ellipse(self):
        actor = getattr(self, 'surf').actor
        actor.property.opacity = .35
        actor.property.color = (.9, .9, .0)
        actor.mapper.scalar_visibility = False
        actor.property.backface_culling = True
        actor.property.specular = 0.1
        actor.property.frontface_culling = True
        actor.actor.position = np.array([self.data['xc'],
                                         self.data['yc'],
                                         self.data['zpos']])
        actor.actor.orientation = np.array([0, 0, self.data['angle']])

    def update_zpos(self, zpos):
        objactor = getattr(self, 'surf')
        x, y, _ = objactor.actor.actor.position
        objactor.actor.actor.set(position=[x, y, zpos])

    def transform(self, transform_filter):
        actor = self.surf.actor.actor
        actor.user_transform = transform_filter


class MombudPicker(HasTraits):
    """
    A Traits and Mayavi based object for picking moms and bud and setting
    orientation to the cell axis
    """
    name = Str()
    data_src3d = Instance(MitoSkel)
    momellipse = Instance(CellEllipse)
    budellipse = Instance(CellEllipse)

    # The first engine. As default arguments (an empty tuple) are given,
    # traits initializes it.
    engine1 = Instance(Engine, args=())
    scene1 = Instance(MlabSceneModel)  # initiliazes _scene1_default
    scene2 = Instance(MlabSceneModel)  # initiliazes _scene2_default
    arrow_src = Instance(ArrowClass)   # initiliazes _arrow_src_default

    neck_pos = None
    base_pos = None
    tip_pos = None
    picktype = Trait( 'mom', {'mom': 'base',
                              'neck': 'neck',
                              'bud': 'tip'})

    button_save = Button('SaveOutput')
    button_arrow = Button('Arrow')
    button_transform = Button('Transform')
    cursors = Dict({'base': None, 'tip': None, 'neck': None})
    spheres = Dict({'base': None, 'tip': None, 'neck': None})
    # colors defined from xkcd pallete
    cur_col = Dict(
        {part: col for col, part in zip(colors, ['tip', 'base', 'neck'])})

    engine_view = Instance(EngineView)

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
                        Group('_','_',
                             Item('z_position'),
                             Item('picktype', style = 'custom',   label = 'Picker Type')),
                        show_labels=False
                        ),
                  Group(
                       Item('scene2',
                            editor=SceneEditor(), height=600,
                            width=600, show_label=False),
                       '_',
                       'button_save',
                       'button_arrow',
                       'button_transform',
                       show_labels=False,
                       ),
                       )
                )
#
#    view = View(HSplit(
#                 Group(
#                       Item('engine_view',
#                            style='custom',
#                            resizable=True),
#                       show_labels=False
#                  ),
#                  Group(
#                       Item('scene1',
#                            editor=SceneEditor(scene_class=MayaviScene),
#                            height=600,
#                            width=600),
#                       'z_position',
#                       'button1',
#                       'button2',
#                       'button3',
#                       show_labels=False,
#                  ),
#                  Group(
#                       Item('scene2',
#                            editor=SceneEditor(), height=600,
#                            width=600, show_label=False),
#                       '_',
#                       'button4',
#                       'button5',
#                       'transform',
#                       show_labels=False,
#                  ),
#                ),
#                resizable=True,
#               )
    # pylint: enable=C0330
    def __init__(self, label='z_position', **traits):
        # init the parent class HasTraits
        HasTraits.__init__(self, **traits)
        self.arrow_src = ArrowClass()

        # set initial z-position of cell based on vtk bounds, also
        # initialize Range to this position
        zmin, zmax = self.data_src3d.data_src.data.bounds[4:]
        self.momellipse.data['zpos'] = zinit = np.mean((zmin, zmax))
        self.budellipse.data['zpos'] = zinit = np.mean((zmin, zmax))
        trait = Range(zmin, zmax, zinit)
        # note, this adds the Range trait to the object instance!
        self.add_trait(label, trait)

    # ------------------------------------------------------------------------
    # Default values
    # ------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------
    # Traits callback
    # ------------------------------------------------------------------------
    @on_trait_change('scene1.activated')
    def _display_scene1(self):
        self.scene1.mayavi_scene.name = 'scene1'
        self.engine_view = EngineView(engine=self.scene1.engine)
        self.scene1.picker.show_gui = False

        # cursor to mark mom/neck/bud locations
        for key in self.cursors:
            self.cursors[key] = mlab.points3d(0, 0, 0,
                                              mode='2dcross',
                                              scale_factor=.25,
                                              color=palette[self.cur_col[key]],
                                              name=key)
            self.cursors[key].set(visible=False)

        # select which scalartype to show on the skeleton
        self.data_src3d.data_src.point_scalars_name = 'DY_raw'
        self.data_src3d.viz_skel(figure=self.scene1.mayavi_scene)

        # add arrow actor to pipeline
        self.arrow_src.viz_arrow(name='arrow',
                                 opacity=0.35,
                                 figure=self.scene1.mayavi_scene)
        self.arrow_src.surf.set(visible=False)

        # draw and mom/bud ellipse surface and adjust the positions
        self.momellipse.make_surf(figure=self.scene1.mayavi_scene)
        self.momellipse.adjust_ellipse()
        self.budellipse.make_surf(figure=self.scene1.mayavi_scene)
        self.budellipse.adjust_ellipse()

        self._labelscene(self.name, 'scene1')
        self.scene1.mlab.view(0, 0)
        self.scene1.scene.background = (0, 0, 0)
        self.scene1.mayavi_scene.on_mouse_pick(self.pickcb)

    def pickcb(self, obj):
        print self.picktype
        x ,y, z = obj.pick_position
#        self._update_curs(self.picktype)
        self._update_curs(self.picktype_)
        self._drawarrow()
#        if self.tip and self.neck:

#        print "%6.4f, %6.4f, %6.4f" % (x,y,z)

    @on_trait_change('scene2.activated')
    def _display_scene2(self):
        self.scene2.mayavi_scene.name = 'scene2'
        self.scene2.scene.background = (0, 0, 0)

        # instanstiate a copy of mitoskel source data
        vtk_src_copy = VTKDataSource(data=self.data_src3d.data_src.data)
        # We will be transforming this copy!!
        self.src_copy = MitoSkel(data_src=vtk_src_copy)
        self.src_copy.viz_skel(figure=self.scene2.mayavi_scene)

        # instanstiate a copy of ellipses objects
        self.mom_t = CellEllipse(name='mom_t',
                                 dataframe=self.momellipse.dataframe)
        self.bud_t = CellEllipse(name='bud_t',
                                 dataframe=self.budellipse.dataframe)
        self.mom_t.make_surf(figure=self.scene2.mayavi_scene)
        self.mom_t.adjust_ellipse()
        self.mom_t.update_zpos(self.z_position)
        self.bud_t.make_surf(figure=self.scene2.mayavi_scene)
        self.bud_t.adjust_ellipse()
        self.bud_t.update_zpos(self.z_position)

        # initialize sphere glyph markers
        for key in self.spheres:
            center = (0., 0., 0.)
            src = tvtk.SphereSource(center=center, radius=.15)
            self.spheres[key] = mlab.pipeline.surface(
                src.output,
                color=palette[self.cur_col[key]],
                name='%s_trnf' % key)
            self.spheres[key].set(visible=False)

        # Keep the view always pointing up
        self.scene2.scene.interactor.interactor_style = \
            tvtk.InteractorStyleTerrain()
        self._labelscene('Transformed View', 'scene2')
        self.scene2.mlab.view(0, 0)

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
        # if no cursors have been picked, the corres. attribute will be None
#        if getattr(self, '%s_pos' % part) is None:
##            setattr(self, '%s_pos' % part,
##                    Array(value=(0, 0, 0), shape=(3,)))
        setattr(self, '%s_pos' % part,
                self.scene1.picker.pointpicker.pick_position)
        array = getattr(self, '%s_pos' % part)
        self.cursors[part].actor.actor.set(position=array)

        if self.cursors[part].visible is False:
            self.cursors[part].set(visible=True)  # make cursors visible now
        self.spheres[part].actor.actor.position = array  # set pos. for spheres

    def _drawarrow(self):
        if self.base_pos and self.tip_pos and self.neck_pos:
            tr, _, _ = arrowvect(self.base_pos, self.tip_pos, self.neck_pos)
            self.arrow_src.transform(tr)
            self.arrow_src.surf.set(visible=True)

    def _savecsv(self):
        output = op.join(datadir, '%s.csv' % self.name)
        f = open(output, 'w')
        f.write('%s\n' % self.name)
        for part in ['neck', 'base', 'tip']:
            out = getattr(self, '%s_pos', part)
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

    @on_trait_change('z_position')
    def _update_range(self):
        self.momellipse.update_zpos(self.z_position)
        self.budellipse.update_zpos(self.z_position)

    @on_trait_change('button_save')
    def _savecoords(self):
        if self.base_pos and self.tip_pos and self.neck_pos:
            self._savecsv()
            self._save_transform_vtk('trans_obj')
        else:
            print "please finish selecting all three points!"

    @on_trait_change('button_arrow')
    def _redraw_arrow(self):
        if self.base_pos and self.tip_pos and self.neck_pos:
            self._drawarrow()
        else:
            print "please finish selecting all three points!"

    @on_trait_change('button_transform')
    def _draw_transformed(self):
        if self.base_pos and self.tip_pos and self.neck_pos:
            _, rot, _ = arrowvect(self.base_pos, self.tip_pos, self.neck_pos)
            tr_filt = tvtk.Transform()
            rot.transpose()
            tr_filt.translate(np.negative(self.base_pos))
            tr_filt.post_multiply()  # translate, THEN rotate
            tr_filt.concatenate(rot)
            tr_filt.translate([-1, 0, 0])

            # this is the transformed VTK object of interest
            self.trans_obj = tvtk.TransformPolyDataFilter(
                input=self.src_copy.data_src.data,
                transform=tr_filt).output

            self.src_copy.transform(tr_filt)

            # transf mom bud shells
            self.mom_t.transform(tr_filt)
            self.mom_t.update_zpos(self.z_position)
            self.bud_t.transform(tr_filt)
            self.bud_t.update_zpos(self.z_position)

            # transf cursor pts
            for key in self.spheres:
                self.spheres[key].actor.actor.user_transform = tr_filt
                self.spheres[key].set(visible=True)
        else:
            print "please finish selecting all three points!"
        mlab.view(0, 0, figure=self.scene1.mayavi_scene)
        mlab.view(0, 0, figure=self.scene2.mayavi_scene)

##############################################################################
if __name__ == "__main__":
    # Read in the celltracing data into a Pandas Dataframe
    DataSize = pd.read_table(op.join(datadir, 'Results.txt'))
    df = DataSize.ix[:, 1:]
    df['cell'] = df.ix[:, 'Label'].apply(lambda x: x.partition(':')[2])
    counter = df.groupby('cell').Label.count()
    hasbuds = df[df.cell.isin(counter[counter > 1].index.values)]

    mlab.close(all=True)
#    vtkF = wr.swalk(datadir, '*csv', start=0, stop=-4)
    Dcells = {key: None for key in hasbuds.cell.values}
    for i in hasbuds.cell.unique()[20:21]:
        filename = i
        # setup VTK input dataset
        vtkob = setup_data(op.join(datadir,
                                   'normalizedVTK/Norm_%s_skeleton.vtk' %
                                   filename))

        # setup cell ellipse objects
        df2 = getelipspar(filename, df)
        mom = CellEllipse(name='mom', dataframe=df2)
        bud = CellEllipse(name='bud', dataframe=df2)

        mitoskel = MitoSkel(data_src=vtkob)

        # instanstiate the GUI
        m = MombudPicker(name=filename,
                         data_src3d=mitoskel,
                         momellipse=mom,
                         budellipse=bud)
        m.configure_traits()
