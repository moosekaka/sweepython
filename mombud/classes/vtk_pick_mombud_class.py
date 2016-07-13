# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 17:54:10 2015
       module to pick mom, bud ,neck positions
"""
import os
import os.path as op
from traits.api import HasTraits, Instance,\
    on_trait_change, Dict, Range, Button, Str, Trait
from traitsui.api import View, Item, Group, HSplit
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from tvtk.api import tvtk
from tvtk.tvtk_classes.arrow_source import ArrowSource
from mayavi.sources.vtk_data_source import VTKDataSource
from mayavi import mlab
from mayavi.core.api import Source, Engine
from mayavi.core.ui.api import SceneEditor, MlabSceneModel, EngineView, \
                               MayaviScene
from mombud.functions import vtkvizfuncs as vz
# pylint: disable=C0103, E1136
datadir = op.join(os.getcwd(), 'mutants')


def vtkopen(fpath):
    """
    wrapper to open polydata files
    """
    reader = tvtk.PolyDataReader(file_name=fpath)
    reader.update
    data = reader.output
    return data


def cellpos(vtkdata, base, tip, neck):
    """
    Return DataFrame of cell along mom-bud axis coords.

    Returns
    -------
    celldf : DataFrame
        Columns `DY`, `x`, `wholecell_xaxis`, `type`, `indcell_xaxis`
    """
    data = vtkdata
    npx = data.points.to_array()
    # indices of npx that would sort npx according to the x-axis
    xind = npx[:, 0].argsort()
    dy = data.point_data.get_array(u"DY_minmax").to_array()

    #  individual skeletons xyz and Δψ
    celldf = pd.DataFrame({'x': npx[:, 0][xind],
                           'DY': dy[xind]})
    xn, xb, xt = [x[0] for x in [neck, base, tip]]

    celldf['wholecell_xaxis'] = ((celldf.ix[:, 'x'] -
                                  celldf.ix[0, 'x']) / (xt - xb))

    celldf['type'] = np.where(celldf['x'] > xn, 'bud', 'mom')

    gr = celldf.groupby('type').groups
    for name in gr:
        if name == 'bud':
            celldf.loc[gr[name], 'ind_cell_axis'] = (
                (celldf.loc[gr[name], 'x'] - xn) / (xt - xn))

        else:
            celldf.loc[gr[name], 'ind_cell_axis'] = (
                (celldf.loc[gr[name], 'x'] - xb) / (xn - xb))

    return celldf


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
        """
        Calls mayavi surface modules to visualize the arrow
        """
        obj = mlab.pipeline.surface(self.vtk_src, **kwargs)
        setattr(self, 'surf', obj)

    def transform(self, transform_filter):
        """
        transform the actor using the transform_filter input
        """
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
        """
        Calls mayavi tube and surface modules to visualize the VTK data
        """
        # have to check for allowable kws in kwargs because of quirk in
        # traits factory function, doesn't seem to handle unused kwargs well
        tubekws = ['figure', 'name', 'tube_radius', 'tube_sides']
        kws1 = {k: v for k, v in kwargs.iteritems() if k in tubekws}
        tube = mlab.pipeline.tube(self.data_src, **kws1)

        surfkws = ['color', 'figure', 'name', 'opacity']
        kws2 = {k: v for k, v in kwargs.iteritems() if k in surfkws}
        skelsurf = mlab.pipeline.surface(tube, **kws2)

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
        """
        transform the actor using the transform_filter input
        """
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
        self.data = vz.setup_ellipsedata(self.name.split('_')[0],
                                         self.dataframe)
        self.src = vz.getellipsesource(self.data['major'],
                                       self.data['minor'])

    def make_surf(self, **kwargs):
        """
        calls mlab surface module on obj
        """
        obj = mlab.pipeline.surface(self.src,
                                    name=self.name,
                                    **kwargs)
        setattr(self, 'surf', obj)

    def adjust_ellipse(self):
        """
        sets the ellipse position to the cell tracing positions, formats
        the ellipse actor
        """
        actor = getattr(self, 'surf').actor
        actor.property.opacity = .35
        actor.property.color = vz.rgbcol('puke yellow')
        actor.mapper.scalar_visibility = False
        actor.property.backface_culling = True
        actor.property.specular = 0.1
        actor.property.frontface_culling = True
        actor.actor.position = np.array([self.data['xc'],
                                         self.data['yc'],
                                         self.data['zpos']])
        actor.actor.orientation = np.array([0, 0, self.data['angle']])

    def update_zpos(self, zpos):
        """
        updates the z-position of the actor
        """
        objactor = getattr(self, 'surf')
        x, y, _ = objactor.actor.actor.position
        objactor.actor.actor.set(position=[x, y, zpos])

    def transform(self, transform_filter):
        """
        transform the actor using the transform_filter input
        """
        actor = self.surf.actor.actor
        actor.user_transform = transform_filter


class MombudPicker(HasTraits):
    """
    A Traits and Mayavi based object for picking moms and bud and setting
    orientation to the cell axis
    """
    name = Str
    data_src3d = Instance(MitoSkel)
    momellipse = Instance(CellEllipse)
    budellipse = Instance(CellEllipse)

    # The first engine. As default arguments (an empty tuple) are given,
    # traits initializes it.
    engine1 = Instance(Engine, args=())
    scene1 = Instance(MlabSceneModel)  # initiliazes _scene1_default
    scene2 = Instance(MlabSceneModel)  # initiliazes _scene2_default
    arrow_src = Instance(ArrowClass)   # initiliazes _arrow_src_default

    # xxx_pos are used to hold the picked coordinates
    neck_pos = None
    base_pos = None
    tip_pos = None
    picktype = Trait('mom', {'mom': 'base',
                             'neck': 'neck',
                             'bud': 'tip'})
    cursors = Dict({'base': None, 'tip': None, 'neck': None})
    spheres = Dict({'base': None, 'tip': None, 'neck': None})
    # TraitsUI buttons interface
    button_save = Button('SaveOutput')
    button_transform = Button('Transform')
    button_graph = Button('Graph')

    # default colors for labels
    def_cols = dict(colors=['light blue', 'bright green', 'red'],
                    labels=['base', 'tip', 'neck'])

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
                             width=700),
                        Group('_',
                             Item('z_position'),
                             Item('picktype',
                                  style='custom',
                                  label='Picker Type')),
                        show_labels=False
                        ),
                        Group(
                            Item('scene2',
                                 editor=SceneEditor(), height=600,
                                 width=500, show_label=False),
                            'button_save',
                            'button_transform',
                            'button_graph',
                            show_labels=False,
                            ),
                       ),
                       resizable=True,
                )
    # pylint: enable=C0330

    def __init__(self, label='z_position', **traits):
        # init the parent class HasTraits
        HasTraits.__init__(self, **traits)
        self.arrow_src = ArrowClass()
        self.cur_col, self.palette = vz.generate_color_labels(**self.def_cols)

        # set initial z-position of cell based on vtk bounds, also
        # initialize Range to this position
        zmin, zmax = self.data_src3d.data_src.data.bounds[4:]
        self.momellipse.data['zpos'] = zinit = np.mean((zmin, zmax))
        self.budellipse.data['zpos'] = zinit = np.mean((zmin, zmax))
        trait = Range(zmin, zmax, zinit)

        # This adds the Range trait to the object instance!
        self.add_trait(label, trait)

        #fig aesthetics kws
        self.kws = dict(tube_radius='0.07', tube_sides=32,
                        opacity=1.0, transparent=False)

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
            label = self.cur_col[key]
            self.cursors[key] = mlab.points3d(0, 0, 0,
                                              mode='2dcross',
                                              scale_factor=.25,
                                              color=self.palette[label],
                                              name=key)
            self.cursors[key].set(visible=False)

        # select which scalartype to show on the skeleton
        self.data_src3d.data_src.point_scalars_name = 'DY_raw'
        self.data_src3d.viz_skel(figure=self.scene1.mayavi_scene, **self.kws)

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

        # scene formating
        self._labelscene(self.name, 'scene1')
        self.scene1.scene.background = (0, 0, 0)
        self.scene1.mlab.view(0, 0)
        self.scene1.mayavi_scene.on_mouse_pick(self._pickcb)

    def _pickcb(self, obj):
        self._update_curs(self.picktype_)
        self._drawarrow()

    @on_trait_change('scene2.activated')
    def _display_scene2(self):
        self.scene2.mayavi_scene.name = 'scene2'
        self.scene2.scene.background = (0, 0, 0)

        # instanstiate a copy of mitoskel source data
        vtk_src_copy = VTKDataSource(data=self.data_src3d.data_src.data)

        # We will be transforming this copy!!
        self.src_copy = MitoSkel(data_src=vtk_src_copy)
        self.src_copy.viz_skel(figure=self.scene2.mayavi_scene, **self.kws)

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
            src = tvtk.SphereSource(center=center, radius=.15,
                                    theta_resolution=32,
                                    phi_resolution=32)
            self.spheres[key] = mlab.pipeline.surface(
                src.output,
                color=self.palette[self.cur_col[key]],
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
        Logic to update the cursor points `part` and the corres. sphere glyph
        on the transformed view, using the mouse Spoint picker
        """
        # if no cursors have been picked, the xxx_pos. attribute will be None
        setattr(self, '%s_pos' % part,
                self.scene1.picker.pointpicker.pick_position)
        array = getattr(self, '%s_pos' % part)
        self.cursors[part].actor.actor.set(position=array)
        if not self.cursors[part].visible:
            self.cursors[part].set(visible=True)  # make cursors visible now
        self.spheres[part].actor.actor.position = array  # set pos. for spheres

    def _drawarrow(self):
        # only draws an arrow when all three points are picked
        if self.base_pos and self.tip_pos and self.neck_pos:
            tr, _, _ = vz.arrowvect(self.base_pos,
                                    self.tip_pos,
                                    self.neck_pos)
            self.arrow_src.transform(tr)
            self.arrow_src.surf.set(visible=True)

    def _savecsv(self):
        output = op.join(datadir, '%s.csv' % self.name)
        with open(output, 'w') as f:
            f.write('%s\n' % self.name)
            for part in ['neck', 'base', 'tip']:
                out = getattr(self, '%s_pos' % part)
                f.write('{},{},{},{}\n'.format(part, *tuple(out)))
            f.write('centerpt,{}\n'.format(self.z_position))
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

    @on_trait_change('button_graph')
    def _plotgraph(self):
        if self.base_pos and self.tip_pos and self.neck_pos:
            self.grph = cellpos(self.data_src3d.data_src.data,
                                self.base_pos,
                                self.tip_pos,
                                self.neck_pos)
            _g = self.grph
            _, ax1 = plt.subplots(1, 1)
            with sns.plotting_context('talk'):
                sns.violinplot(x='type',
                               y='DY',
                               hue='type',
                               data=_g,
                               ax=ax1)

    @on_trait_change('button_transform')
    def _draw_transformed(self):
        # allows transformation only when three points picked
        if self.base_pos and self.tip_pos and self.neck_pos:
            _, rot, _ = vz.arrowvect(self.base_pos,
                                     self.tip_pos,
                                     self.neck_pos)

            # this filter transforms the cell-axis to the x-axis unit vector
            tr_filt = vz.inverse_tr(rot, self.base_pos)

            # this is the transformed VTK object of interest
            self.trans_obj = tvtk.TransformPolyDataFilter(
                input=self.src_copy.data_src.data,
                transform=tr_filt).output

            # this is just to visualize the VTK actor
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
        mlab.view(0, 0, distance='auto', figure=self.scene1.mayavi_scene)
        mlab.view(0, 0, distance='auto', figure=self.scene2.mayavi_scene)


def main(*args):
    """
    Mom bud cell axis picking tool

    Args
    ----
    start, stop, step : int
        arguments for slice(), i.e. [start:stop:step] to control the range
        cells to be picked
    """
    # Read in the celltracing data into a Pandas Dataframe
    DataSize = pd.read_csv(op.join(datadir, 'csv', 'Results_combined.csv'))
    df = DataSize.ix[:, 1:]
    df['cell'] = df.ix[:, 'Label'].apply(lambda x: x.partition(':')[2])
    counter = df.groupby('cell').size()
    hasbuds = df[df.cell.isin(counter[counter > 1].index.values)]
    hasbuds = hasbuds[hasbuds.cell.str.contains('071016')]
    for i in hasbuds.cell.unique()[slice(*args)]:
        filename = i
        foldername = filename.partition('_')[0]

        # setup VTK input dataset
        vtkob = vz.setup_vtk_source(
            op.join(datadir,
                    'normalizedVTK',
                    foldername,
                    'Norm_%s_skeleton.vtk' % filename))

        # setup cell ellipse objects
        df2 = vz.getelipspar(filename, df)
        mom = CellEllipse(name='mom', dataframe=df2)
        bud = CellEllipse(name='bud', dataframe=df2)

        # setup VTK skel data
        mitoskel = MitoSkel(data_src=vtkob)

        # instanstiate the GUI
        m = MombudPicker(name=filename,
                         data_src3d=mitoskel,
                         momellipse=mom,
                         budellipse=bud)
        m.configure_traits()

# ===========================================================================
if __name__ == "__main__":
    mlab.close(all=True)
    main(82,85)
