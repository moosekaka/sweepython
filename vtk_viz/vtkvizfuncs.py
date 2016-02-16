# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 17:54:10 2015
       functions to render vtk output from pipeline
"""
import numpy as np
from mayavi import mlab
from mayavi.sources.api import ParametricSurface
from tvtk.api import tvtk
import config
# pylint: disable=C0103


def countFuncCall():
    """counter function
    """
    config.counter += 1
    return config.counter


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


def cellplot(fig, filename, scalartype='DY_raw', **kwargs):
    """
    Draw vtk one whole cell

    Parameters
    ----------
    fig : str
        Name of cell

    filelist : list
        list of vtk files to iterate over and draw

    filename : str
        filename

    scalartype : str
        point data type to plot on skeleton, can be of type:
        DY_minmax,
        WidthEq,
        DY_raw,
        rRFP,
        rGFP,
        bkstRFP,
        bkstGFP,

    Returns
    -------
    src : vtk object
        mayavi handle to the vtk obj
    surfTube : mayavi surface object
        mayavi pipeline surface
    """
    rad = kwargs.pop('rad', .05)

    src = mlab.pipeline.open(filename)
    tube = mlab.pipeline.tube(src, figure=fig)
    tube.filter.radius = rad
    surfTube = mlab.pipeline.surface(tube)
    surfTube.actor.mapper.scalar_visibility = True
    mod_mngr = tube.children[0]
    mmgr = mod_mngr.scalar_lut_manager
    mmgr.scalar_bar.title = scalartype
    mmgr.data_name = scalartype
    src.point_scalars_name = scalartype
    # mmgr.show_legend = True
    mmgr.reverse_lut = True
    mmgr.lut_mode = 'RdBu'
    mmgr.number_of_labels = 5
    mmgr.scalar_bar.label_format = '%4.f'
    mmgr.label_text_property.font_size = 12
    mmgr.scalar_bar_representation.position = [.85, .25]
    mmgr.scalar_bar_representation.position2 = [.1, .4]
    tube.filter.number_of_sides = 32
    mlab.view(0, 0, 180)
    mlab.view(distance='auto')
    return src, surfTube



def rendsurf(vtksrc, **kwargs):
    """   add surface.vtk file to pipeline"""

    color = kwargs.pop('color', (0.9, .8, .1))
    alpha = kwargs.pop('alpha', .15)

    src2 = mlab.pipeline.open(vtksrc)
    surface1 = mlab.pipeline.surface(src2)
    surface1.actor.property.opacity = alpha
    surface1.actor.mapper.scalar_visibility = False
    surface1.actor.property.color = color


def labelbpoints(graph, **kwargs):
    """  brnch points and end points
    """
    bcol = kwargs.pop('bcol', (1, .2, 1.0))
    ecol = kwargs.pop('ecol', (.1, 1, 1.0))
    size = kwargs.pop('size', .1)

    Nodes = [nattr['coord'] for _, nattr
             in graph.nodes(data=True)
             if nattr['degree'] > 1]
    Ends = [nattr['coord'] for _, nattr
            in graph.nodes(data=True)
            if nattr['degree'] == 1]

    xyz = np.array(Nodes)
    xyz2 = np.array(Ends)
    points = mlab.pipeline.scalar_scatter(  # branchpoints
        xyz[:, 0], xyz[:, 1], xyz[:, 2])
    points2 = mlab.pipeline.scalar_scatter(  # end points
        xyz2[:, 0], xyz2[:, 1], xyz2[:, 2])
    bpts = mlab.pipeline.glyph(points)
    bpts.glyph.glyph_source.glyph_source.radius = size
    epts = mlab.pipeline.glyph(points2)
    epts.glyph.glyph_source.glyph_source.radius = size

    bpts.actor.property.color = bcol
    bpts.actor.mapper.scalar_visibility = 0
    bpts.actor.property.opacity = 0.8
    bpts.glyph.glyph_source.glyph_source.phi_resolution = 16
    bpts.glyph.glyph_source.glyph_source.theta_resolution = 16

    epts.actor.property.color = ecol
    epts.actor.mapper.scalar_visibility = 0
    epts.actor.property.opacity = 0.1
    epts.glyph.glyph_source.glyph_source.phi_resolution = 16
    epts.glyph.glyph_source.glyph_source.theta_resolution = 16


def getelipspar(filename, df):
    """ parameters for ellipse from cell tracing """
    dftemp = df[df.cell == filename]
    dftemp.reset_index(1, drop=True, inplace=True)
    return dftemp


def drawelips(strg, df2, zpos=0):
    """draw ellipsoid based on getelipspar
    """
    x, y = df2.ix[strg, 'center']
    mlab.points3d(np.array(x), np.array(y), np.array(zpos),
                  mode='2dcross',
                  scale_factor=.3,
                  color=(.5, .7, 0.),)
    source = ParametricSurface()
    source.function = 'ellipsoid'
    source.parametric_function.set(x_radius=df2.ix[strg, 'Major']*.055/2,
                                   y_radius=df2.ix[strg, 'Minor']*.055/2,
                                   z_radius=df2.ix[strg, 'Minor']*.055/2)
    ee = mlab.pipeline.surface(source)
    actor = ee.actor
    actor.property.opacity = .35
    actor.property.color = (.9, .9, .0)
    actor.mapper.scalar_visibility = False
    actor.property.backface_culling = True
    actor.property.specular = 0.1
    actor.property.frontface_culling = True
    actor.actor.position = np.array([x, y, zpos])
    actor.actor.orientation = np.array([0, 0, df2.ix[strg, 'Angle']])
    return ee, source


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
    m3 = (t3*x1-x3*t1-x3*t2)/(x2*t1+t2*x2)
    D = np.sqrt((t3/(t1 + t2))**2 +
                ((t3*x1 - x3*t1 - x3*t2)/(x2*t1 + t2*x2))**2 + 1)
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


def xcross(p, **kwargs):
    """ draw a X cursor at p
    """
    mlab.points3d(p[0], p[1], p[2],
                  mode='2dcross',
                  scale_factor=.3,
                  color=(.5, .7, 0.),
                  **kwargs)
