# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 17:54:10 2015
       functions to render vtk output from pipeline
"""
import os
import numpy as np
from mayavi import mlab
from mayavi.sources.api import ParametricSurface
from tvtk.api import tvtk
import random
from collections import defaultdict
import config
# pylint: disable=C0103
vtkF = defaultdict(dict)


def countFuncCall():
    """counter function
    """
    config.counter += 1
    return config.counter


def cellplot(fig, filelist, filename):
    """draw vtk one whole cell
    """
    src = mlab.pipeline.open(filelist[filename])
    tube = mlab.pipeline.tube(src, figure=fig)
    tube.filter.radius = .07
    src.point_scalars_name = 'DY_raw'
    surfTube = mlab.pipeline.surface(tube)
    surfTube.actor.mapper.scalar_visibility = True
    return fig, src, surfTube


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


def rendsurf(vtksrc):
    """   add surface.vtk file to pipeline"""
    surfpath = vtksrc.base_file_name.partition('NormFiles')
    path2 = surfpath[2].partition('Norm_')
    surfpath2 = os.path.join(surfpath[0],
                             'surfaceFiles',
                             path2[2][:3],
                             '%s_surface.vtk' % path2[2][4:-13])
    src2 = mlab.pipeline.open(surfpath2)
    surface1 = mlab.pipeline.surface(src2)
    surface1.actor.property.opacity = .8
    surface1.actor.mapper.scalar_visibility = False
    surface1.actor.property.color = (.4, .4, .9)


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
