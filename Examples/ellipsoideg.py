# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 16:31:56 2015

@author: sweel
"""

from mayavi.api import Engine
from mayavi.sources.api import ParametricSurface
from mayavi.modules.api import Surface
from mayavi import mlab
from tvtk.tools import visual

import numpy as np


def Arrow_From_A_to_B(x1, y1, z1, x2, y2, z2, scale=None):
    ar1 = visual.arrow(x=x1, y=y1, z=z1)
    ar1.length_cone = 0.4

    arrow_length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    if scale is None:
        ar1.actor.scale = [arrow_length, arrow_length, arrow_length]
    else:
        ar1.actor.scale = scale
    ar1.pos = ar1.pos / arrow_length
    ar1.axis = [x2 - x1, y2 - y1, z2 - z1]
    return ar1


engine = Engine()
engine.start()
scene = engine.new_scene()
scene.scene.disable_render = True  # for speed

visual.set_viewer(scene)

surfaces = []
for i in range(2):
    source = ParametricSurface()
    source.function = 'ellipsoid'
    engine.add_source(source)

    surface = Surface()
    source.add_module(surface)

    actor = surface.actor  # mayavi actor, actor.actor is tvtk actor
    # actor.property.ambient = 1 # defaults to 0 for some reason, ah don't
    # need it, turn off scalar visibility instead
    actor.property.opacity = 0.7
    actor.property.color = (0, 0, 1)  # tuple(np.random.rand(3))
    # don't colour ellipses by their scalar indices into colour map
    actor.mapper.scalar_visibility = False
    # gets rid of weird rendering artifact when opacity is < 1
    actor.property.backface_culling = True
    actor.property.specular = 0.1
    #actor.property.frontface_culling = True
    actor.actor.orientation = np.array([1, 0, 0]) * 360  # in degrees
    actor.actor.origin = np.array([0, 0, 0])
    actor.actor.position = np.array([0, 0, 0])
    actor.actor.scale = np.array([0.26490647, 0.26490647, 0.92717265])
    actor.enable_texture = True
    actor.property.representation = ['wireframe', 'surface'][i]
    surfaces.append(surface)

Arrow_From_A_to_B(0, 0, 0, 0.26490647, 0, 0, np.array([0.26490647, 0.4, 0.4]))
Arrow_From_A_to_B(0, 0, 0, 0, 0.26490647, 0, np.array([0.4, 0.26490647, 0.4]))
Arrow_From_A_to_B(0, 0, 0, 0, 0, 0.92717265, np.array([0.4, 0.4, 0.92717265]))

source.scene.background = (1.0, 1.0, 1.0)
scene.scene.disable_render = False  # now turn it on

# set the scalars, this has to be done some indeterminate amount of time
# after each surface is created, otherwise the scalars get overwritten
# later by their default of 1.0
for i, surface in enumerate(surfaces):
    vtk_srcs = mlab.pipeline.get_vtk_src(surface)
    print('len(vtk_srcs) = %d' % len(vtk_srcs))
    vtk_src = vtk_srcs[0]
    try:
        npoints = len(vtk_src.point_data.scalars)
    except TypeError:
        print('hit the TypeError on surface i=%d' % i)
        npoints = 2500
    vtk_src.point_data.scalars = np.tile(i, npoints)

# on pick, find the ellipsoid with origin closest to the picked coord,
# then check if that coord falls within that nearest ellipsoid, and if
# so, print out the ellispoid id, or pop it up in a tooltip

mlab.show()
