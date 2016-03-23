# -*- coding: utf-8 -*-
"""
# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
# Copyright (c) 2008, Enthought, Inc.
# License: BSD Style.
"""

from mayavi import mlab

# To access any VTK object, we use 'tvtk', which is a Python wrapping of
# VTK replacing C++ setters and getters by Python properties and
# converting numpy arrays to VTK arrays when setting data.
from tvtk.api import tvtk
import numpy as np
from transformations import transformations as tf

def colored_lines():
  # Create three points. Join (Origin and P0) with a red line and
  # (Origin and P1) with a green line
  origin = [0.0, 0.0, 0.0]
  p0 = [1.0, 0.0, 0.0]
  p1 = [0.0, 1.0, 0.0]

  # Create a vtkPoints object and store the points in it
  pts = tvtk.Points()
  pts.insert_next_point(origin)
  pts.insert_next_point(p0)
  pts.insert_next_point(p1)

  # Setup two colors - one for each line
  red = [255, 0, 0]
  green = [0, 255, 0]

  # Setup the colors array
  colors = tvtk.UnsignedCharArray()
  colors.number_of_components = 3
  colors.name = "Colors"

  # Add the colors we created to the colors array
  colors.insert_next_tuple_value(red)
  colors.insert_next_tuple_value(green)

  # Create the first line (between Origin and P0)
  line0 = tvtk.Line()
  line0.point_ids.set_id(0,0) # the second 0 is the index of the Origin in the vtkPoints
  line0.point_ids.set_id(1,1) # the second 1 is the index of P0 in the vtkPoints

  # Create the second line (between Origin and P1)
  line1 = tvtk.Line()
  line1.point_ids.set_id(0,0) # the second 0 is the index of the Origin in the vtkPoints
  line1.point_ids.set_id(1,2) # 2 is the index of P1 in the vtkPoints

  # Create a cell array to store the lines in and add the lines to it
  lines = tvtk.CellArray()
  lines.insert_next_cell(line0)
  lines.insert_next_cell(line1)

  # Create a polydata to store everything in
  linesPolyData = tvtk.PolyData()

  # Add the points to the dataset
  linesPolyData.points =pts

  # Add the lines to the dataset
  linesPolyData.lines = lines

  # Color the lines - associate the first component (red) of the
  # colors array with the first component of the cell array (line 0)
  # and the second component (green) of the colors array with the
  # second component of the cell array (line 1)
  linesPolyData.cell_data.scalars =colors

  return linesPolyData

"""
make little polydata indicating a pose, sphere at the pose origin,
cylinder in z-axis direction flag at the end of cylinder in x-direction
"""
def pose_indicator(scale=0.1):
  # setting up the sphere at pose origin
  sphere = tvtk.SphereSource()
  sphere.radius = 0.1
  sphere.theta_resolution = sphere.theta_resolution  * 2
  sphere.phi_resolution = sphere.phi_resolution  * 2
  sphere.update()
  sphere_poly = sphere.output

  colors = np.empty((sphere_poly.number_of_cells, 3), dtype=np.int)
  colors[:,0] = 0 # red
  colors[:,1] = 255 # green
  colors[:,2] = 0 # blue
  sphere_color = tvtk.UnsignedCharArray()
  sphere_color.from_array(colors)
  sphere_color.name = 'Colors'

  sphere_poly.cell_data.scalars = sphere_color

  # setting up cylinder in z direction
  line = tvtk.LineSource()
  line.point1 = [0.0, 0.0, 0.0]
  line.point2 = [0.0, 0.0, 1.0]

  tube_filter = tvtk.TubeFilter(input=line.output)
  tube_filter.capping = 1
  tube_filter.radius = 0.05
  tube_filter.number_of_sides = 8
  tube_filter.update()
  tube_poly = tube_filter.output

  colors = np.empty((tube_poly.number_of_cells, 3), dtype=np.int)
  colors[:,0] = 0 # red
  colors[:,1] = 0 # green
  colors[:,2] = 255 # blue
  tube_color = tvtk.UnsignedCharArray()
  tube_color.from_array(colors)
  tube_color.name = 'Colors'

  tube_poly.cell_data.scalars = tube_color

  # setting up plane in x-direction at top of marker
  size = 0.25
  plane = tvtk.PlaneSource()
  plane.origin = [0.0, 0.0, 1.0]
  plane.point1 = [0.0, 0.0, 1.0 - size]
  plane.point2 = [size, 0.0, 1.0]
  plane.update()
  plane_poly = plane.output

  colors = np.empty((plane_poly.number_of_cells, 3), dtype=np.int)
  colors[:,0] = 255 # red
  colors[:,1] = 0 # green
  colors[:,2] = 0 # blue
  plane_color = tvtk.UnsignedCharArray()
  plane_color.from_array(colors)
  plane_color.name = 'Colors'

  plane_poly.cell_data.scalars = plane_color

  # combine into one polydata object
  combined = tvtk.AppendPolyData()
  combined.add_input(sphere_poly)
  combined.add_input(tube_poly)
  combined.add_input(plane_poly)
  # combined.update()
  combined_poly = combined.output

  # scale combined output
  scale_mat = np.eye(4)
  scale_mat[0,0] = scale
  scale_mat[1,1] = scale
  scale_mat[2,2] = scale

  scale = tvtk.Transform()
  scale.set_matrix(scale_mat.flatten())

  scaler = tvtk.TransformPolyDataFilter(input=combined_poly)
  scaler.transform = scale

  scaled_poly = scaler.output
  return scaled_poly


"""
loads an obj mesh for ploting with plot_mesh
returned is a tvtk.PolyData
"""
def load_obj(filename):
  obj = tvtk.OBJReader()
  obj.file_name = filename
  mesh = obj.output
  return mesh

"""
displays polydata in a given pose, optional adds a coord frame to the origin with a given
scale (scale == 0.0 --> no coord frame displayed)
"""
def plot_polydata(polydata, T_obj2world=np.eye(4), axes_scale=0.0, opacity=0.2, color=(0,0,1), figure=None):
  mapper = tvtk.PolyDataMapper(input=polydata)

  p = tvtk.Property(opacity=opacity, color=color)

  actor = tvtk.Actor(mapper=mapper, property=p)

  # set pose
  mat = tvtk.Matrix4x4()
  mat.from_array(T_obj2world)
  actor.user_matrix = mat

  if not figure:
    figure = mlab.figure()

  figure.scene.add_actor(actor)

  if axes_scale > 0.0:
    plot_coordframe(axes_scale, figure)

  return figure

"""
plots a coord frame at the scenes origin with definable scale
"""
def plot_coordframe(scale=1.0, figure=None):
  axes = tvtk.AxesActor()
  axes.axis_labels = 0
  axes.total_length = [scale, scale, scale]

  if not figure:
    figure = mlab.figure()

  figure.scene.add_actor(axes)

  return figure

def main():
  model = load_obj('./mug.obj')
  p = pose_indicator(0.03)

  angles = np.random.rand(10) * 2 * np.pi
  axis = np.array([0,1,0])

  quaternions = np.empty((len(angles), 4), dtype=np.double)
  quaternions[:,0] = np.cos(angles / 2.0)
  quaternions[:,1:4] =  axis[None, :]
  quaternions[:,1:4] = quaternions[:,1:4] * np.sin(angles / 2.0)[:,None]

  f = mlab.figure()

  for q in quaternions:
    angles = tf.euler_from_quaternion(q)
    mat = tf.compose_matrix(angles=angles, translate=np.array([0,0,0]))
    plot_polydata(p, mat, opacity=1.0, figure=f)

  plot_polydata(model, figure=f)

  plot_coordframe(0.1, f)
  mlab.show(stop=True)

if __name__ == '__main__':
  main()