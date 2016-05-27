# -*- coding: utf-8 -*-
"""
functions to render vtk output from pipeline
"""
import math
import vtk
import numpy as np
from mayavi import mlab
from mayavi.sources.api import ParametricSurface
from mayavi.sources.vtk_data_source import VTKDataSource
from tvtk.api import tvtk
import config
import networkx as nx
# pylint: disable=C0103


def nicegrph(graph, axinput, grphtype='neato'):
    # ‘neato’|’dot’|’twopi’|’circo’|’fdp’|
    cols = []
    sizes = []

    for n, attr in graph.nodes(data=True):
        if attr['degree'] == 1:
            cols.append('#3366FF')
            sizes.append(30)
        else:
            cols.append('#FF5050')
            sizes.append(50)
    nx.draw_graphviz(graph,
                     prog=grphtype,
                     ax=axinput,
                     node_size=sizes,
                     node_color=cols)


def callreader(filepath):
    """
    Convenience wrapper for vtk reader call
    ***THIS RETURNS A VTK version polydata, for a mayavi source, use
    setup_vtk_source***
    """
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(filepath)
    reader.Update()
    polydata = reader.GetOutput()
    return polydata


def setup_vtk_source(fpath):
    """
    Given a VTK file name `fpath`, this creates a tvtk reader
    and returns a mayavi pipeline ready source object.
    created.

    Parameters:
    -----------
    fname: Str
        VTK file path

    Returns:
    --------
    src: VTKDataSource
        mayavi pipeline ready source object

    """
    dat = tvtk.PolyDataReader(file_name=fpath)
    dat.update()   # very IMPORTANT!!!
    src = VTKDataSource(data=dat.output)
    return src


def labellines(vtksrc):
    """
    Plot cell/line IDs"""
    dataset = vtksrc.outputs[0]
    for line in range(dataset.number_of_lines):
        cellhand = dataset.get_cell(line)
        if cellhand.length2 < 66.5 and cellhand.length2 > .5:
            tL = int(math.ceil(cellhand.number_of_points / 2))
            x = cellhand.points[tL][0]
            y = cellhand.points[tL][1]
            z = cellhand.points[tL][2]
            print"Line %2d has length: %8.4f" % (line,
                                                 cellhand.number_of_points)
            mlab.text3d(x, y, z, '%s' % line, scale=0.15)


def edgeplot(fig, vtksrc, cellid, scalartype='DY_raw'):
    """
    Draw one edge of the vtk cell. Uses tvtk functions, not VTK.

    """
    dataset = tvtk.to_tvtk(vtksrc)
    line = dataset.get_cell(cellid)
    vals = dataset.point_data.get_array(scalartype)
    pts = line.points.to_array()
    ptid = line.point_ids
    scalvals = [vals[int(pid)] for pid in ptid]

    src = mlab.pipeline.line_source(pts[:, 0],
                                    pts[:, 1],
                                    pts[:, 2],
                                    scalvals)
    tube = mlab.pipeline.tube(src, figure=fig)
    tube.filter.radius = .03
    surfTube = mlab.pipeline.surface(tube)
    surfTube.actor.mapper.scalar_visibility = True
    mod_mngr = tube.children[0]
    mmgr = mod_mngr.scalar_lut_manager
    mmgr.show_legend = True
    mmgr.reverse_lut = True
    mmgr.use_default_range = False
    mmgr.lut.set(range=[dataset.scalar_range[0],
                        dataset.scalar_range[1]])
#    mmgr.lut.set(range=[min(scalvals),
#                        max(scalvals)])
    mmgr.lut_mode = 'RdBu'
    mmgr.number_of_labels = 5
    mmgr.scalar_bar.label_format = '%4.f'
    mmgr.label_text_property.font_size = 12
    mmgr.scalar_bar_representation.position = [.85, .25]
    mmgr.scalar_bar_representation.position2 = [.1, .4]
    tube.filter.number_of_sides = 32


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
        bkstGFP
    rad : flt
        tube radius

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
    """   add surface.vtk file to pipeline

    Parameters
    ----------
    color: flt
        color of surface
    alpha : flt
        opacity
    """

    color = kwargs.pop('color', (0.9, .8, .1))
    alpha = kwargs.pop('alpha', .15)

    src2 = mlab.pipeline.open(vtksrc)
    surface1 = mlab.pipeline.surface(src2)
    surface1.actor.property.opacity = alpha
    surface1.actor.mapper.scalar_visibility = False
    surface1.actor.property.color = color


def labelbpoints(graph, **kwargs):
    """  brnch points and end points

    Parameters
    ----------
    bsize : flt
        size of branchpoint
    esize : flt
        size of endpoint
    bcol: tuple
        color of branchpoint
    ecol : tuple
        color of endpoint

    """
    bcol = kwargs.pop('bcol', (1, .2, 1.0))
    ecol = kwargs.pop('ecol', (.1, 1, 1.0))
    bsize = kwargs.pop('bsize', .15)
    esize = kwargs.pop('esize', .1)

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
    bpts.glyph.glyph_source.glyph_source.radius = bsize
    epts = mlab.pipeline.glyph(points2)
    epts.glyph.glyph_source.glyph_source.radius = esize

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


#def getelipspar(filename, df):
#    """ parameters for ellipse from cell tracing """
#    dftemp = df[df.cell == filename]
#    dftemp.reset_index(1, drop=True, inplace=True)
#    return dftemp


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


#def arrowvect(B, A, C):
#    """draws a vector based on base, B and tip, A.
#    calculates the transformation matrix trans and returns it along with the
#    rotation matrix
#    """
#    normalizedX = np.zeros(3)
#    normalizedY = np.zeros(3)
#    normalizedZ = np.zeros(3)
#    AP = np.zeros(3)
#    math = tvtk.Math()
#    math.subtract(A, B, normalizedX)  # normalizedX is the arrow unit vector
#    math.subtract(C, B, AP)
#    length = math.norm(normalizedX)
#    math.normalize(normalizedX)
#    math.normalize(AP)  # another unit vector used to fix the local x-y plane
#
#    x1, x2, x3 = normalizedX
#    t1, t2, t3 = AP
#    l3 = -t3/(t1+t2)
#    m3 = (t3*x1 - x3*t1 - x3*t2) / (x2*t1 + t2*x2)
#    D = np.sqrt((t3 / (t1 + t2))**2 +
#                ((t3*x1 - x3*t1 - x3*t2) / (x2*t1 + t2*x2))**2 + 1)
#    z1 = l3/D
#    z2 = m3/D
#    z3 = 1/D
#    normalizedZ = np.array([z1, z2, z3])
#    math.cross(normalizedZ, normalizedX, normalizedY)
#    matrix = tvtk.Matrix4x4()
#    matrix.identity()
#    for el in range(3):  # rotation matrix to x-axis
#        matrix.set_element(el, 0, normalizedX[el])
#        matrix.set_element(el, 1, normalizedY[el])
#        matrix.set_element(el, 2, normalizedZ[el])
#    trans = tvtk.Transform()
#    trans.translate(B)  # translate origin to base of arrow
#    trans.concatenate(matrix)  # rotate around the base of arrow
#    trans.scale(length, length, length)
#    return trans, matrix, length


def xcross(p, **kwargs):
    """ draw a X cursor at p
    """
    mlab.points3d(p[0], p[1], p[2],
                  mode='2dcross',
                  scale_factor=.3,
                  color=(.5, .7, 0.),
                  **kwargs)


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


def dircos_kernel(X, Y):
    """
    Calculates direction cosines compononent Z from vectors X, Y
    """
    x1, x2, x3 = X
    t1, t2, t3 = Y
    l3 = -t3 / (t1 + t2)
    m3 = (t3 * x1 - x3 * t1 - x3 * t2) / (x2 * t1 + t2 * x2)
    D = np.sqrt((t3 / (t1 + t2))**2 +
                ((t3 * x1 - x3 * t1 - x3 * t2) / (x2 * t1 + t2 * x2))**2 + 1)
    z1 = l3 / D
    z2 = m3 / D
    z3 = 1 / D
    return np.array([z1, z2, z3])


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
    normalizedX = np.subtract(tip, base)   # normalizedX = arrow unit vector
    AP = np.subtract(neck, base)
    length = np.linalg.norm(normalizedX)
    normalizedX = normalizedX / length
    # another unit vector used to fix the local x-y plane
    AP = AP / np.linalg.norm(AP)
    normalizedZ = dircos_kernel(normalizedX, AP)
    normalizedY = np.cross(normalizedZ, normalizedX)
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
