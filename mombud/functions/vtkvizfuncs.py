# -*- coding: utf-8 -*-
"""
functions to render vtk output from pipeline
"""
import math
from collections import defaultdict
import vtk
import numpy as np
from mayavi import mlab
from mayavi.sources.api import ParametricSurface
from mayavi.sources.vtk_data_source import VTKDataSource
from tvtk.api import tvtk
import networkx as nx
from seaborn import xkcd_palette as scolor
from wrappers import UsageError
# pylint: disable=C0103


def rgbcol(colorname):
    """
    returns RGB float values based on xkcd color name input string
    """
    return scolor([colorname])[0]


def generate_color_labels(colors=None, labels=None):
    """
    Generates rgb values when passed a dictionary of xkcd color list
    and labels list."""

    if 'colors' is None:
        raise UsageError('must specify color list')
    if 'labels' is None:
        raise UsageError('must specify labels list')

    dic_label_colors = dict(zip(labels, colors))
    rgb_vals = dict(zip(colors, scolor(colors)))

    return dic_label_colors, rgb_vals


def nicegrph(graph, axinput, **kwargs):
    """
    Draw graph using graphviz

    **kwargs
    --------
    bcol, ecol : str
        branchpoint and endpoint colors

    grphtype : str ("neato")
        layout for graphviz, options:\n
        `neato` | `dot` | `twopi` | `circo` | `fdp`

    """
    cols = []
    sizes = []
    b_color = kwargs.pop('bcol', rgbcol('reddish pink'))
    e_color = kwargs.pop('ecol', rgbcol('blue'))
    grphtype = kwargs.pop('grphtype', 'neato')

    for _, attr in graph.nodes(data=True):
        if attr['degree'] == 1:
            cols.append(e_color)
            sizes.append(30)
        else:
            cols.append(b_color)
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
    fname : Str
        VTK file path

    Returns:
    --------
    src : VTKDataSource
        mayavi pipeline ready source object

    """
    dat = tvtk.PolyDataReader(file_name=fpath)
    dat.update()   # very IMPORTANT!!!
    src = VTKDataSource(data=dat.output)
    if hasattr(src, 'point_scalars_name'):
        src.point_scalars_name = 'DY_raw'  # make this default scalar
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
    return scalvals, mmgr


def cellplot(fig, filename, **kwargs):
    """
    Draw vtk one whole cell

    Parameters
    ----------
    fig : str
        Name of cell

    filename : str
        filename


    **kwargs
    --------
    rad : flt
        sets tube radius

    scalartype : str (`DY_raw`)
        point data type to plot on skeleton, can be of type:\n
        `DY_minmax` | `WidthEq` | `DY_raw` | `rRFP` | `rGFP` | `bkstRFP` |\
        `bkstGFP`

    scalar_visible : Bool (True)
        turns on/off heatmap scalars of vkt skel

    legend : Bool (True)
        turns on/off legend for heatmap

    Returns
    -------
    src: vtk object
        mayavi handle to the vtk obj
    surfTube: mayavi surface object
        mayavi pipeline surface
    """
    # config switch
    scalartype = kwargs.pop('scalartype', 'DY_raw')
    scalar_visible = kwargs.pop('scalar_visible', True)
    src = mlab.pipeline.open(filename)
    src.point_scalars_name = scalartype
    # tube object
    tube = mlab.pipeline.tube(src, figure=fig)
    tube.filter.radius = kwargs.pop('rad', .05)
    tube.filter.number_of_sides = 32
    # surface object
    surfTube = mlab.pipeline.surface(tube)
    surfTube.actor.mapper.scalar_visibility = scalar_visible
    # LUT config
    mmgr = surfTube.module_manager.scalar_lut_manager
    mmgr.scalar_bar.title = scalartype
    mmgr.data_name = scalartype
    mmgr.show_legend = kwargs.pop('legend', True)
    mmgr.reverse_lut = True
    mmgr.lut_mode = 'RdBu'
    mmgr.number_of_labels = 4
    mmgr.label_text_property.font_size = 10
    mmgr.title_text_property.font_size = 10
    mmgr.scalar_bar.label_format = '%4.1f'
    mmgr.scalar_bar_representation.position = [.85, .2]
    mmgr.scalar_bar_representation.position2 = [.1, .4]
    return src, surfTube


def rendsurf(vtksrc, **kwargs):
    """   add surface.vtk file to pipeline

    Parameters
    ----------
    color : str
        xkcd color of surface
    alpha : flt
        opacity

    Returns : surface1
        mlab surface actor obj
    """

    color = rgbcol(kwargs.pop('color', 'yellow'))
    alpha = kwargs.pop('alpha', .15)

    src2 = mlab.pipeline.open(vtksrc)
    surface1 = mlab.pipeline.surface(src2)
    surface1.actor.property.opacity = alpha
    surface1.actor.mapper.scalar_visibility = False
    surface1.actor.property.color = color
    return surface1


def labelbpoints(graph, **kwargs):
    """  brnch points and end points

    Parameters
    ----------
    graph : network graph
        graph object

    kwargs:
    ------
    bsize : flt
        size of branchpoint
    esize : flt
        size of endpoint
    bcol : str
        xkcd color of branchpoint
    ecol : str
        xkcd color of endpoint
    """

    Nodes = [nattr['coord'] for _, nattr
             in graph.nodes(data=True)
             if nattr['degree'] > 1]
    Ends = [nattr['coord'] for _, nattr
            in graph.nodes(data=True)
            if nattr['degree'] == 1]

    D = defaultdict(lambda: defaultdict(dict))
#    colors = ["burnt orange", "pale blue"]
#    palette = dict(zip(colors, scolor(colors)))

    D['bpts']['obj'] = Nodes
    D['bpts']['size'] = kwargs.pop('bsize', .2)
    D['bpts']['col'] = rgbcol(kwargs.pop('bcol', 'hot pink'))
    D['epts']['obj'] = Ends
    D['epts']['size'] = kwargs.pop('esize', .15)
    D['epts']['col'] = rgbcol(kwargs.pop('ecol', 'pale blue'))
    for key in D:
        xyz = np.array(D[key]['obj'])
        points = mlab.pipeline.scalar_scatter(  # branchpoints
            xyz[:, 0], xyz[:, 1], xyz[:, 2])
        D[key]['glyph'] = mlab.pipeline.glyph(points)
        D[key]['glyph'].glyph.glyph_source.glyph_source.radius = D[key]['size']
        D[key]['glyph'].actor.property.color = D[key]['col']
        D[key]['glyph'].glyph.glyph_source.glyph_source.phi_resolution = 16
        D[key]['glyph'].glyph.glyph_source.glyph_source.theta_resolution = 16
        D[key]['glyph'].actor.mapper.scalar_visibility = 0


def getelipspar(fname, dfin, **kwargs):
    """
    Return parameters for ellipse from cell tracing dataframe

    Parameters:
    -----------
    fname : Str
        cell name id, e.g. `MFB1_032016_002_RFPstack_000`
    dfin : DataFrame
        Dataframe file from pandas.csv read of a cell tracing data
    useold : Str
        switch to specify the coordinate transform between VTK and
        imageJ for BF cell tracings, default=False

    Returns:
    --------
    dfout : DataFrame
        dataframe of ellipse parameters (major, minor radius, center etc.)
    """
    # selection returns a view by default, we want a copy!
    dfout = dfin[dfin.cell == fname]
    dfout = dfout.reset_index(drop=True)
    dfout['Major'] = dfout.Major * .055/2
    dfout['Minor'] = dfout.Minor * .055/2
    dfout['vol'] =  \
        4 / 3 * np.pi * dfout.Major * dfout.Minor**2
    dfout = dfout.sort_values(by='vol')
    dfout.index = ['bud', 'mom']
    useold = kwargs.pop('useold', False)
    # reverse the y-coordinate system (VTK vs ImageJ), note for new data,
    # assumes 250 x 250 pixel BF
    if not useold:
        dfout['center'] = zip((dfout.X) * .055, (250 - dfout.Y) * .055)
    # old dataset, used different centering scheme for BF
    else:
        dfout['center'] = zip((dfout.X - 25) * .055, (225 - dfout.Y) * .055)
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


def inverse_tr(rotmatrix, dist):
    """
    returns the inverse of arrowvect transform, i.e. transform the original
    orientation of arrowvect back to a x-axis unit vector
    """
    tr_filt = tvtk.Transform()
    rotmatrix.transpose()
    tr_filt.translate(np.negative(dist))
    tr_filt.post_multiply()  # translate, THEN rotate
    tr_filt.concatenate(rotmatrix)
    tr_filt.translate([-1., 0, 0])  # account for lenght of unit vector
    return tr_filt
