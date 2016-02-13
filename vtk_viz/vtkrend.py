# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 16:05:11 2015

@author: sweel
"""

"""
       visualize the skel and surface vtk file
"""
import matplotlib.pyplot as plt
import math
import os
import numpy as np
from mayavi import mlab
import cPickle as pickle
import fnmatch
from collections import defaultdict
# pylint: disable=C0103
vtkF = defaultdict(dict)
G = {}
plt.close('all')
mlab.close(all=True)


def cellplot(fig, filelist, filename):
    """draw vtk one whole cell
    """
    src = mlab.pipeline.open(filelist[filename])
    tube = mlab.pipeline.tube(src, figure=fig)
    tube.filter.radius = .03
    surfTube = mlab.pipeline.surface(tube)
    surfTube.actor.mapper.scalar_visibility = True
    #    LUT of the tube controlled as module manager
    mod_mngr = tube.children[0]
    mmgr = mod_mngr.scalar_lut_manager
    mmgr.scalar_bar.title = 'Intensity'
    mmgr.data_name = 'Intensity'
    src.point_scalars_name = 'Intensity'
#    mmgr.show_legend = True
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
    return fig, src


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
    surface1.actor.property.opacity = .13
    surface1.actor.mapper.scalar_visibility = False
    surface1.actor.property.color = (.9, .75, .0)


def labelbpoints(graph):
    """  brnch points and end points
    """
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
    bpts.glyph.glyph_source.glyph_source.radius = 0.05
    epts = mlab.pipeline.glyph(points2)
    epts.glyph.glyph_source.glyph_source.radius = 0.05
    #   coloring and opacity for the nodes
    bpts.actor.property.color = (1, .2, 1.0)
    bpts.actor.mapper.scalar_visibility = 0
    bpts.actor.property.opacity = 0.8
    bpts.glyph.glyph_source.glyph_source.phi_resolution = 16
    bpts.glyph.glyph_source.glyph_source.theta_resolution = 16

    epts.actor.property.color = (.1, 1, 1.0)
    epts.actor.mapper.scalar_visibility = 0
    epts.actor.property.opacity = 0.1
    epts.glyph.glyph_source.glyph_source.phi_resolution = 16
    epts.glyph.glyph_source.glyph_source.theta_resolution = 16


def labellines(vtksrc):
    """ plot cell/line IDs"""
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


def edgeplot(fig, vtksrc, cellid):
    """draw one edge of the vtk cell
    """
    dataset = vtksrc.outputs[0]
    line = dataset.get_cell(cellid)
    vals = dataset.point_data
    pts = line.points.to_array()
    ptid = line.point_ids
    scalvals = [vals.scalars[pid] for pid in ptid]
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
    mmgr.lut_mode = 'RdBu'
    mmgr.number_of_labels = 5
    mmgr.scalar_bar.label_format = '%4.f'
    mmgr.label_text_property.font_size = 12
    mmgr.scalar_bar_representation.position = [.85, .25]
    mmgr.scalar_bar_representation.position2 = [.1, .4]
    tube.filter.number_of_sides = 32

# =============================================================================
# filelist and graph list
# =============================================================================
for root, dirs, files in os.walk(os.getcwd()):
    for i in files:
        if fnmatch.fnmatch(i, '*skeleton.vtk'):
            media = root.rsplit('\\', 1)[1]
            vtkF[media][i[5:-13]] = os.path.join(root, i)
#        if fnmatch.fnmatch(i, '*grph.pkl'):
#            G.setdefault(root.rsplit('\\', 1)[1], []).append(
#                os.path.join(root, i))

filekeys = {item: vtkF[media][item] for media
            in sorted(vtkF.keys()) for item
            in sorted(vtkF[media].keys())}

for i in G:
    with open(G[i][0], 'rb') as inpt:
        temp = pickle.load(inpt)[2]
        G[i].append(temp)

filekey = 'YPE_042715_018_RFPstack_052'
labs = filekey[:3]
#grphs = {i.graph['cell']: G[labs][1][h] for h, i in enumerate(G[labs][1])}
#grphsrc = grphs[filekey]

# =============================================================================
# Draw cell using cellplot and edgeplot
# =============================================================================
figone = mlab.figure(figure=filekey,
                     size=(600, 400),
                     bgcolor=(0., 0., 0.))
fighand, vtkobj = cellplot(figone, filekeys, filekey)
rendsurf(vtkobj)
mlab.savefig('%s.png' % filekey)
#labelbpoints(grphsrc)
#labellines(vtkobj)

#figtwo = mlab.figure(figure='tubefig',
#                     size=(800, 600),
#                     bgcolor=(0., 0., 0.))
#edgeplot(figtwo, vtkobj, 40)
