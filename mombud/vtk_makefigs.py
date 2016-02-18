# -*- coding: utf-8 -*-
"""
       visualize the skel and surface vtk file
"""
import math
import os
import os.path as op
from collections import defaultdict
import fnmatch
import matplotlib.pyplot as plt
from mayavi import mlab
from _make_networkx import makegraph as mg
import vtk
from mombud.vtk_viz import vtkvizfuncs as vf

# pylint: disable=C0103
vtkF = defaultdict(dict)
vtkS = defaultdict(dict)
plt.close('all')
mlab.close(all=True)
datadir = op.join(os.getcwd(), 'data')


def callreader(filepath):
    """
    Convenience wrapper for vtk reader call
    """
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(filepath)
    reader.Update()
    polydata = reader.GetOutput()
    return polydata


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
if __name__ == '__main__':
    for root, dirs, files in os.walk(op.join(datadir, 'pipelineFigs')):
        for i in files:
            if fnmatch.fnmatch(i, '*skeleton.vtk'):
                media = root.rsplit('\\', 1)[1]
                vtkF[media][i[5:-13]] = os.path.join(root, i)

    for root, dirs, files in os.walk(op.join(datadir, 'surfaceFiles')):
        for i in files:
            if fnmatch.fnmatch(i, '*surface.vtk'):
                media = root.rsplit('\\', 1)[1]
                vtkS[media][i[:-12]] = os.path.join(root, i)

    filekeys = {item: vtkF[media][item] for media
                in sorted(vtkF.keys()) for item
                in sorted(vtkF[media].keys())}

#    for key in sorted(filekeys.keys())[400:403]:
#        data = callreader(vtkF[key[:3]][key])
#        node_data, edge_data, nxgrph = mg(data, files)
#        figone = mlab.figure(figure=key,
#                             size=(800, 600),
#                             bgcolor=(.15, .15, .15))
#        vtkobj, _ = vf.cellplot(figone, filekeys[key])
#        vf.rendsurf(vtkS[key[:3]][key[4:]])
#        vf.labelbpoints(nxgrph)
    filekey = 'YPE_042715_018_RFPstack_052'
    labs = filekey[:3]
    data = callreader(vtkF[labs][filekey])
    node_data, edge_data, nxgrph = mg(data, files)

    figone = mlab.figure(figure=filekey,
                         size=(800, 600),
                         bgcolor=(.25, .25, .25))

    vtkobj, _ = vf.cellplot(figone, filekeys[filekey])
    vf.rendsurf(vtkS[filekey[:3]][filekey[4:]])
    vf.labelbpoints(nxgrph)

    #mlab.savefig(op.join(datadir, '%s.png' % filekey))
    #labellines(vtkobj)

    #figtwo = mlab.figure(figure='tubefig',
    #                     size=(800, 600),
    #                     bgcolor=(0., 0., 0.))
    #edgeplot(figtwo, vtkobj, 40)
