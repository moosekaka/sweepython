# -*- coding: utf-8 -*-
"""
Batch visualize skel and surface vtk files
"""
import os
import os.path as op
import matplotlib.pyplot as plt
from mayavi import mlab
from pipeline.make_networkx import makegraph as mg
from mombud.functions import vtkvizfuncs as vf
from wrappers import ddwalk, UsageError
# pylint: disable=C0103


def main(*args, **kwargs):
    """
    make a list of figures for mito networks

    Args
    ----
    start, stop, step : int
        arguments for slice(), i.e. [start:stop:step]


    kwargs
    ------
    save : Bool
        save toggle
    offscreen : Bool
        disables on screen rendering, won't display plot
    size, bgcolor : Tuple
        size and background color of Figure
    bsize, esize : Float
        branch and end points size
    """

    plt.close('all')
    mlab.close(all=True)
    inputdir = op.join(os.getcwd(), 'mutants')
    datadir = op.join(os.getcwd(), 'mutants', 'normalizedVTK')

    try:
        vtkF = ddwalk(datadir, '*skeleton.vtk', start=5, stop=-13)
        vtkS = ddwalk(op.join(inputdir, 'surfaceFiles'),
                      '*surface.vtk', stop=-12)
    except UsageError:
        raise

    filekeys = {item: vtkF[media][item] for media
                in sorted(vtkF.keys()) for item
                in sorted(vtkF[media].keys())}

    figone = mlab.figure(size=kwargs.get('size', (1200, 800)),
                         bgcolor=kwargs.get('bgcolor', (.05, .05, .05)))
    if kwargs.get('offscreen', False):
        figone.scene.off_screen_rendering = True

#    slicerange = kwargs.get('slicerange', (None, None, 100))
    for key in sorted(filekeys.keys())[slice(*args)]:
        mlab.clf(figure=figone)
        temp = key.partition("_")
        etype = temp[0]
        data = vf.callreader(vtkF[etype][key])
        node_data, edge_data, nxgrph = mg(data, key)
        vtkobj, _ = vf.cellplot(figone, filekeys[key])
        vf.rendsurf(vtkS[etype][key])
        vf.labelbpoints(nxgrph,
                        bsize=kwargs.get('bsize', 0.08),
                        esize=kwargs.get('esize', 0.08))

        if kwargs.get('save', False):
            mlab.savefig(op.join(datadir, '%s.png' % key))

# filelist and graph list
if __name__ == '__main__':
    main(1, 5, save=True, offscreen=True)
