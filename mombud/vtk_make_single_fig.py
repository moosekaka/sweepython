# -*- coding: utf-8 -*-
"""
Plot mitoskel network in with various scalar values
"""
import sys
import os
import os.path as op
import matplotlib.pyplot as plt
from mayavi import mlab
from pipeline.make_networkx import makegraph as mg
from mombud.functions import vtkvizfuncs as vz
from wrappers import swalk, UsageError
# pylint: disable=C0103


def main(scalar_select=None, save=False, **kwargs):
    """
    Plots various scalar types for a single mito network

    Args:
    -----

    scalar_select : list
        types of scalars to be plotted, choose from this list:\n
        ['bgfp', 'brfp', 'dy', 'dyr', 'gfp', 'rfp', 'weq']

    kwargs
    ------

    Defaults for optional keywords :

        `size` =(1200,800), `bgcolor` =(0.05, 0.05, 0.05),
        `rad`= 0.08
    """
    try:
        scal_dict = dict(dyr='DY_raw', dy='DY_minmax',
                         rfp='rRFP', gfp='rGFP',
                         brfp='bkstRFP', bgfp='bkstGFP', weq='WidthEq')
        if scalar_select is None:
            scalar_select = scal_dict.keys()

        plt.close('all')
        mlab.close(all=True)
        datadir = op.join(os.getcwd(), 'mutants',
                          'transformedData', 'filtered')
        filekey = 'NUM1_032016_011_RFPstack_030'

        vtkF = swalk(datadir, '*.vtk', start=0, stop=-4)

        data = vz.callreader(vtkF[filekey])
        _, _, nxgrph = mg(data, filekey)
        figone = mlab.figure(figure=filekey,
                             size=kwargs.get('size', (1200, 800)),
                             bgcolor=kwargs.get('bgcolor', (.0, .0, .0)))

        # graphviz of mitograph
        f, ax = plt.subplots()
        vz.nicegrph(nxgrph, ax)
        if save:
            f.savefig(op.join(datadir, 'pipeline_grph.png'))

        for key in scalar_select:
            mlab.clf(figure=figone)
            vz.cellplot(figone,
                        vtkF[filekey],
                        scalartype=scal_dict[key],
                        scalar_visible=True,
                        rad=kwargs.get('rad', 0.08),
                        legend=True)
            vz.labelbpoints(nxgrph,
                            bsize=.12,
                            esize=0.06,
                            bcol='light magenta',
                            ecol='cyan')
            mlab.view(0, 0, distance='auto')
            if save:
                mlab.savefig(op.join(datadir,
                                     'pipeline_%s.png' % scal_dict[key]))
            print "Done!"
            return 0

    except UsageError as e:
        print e
        return 1

# filelist and graph list
if __name__ == '__main__':
    sys.exit(main())
#    main(scalar_select=['dy', 'bgfp', 'weq'], save=True)
