"""
Pipeline to normalize'raw' vtk files and make mito network graph
"""
import sys
import os
import os.path as op
import cPickle as pickle
import wrappers as wr
from pipeline import pipefuncs as pf

# pylint: disable=C0103
datadir = op.join(os.getcwd())

# filelist and graph list
if __name__ == '__main__':
    try:
        with open(op.join(datadir, 'fileMetas.pkl'), 'rb') as inpt:
            filemetas = pickle.load(inpt)
    except IOError:
        print "Error: Make sure you have file metadatas in working directory"

    try:
        vtkSkel = wr.ddwalk(op.join(datadir, 'SkelVTK'),
                            '*skeleton.vtk', stop=-13)
        vtkVolRfp = wr.ddwalk(op.join(datadir, 'resampledFiles'),
                              '*RF*resampled.vtk', stop=-14)
        vtkVolGfp = wr.ddwalk(op.join(datadir, 'resampledFiles'),
                              '*GF*resampled.vtk', stop=-14)
    except Exception:
        print "Error: check your filepaths"
        sys.exit()

    filekeys = {item: vtkSkel[lab][item] for lab
                in sorted(vtkSkel.keys()) for item
                in sorted(vtkSkel[lab].keys())}

    for lab in sorted(vtkSkel.keys())[-1:]:
        for key in sorted(vtkSkel[lab].keys())[-1:]:
            print 'processing %s: %s ' % (lab, key)
            data = pf.add_scalars(vtkSkel[lab][key],
                                  vtkVolRfp[lab][key],
                                  vtkVolGfp[lab][key.replace('RFP', 'GFP')])
            filename = op.join(datadir, 'Norm_%s_skeleton.vtk' % key)
            nm, rw, rb, gb, wq = pf.normSkel(data,
                                             filemetas['_'.join((lab, key))])
            calc = {'DY_minmax': nm,
                    'DY_raw': rw,
                    'bkstRFP': rb,
                    'bkstGFP': gb,
                    'WidthEq': wq}
            pf.writevtk(data, filename, **calc)
