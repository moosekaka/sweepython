import os
import os.path as op
import cPickle as pickle
from wrappers import ddwalk, UsageError
from post_mitograph import mkdir_exist
from pipeline import pipefuncs as pf
from pipeline import make_networkx as mn
# pylint: disable=C0103
datadir = op.join(os.getcwd())


def main():
    """
    Pipeline to normalize'raw' vtk files and make mito network graph
    """
    try:
        with open(op.join(datadir, 'mutants', 'fileMetas.pkl'), 'rb') as inpt:
            filemetas = pickle.load(inpt)
    except IOError:
        print "Error: Make sure you have file metadatas in working directory"

    try:
        vtkSkel = ddwalk(op.join(datadir, 'mutants', 'SkelVTK'),
                         '*RFP*skeleton.vtk', stop=-13)
        vtkVolRfp = ddwalk(op.join(datadir, 'mutants', 'resampledFiles'),
                           '*RF*resampled.vtk', stop=-14)
        vtkVolGfp = ddwalk(op.join(datadir, 'mutants', 'resampledFiles'),
                           '*GF*resampled.vtk', stop=-14)
    except UsageError:
        raise

    for lab in sorted(vtkSkel.keys())[:]:
        folder = op.join(datadir, 'mutants', 'normalizedVTK', lab)
        mkdir_exist(folder)

        for key in sorted(vtkSkel[lab].keys())[:]:
            print 'processing %s: %s ' % (lab, key)
            data = pf.add_scalars(vtkSkel[lab][key],
                                  vtkVolRfp[lab][key],
                                  vtkVolGfp[lab][key.replace('RFP', 'GFP')])
            filename = op.join(folder, 'Norm_%s_%s_skeleton.vtk' % (lab, key))
            (nm, rw, rb, gb,
             wq) = pf.normSkel(data, filemetas['_'.join((lab, key[:-4]))])

            calc = {'DY_minmax': nm,
                    'DY_raw': rw,
                    'bkstRFP': rb,
                    'bkstGFP': gb,
                    'WidthEq': wq}
            pf.writevtk(data, filename, **calc)
#            nl, el, nxgrph = mn.makegraph(data, "_".join((lab, key)))

if __name__ == '__main__':
    main()
