# -*- coding: utf-8 -*-
"""
Main module for mito network normalization
@author: sweel_rafelski
"""
import os
import os.path as op
import cPickle as pickle
import re
from collections import defaultdict
from post_mitograph import mkdir_exist
from pipeline import pipefuncs as pf
# pylint: disable=C0103

# data folder should contain subfolders of cell conditions, each condition
# having a skeleton, channel 1 and channel2 resampled VTK file for each cell
# quick check is the number of files in each subfolder is divisible by 3
datafolder = op.join('.', 'mutants', 'pre_normalized')
# output will be saved here
savefolder = op.join('.', 'mutants', 'normalized_vtk')
# modify the 'RFP' and 'GFP' keys to suit
SEARCHDICT = defaultdict(dict,
                         {'resampled': {'RFP': 'ch2',
                                        'GFP': 'ch1'},
                          'skeleton': {'RFP': 'skel'}})
# master dictionary of file paths stored here, with 'skel', 'ch1' and 'ch2' as
# the top level keys
vtks = defaultdict(dict)


def readfolder(folder):
    """
    Helper function to return a defaultdict of VTK file paths and labels
    """
    for subfolder in os.listdir(folder):
        if op.isdir(op.join(folder, subfolder)):
            for files in os.listdir(op.join(folder, subfolder)):
                vtk_type = re.search(r'(skeleton|resampled)', files)
                channel_type = re.search(r'([GR]FP)\w+\d+', files)
                if vtk_type and channel_type:

                    cell_id = '_'.join([subfolder,
                                        channel_type.
                                        string[:channel_type.end()]])
                    # this is just used to determine which type of VTK file
                    # (ie skeleton or voxel/channel file to use)
                    prefix = (SEARCHDICT.
                              get(vtk_type.group(1)).
                              get(channel_type.group(1)))

                    # if no prefix found, means skeleton.vtk is based on
                    # ch1 which is not what we want
                    if prefix:
                        vtks[prefix][cell_id] = op.join(folder,
                                                        subfolder, files)
    return vtks


def main():
    """
    Pipeline to normalize'raw' vtk files and make mito network graph
    """

    try:
        with open(op.join('.', 'mutants', 'filemetas_new.pkl'), 'rb') as inpt:
            filemetas = pickle.load(inpt)
    except WindowsError:
        print "Error: Make sure you have file metadatas in working directory"

    paths = readfolder(datafolder)

    cells = paths['skel']
    keys = sorted(cells.keys(), reverse=True)
    while keys:
        key = keys.pop()
        mkdir_exist(savefolder)
        savename = op.join(savefolder,
                           'Normalized_{}_mitoskel.vtk'.format(key))

        data, v1, v2 = pf.point_cloud_scalars(
            paths['skel'][key],
            paths['ch1'][key.replace('RFP', 'GFP')],
            paths['ch2'][key])
        dict_output = pf.normalize_skel(data, v1, v2,
                                        backgroundfile=filemetas[key[:-4]])
        pf.write_vtk(data, savename, **dict_output)
        print "{} normalized!".format(key)

if __name__ == '__main__':
    main()
