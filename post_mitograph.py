# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 17:21:43 2016

@author: sweel_Rafelski
"""
import sys
import os
import os.path as op
import shutil as sh
import fnmatch as fn
import traceback
import errno
from collections import defaultdict
from wrappers import UsageError


def mkdir_exist(path):
    """
    Makes a folder if it does not already exists
    """
    try:
        os.makedirs(path)
        print "{} created!".format(path)
    except OSError as error:
        if error.errno != errno.EEXIST:
            raise
        print "{}/ already exists".format(path)


def rename_copy(folder, fname):
    """
    Returns the media label and filename with date appended
    """
    fname = fname.partition("_")
    date = folder.rsplit(os.sep)[-1]
    newname = "_".join((date, fname[-1]))
    return (fname[0], newname)


def postmito_process():
    dirs = None
    for root, dirs, files in os.walk(os.getcwd(), topdown=False):
        # search for foldernames starting with six numbers (dates)
        dirs[:] = [op.join(root, d) for d in dirs
                   if fn.fnmatch(d, "[0-9][0-9][0-9][0-9][0-9][0-9]")]

    D = defaultdict(dict)
    D['SkelVTK'] = []
    D['resampleFiles'] = []
    D['surfaceFiles'] = []
    for dircell in dirs:
        for root, dirs, files in os.walk(dircell):
            for f in files:
                if 'skeleton' in f:
                    D['SkelVTK'].append((root, f))
                elif 'resampled' in f:
                    D['resampleFiles'].append((root, f))
                elif 'surface' in f:
                    D['surfaceFiles'].append((root, f))

    basedir = op.dirname(D['SkelVTK'][0][0])
    for key in D:  # key == skelVTK, resampledFiles...
        mkdir_exist(op.join(basedir, key))

        for items in D[key]:  # item == (path_to_file, filename)
            media, filename = rename_copy(*items)
            subfolder = op.join(basedir, key, media)
            mkdir_exist(subfolder)
            print ("copying {}-->\n{}"
                   .format(op.join(*items), op.join(subfolder, filename)))

            # for compatibility with older naming system (refactor in future)
            if 'surface' in items[1]:
                filename = "_".join((media, filename))
            sh.copy(op.join(*items), op.join(subfolder, filename))
    return D


def main():
    """
    Post MitoGraph folder processing
    """
    try:
        try:
            os.chdir(op.expanduser(os.sep.join(
                ('~', 'Documents', 'Github', 'sweepython',
                 'WorkingData', 'mutants'))))
        except WindowsError:
            traceback.print_stack(limit=1)
            raise UsageError("Couldn't find folder, check path")

        postmito_process()
        return 0
#
    except UsageError as e:
        print e
        return 1

if __name__ == '__main__':
    sys.exit(main())
