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
from wrappers import FalseException, UsageError


def mkdir_exist(path):
    try:
        os.makedirs(path)
    except OSError as error:
        if error.errno != errno.EEXIST:
            raise


def renamecopy(fname, folder, word):
    """
    macro to change filename label to date and return path of file
    """
    if fn.fnmatch(filename, word):
        fname = fname.partition("_")
        pth = op.join(op.dirname(folder), fname[0])
        date = folder.rsplit(os.sep)[-1]
        newname = "_".join((date, fname[-1]))
        return (pth, newname)


def postmito_process():
    for root, dirs, files in os.walk(os.getcwd(), topdown=False):
        # search for foldernames starting with six numbers (dates)
        dirs[:] = [op.join(root, d) for d in dirs
                   if fn.fnmatch(d, "[0-9][0-9][0-9][0-9][0-9][0-9]")]

    for dircell in dirs:
        for root, dirs, files in os.walk(dircell):
            for f in files:
                P = renamecopy(f, root, '*resampled*')  # change
                if P:
                    mkdir_exist(P[0])
                    sh.copy(op.join(root, f), op.join(P[0], P[1]))


    for root, dirs, files in os.walk(os.getcwd()):
        for f in files:
            if fn.fnmatchcase(f, '*surface*'):
                folder = root.rsplit(os.sep)[-1]
                old = op.join(root, f)
    #            new = old + ".vtk"
                new = op.join(root, '%s_%s' % (folder, f))
                os.rename(old, new)


def main():
    try:
        try:
            os.chdir(op.expanduser(os.sep.join(
                ('~', 'Documents', 'Github', 'sweepython',
                 'WorkingData', 'mutants'))))
        except WindowsError:
            traceback.print_stack(limit=1)
            raise UsageError("Couldn't find folder, check path")

        path = os.getcwd()
        print "Changed to {}".format(path)

        cmove = movefilesup(path)
        print "Moved {} files!".format(cmove)

        repad(path)

        while True:
            try:
                invar = raw_input("Do you wish to switch the labels for GFP/RFP??\ Press Y for yes")
                if invar == 'Y':
                    switch_labels(path)
                elif invar == 'Q' or invar == 'N':
                    print 'Quitting!'
                    break

            except KeyError as e:
                print "{} is not a valid input, try again (Y| [Q | N])".format(e)
                continue


        return 0

    except UsageError as e:
        print e
        return 1

if __name__ == '__main__':
    sys.exit(main())