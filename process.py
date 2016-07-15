# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 13:15:19 2016

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
# pylint: disable=C0103


def mkdir_exist(path):
    try:
        os.makedirs(path)
    except OSError as error:
        if error.errno != errno.EEXIST:
            raise


def movefilesup(path):
    """
    Move files from `Pos0` folder up to top level of `path` folder
    """
    count = 0
    for root, dirs, _ in os.walk(path, topdown=True):
        hasPos0 = [(root, d) for d in dirs if fn.fnmatch(d, "Pos0")]
        if hasPos0:  # list is not empty
            for r, _, files in os.walk(op.join(*hasPos0[0])):
                for f in files:
                    src = op.join(r, f)
                    dst = op.join(hasPos0[0][0], f)
                    count += 1
                    print "Moving {} -->\n{}".format(src, dst)
                    sh.move(src, dst)
    return count


def split_pad(folder):
    """
    repad the FOLDER names middle index to 00[0-9] format
    """
    for dirs in os.listdir(folder):
        olddir = op.join(folder, dirs)
        oldstr = op.basename(olddir)
        oldstrL = oldstr.split("_")
        newstr = "_".join((oldstrL[0], oldstrL[1].zfill(3)))
        newdir = olddir.replace(oldstr, newstr)
        print "renaming {} -->\n{}".format(olddir, newdir)
        try:
            os.rename(olddir, newdir)
        except WindowsError:
            pass


def repad(pth):
    """
    Helper function for split_pad
    """
    try:
        split_pad(pth)
    except IndexError:
        try:  # folder doesnt have "_", try subfolder
            for dirs in os.listdir(os.getcwd()):
                print "Now in {}".format(dirs)
                split_pad(op.join(pth, dirs))
        except IndexError:
            traceback.print_stack(limit=4)
            raise UsageError("Check folder paths")


def switch_labels(pth):
    for root, dirs, files in os.walk(pth):
        for f in files:
            if 'GFP' in f:
                try:
                    old_g = op.join(root, f)
                    old_r = old_g.replace('GFP', 'RFP')
                    if op.isfile(old_r):
                        print "Switchin labels for {}\n".format(old_g)
                        # save the gfp label to some long number
                        os.rename(old_g, op.join(root, "_27092450646347351L"))
                        # rename the rfp label to GFP
                        os.rename(old_r, old_r.replace('RFP', 'GFP'))
                        # rename the long number label to RFP
                        os.rename(op.join(root, "_27092450646347351L"), old_r)

                    else:
                        traceback.print_stack(limit=2)
                        raise FalseException

                except FalseException:
                    raise UsageError("\nRelabeled file does not exist")


def main():
    try:
        try:
            # change this path to where the preprocessed
            # raw tif image stacks are
            os.chdir(op.expanduser(os.sep.join(
                ('~', 'Documents', 'Github', 'sweepython',
                 'WorkingData', '071016', 'preprocessed'))))
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
