# -*- coding: utf-8 -*-
"""
Search for dependents of function
@author: sweel_rafelski
"""

import os
import os.path as op
import fnmatch
import re
import sys


def searchfunc(fun, rootdir):
    """
    Return modules where fun is used
    """
    test = sys.argv[0].rsplit('/')[-1]
    print ('\nFunction %s used in module(s):\n' % fun)
    for root, dirs, files in os.walk(rootdir):
        for f in files:
            if f != test and fnmatch.fnmatch(f, '*.py'):  # exclude this module
                with open(op.join(root, f)) as script:
                    if script.read().find(fun) != -1:
                        print ('\t %s' % (op.join(root, f)))

if __name__ == '__main__':

    datadir = os.getenv('PYTHONPATH').split(';')
    scrloc = datadir[[num for num, val in enumerate(datadir)
                      if re.search('sweepython', val)][0]]


    searchfunc('vtk_viz', scrloc)

