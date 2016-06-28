# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 15:21:01 2016
Convenience functions for opening and iterating files
@author: sweel_rafelski
"""
import os
import os.path as op
import fnmatch
from collections import defaultdict


class UsageError(Exception):
    """
    Class for user-facing (non-programming) errors
    """
    pass


def swalk(ddir, txt, start=None, stop=None):
    """
    wrapper caller for single level file dict, returns a dict

    Parameters
    ----------
    ddir : str
        path of dir

    txt : str
        text to be matched with fnmatch, shell format

    start : int
      start of slice

    end : int
      end of slice
    """
    vtf = dict()
    for root, _, files in os.walk(ddir):
        for f in files:
            if fnmatch.fnmatch(f, txt):
                vtf[f[start:stop]] = op.join(root, f)
    if len(vtf):
        return vtf
    else:
        raise UsageError('Search for file with ext. {} in dir {} failed'
                         .format(txt, ddir))


def ddwalk(ddir, txt, start=None, stop=None):
    """
    wrapper caller for multi level file dicts, returns a defaultdict

    Parameters
    ----------
    ddir : str
        path of dir

    txt : str
        text to be matched with fnmatch, shell format

    start : int
      start of slice

    end : int
      end of slice

    Returns
    -------
    vtk : defaultdict(dict)
        default two level nested dict
    """

    vtf = defaultdict(dict)
    for root, _, files in os.walk(ddir):
        for f in files:
            if fnmatch.fnmatch(f, txt):
                media = root.rsplit(os.sep, 1)[1]
                vtf[media][f[start:stop]] = op.join(root, f)
    if len(vtf):
        return vtf
    else:
        raise UsageError('Search for file with ext. {} in dir {} failed'
                         .format(txt, ddir))
