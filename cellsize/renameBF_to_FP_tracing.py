# -*- coding: utf-8 -*-
"""
Created on Sun Aug 02 00:14:16 2015

@author: sweel
"""

import pandas as pd
import glob
import os

A = pd.read_table('table.txt',header=None)
rename_list = sorted(pd.Series(A.ix[:,0]).unique())
oldname_list = glob.glob('[0c]*tif')
print len(oldname_list)
print len(rename_list)

for e, f in enumerate(oldname_list):
    os.rename(f, '%s.tif' % rename_list[e])
