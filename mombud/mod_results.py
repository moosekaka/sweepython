# -*- coding: utf-8 -*-
"""
Created on Fri May 13 12:08:48 2016
@author: sweel_rafelski
combine fit ellipse Results measurement into suitable form for mombud analysis
"""

import os
import os.path as op
import fnmatch as fn
import pandas as pd

for root, dirs, files in os.walk(os.getcwd(), topdown=False):
    dirs[:] = [op.join(root, d) for d in dirs if fn.fnmatch(d, "032[0-9]16*")]

P = pd.DataFrame()
for res in dirs:
    date = res.rsplit(os.sep)[-1]
    temp = pd.read_table(res + "/Stacks/Results.xls")
    temp['date'] = date
    temp['type'] = temp['Label'].apply(
        lambda x: x.rpartition(":")[-1].partition("_")[0])
    temp['first_str'] = temp['Label'].apply(
        lambda x: x.split(":")[0] + ':')
    temp['newname'] = temp['Label'].str.rpartition(
        ":").get(2).str.partition("_").get(2)
    temp['newname2'] = temp[['type',
                             'date',
                             'newname']].apply(lambda x: '_'.join(x), axis=1)
    temp['newname2'] = temp[['first_str',
                            'newname2']].apply(lambda x: ''.join(x), axis=1)
    P = P.append(temp, ignore_index=True)


P.index += 1  # change to one based index
P = P[['newname2', 'X', 'Y', 'Major', 'Minor', 'Angle']]
P = P.rename(columns={'newname2': 'Label'})
P.to_csv(op.join(root, "mutants", "Results.txt"), sep='\t')
