# -*- coding: utf-8 -*-
"""
Created on Fri May 13 12:08:48 2016
@author: sweel_rafelski
get background values for filemetas.pkl
"""

import os
import cPickle as pickle
import fnmatch as fn
import pandas as pd

for r, _, files in os.walk(os.getcwd(), topdown=False):
    F = [fil for fil in files if fn.fnmatch(fil, "*[0-9][0-9]*xls")]

P = pd.DataFrame()
for res in F:
    date = res.partition("_")[-1][:-4]
    temp = pd.read_table(res)
    temp['date'] = date
    temp['type'] = temp['Label'].apply(
        lambda x: x.rpartition(":")[-1].partition("_")[0])
    temp['newname'] = temp['Label'].str.rpartition(
        ":").get(2).str.partition("_").get(2)
    temp['newname2'] = temp[['type',
                             'date',
                             'newname']].apply(lambda x: '_'.join(x), axis=1)
    P = P.append(temp, ignore_index=True)

RFP = P[P.Label.str.contains("RFP")].reset_index(drop=True)
GFP = P[P.Label.str.contains("GFP")].reset_index(drop=True)
meta = RFP.join(GFP, lsuffix='_left', rsuffix='_right')
meta[['newname2_left', 'Mean_left', 'Mean_right']].to_pickle('fileMetas.pkl')

D = meta[['newname2_left', 'Mean_left', 'Mean_right']].to_dict()
FPmean = zip(D['Mean_left'].values(), D['Mean_right'].values())
protodict = zip(D['newname2_left'].values(), FPmean)
realdict = {}

for key, val in protodict:
    realdict[key] = val

with open('fileMetas.pkl', 'wb') as output:
    pickle.dump(realdict, output)
