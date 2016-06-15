# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 00:50:45 2016

@author: sweel_Rafelski
"""
import cPickle as pickle
from collections import defaultdict
# pylint: disable=C0103

# _____________________________________________________________________________
if __name__ == '__main__':
    with open('celldata.pkl', 'rb') as inp:
        vtkdf = pickle.load(inp)

    tempdic = defaultdict(dict)
    dd = 0.3
    for k in sorted(vtkdf.keys()):
        print "now on cell {}".format(k)
        neck_position = vtkdf[k]['neckpos']
        cell = vtkdf[k]['df']
        tempdic[k]['bud'] = {}
        tempdic[k]['mom'] = {}
        c2 = cell.reset_index(level='name')
        for i in range(50):
            cellsample = cell.sample(n=cell.shape[0])
            cs = cellsample[['DY', 'DY_abs']].reset_index(level='name')
            c3 = c2[['x', 'type']].merge(cs,
                                         left_index=True,
                                         right_index=True,
                                         indicator=True)
            tempdic[k]['bud'][i] = c3.loc[
                (c3.x >= neck_position) &
                (c3.x < (neck_position + dd))].DY.mean()
            tempdic[k]['mom'][i] = c3.loc[
                (c3.x < neck_position) &
                (c3.x >= (neck_position - dd))].DY.mean()

    with open('boot.pkl', 'wb') as out:
        pickle.dump(tempdic, out)
