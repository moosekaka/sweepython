# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 00:50:45 2016

@author: sweel_Rafelski
"""
from collections import defaultdict
import cPickle as pickle
from mombud.functions import vtk_mbfuncs as vf
# pylint: disable=C0103


def bootNeck(vtkdf, dd=0.3, num_runs=100, save=False, **kwargs):
    """
    Parameters
    ----------

    vtkdf : DataFrame
        Ind cell data datafrmes

    dd : Float
        distance from neck

    num_runs : Int
        num of sample runs for bootstrap

    Kwargs
    ------

    outdatapath : Str
        filepath for saving bootstrap results

    """
    merge = defaultdict(dict)
    for k in sorted(vtkdf.keys()):
        print "now on cell {}".format(k)
        neck_position = vtkdf[k]['neckpos']
        cell = vtkdf[k]['df']
        merge[k]['bud'] = {}
        merge[k]['mom'] = {}
        cleft = cell.reset_index(level='name')
        for i in range(num_runs):
            cright = cell.sample(n=cell.shape[0])
            cright = cright[['DY', 'DY_abs']].reset_index(level='name')
            c3 = cleft[['x', 'type']].merge(cright,
                                            left_index=True,
                                            right_index=True,
                                            indicator=True)
            merge[k]['bud'][i] = c3.loc[
                (c3.x >= neck_position) &
                (c3.x < (neck_position + dd))].DY.mean()
            merge[k]['mom'][i] = c3.loc[
                (c3.x < neck_position) &
                (c3.x >= (neck_position - dd))].DY.mean()

    if save:
        fpath = kwargs.get('outdatapath')
        with open(fpath, 'wb') as out:
            pickle.dump(merge, out)
# _____________________________________________________________________________
if __name__ == '__main__':
    args = {'inpdatpath': 'celldata.pkl',
            'outdatapath': 'bootneck.pkl'}
    data = vf.wrapper(**args)
    bootNeck(data, dd=0.9, num_runs=1, save=False, **args)
