# -*- coding: utf-8 -*-
"""
Created on Wed Jul 08 12:25:50 2015
Plot Autocorrelation curves by population
@author: sweel
"""
import numpy as np
from itertools import islice, izip_longest
from collections import defaultdict


def acf(vec):
    """Statistical method of autocorrelation using numpy corrcoef
    """
    length = len(vec)
    cor = [np.corrcoef(vec[:-i], vec[i:])[0, 1]
           for i in range(1, length)]
    return np.array([1] + cor)


def autocorout(edgedata):
    """Make the autocorrelation coefficient for every edge in a cell

    Parameters
    ----------

    edgeData :
        array of edges points/data from one cell
    """
    X = []
    for j in edgedata:
        #        edgeCor = estimated_autocorrelation(j)
        edgecor = acf(j)
        #        X.append(edgeCor[0])
        X.append(edgecor)
    return([i for i in X if len(i)])


def ac_cell(points, shift):
    """Make the autocorrelation coefficient for every edge in a cell

    Parameters
    ----------
    cell :
        cell/file Name
    edgeData :
        array of edges points/data from one cell
    """
    z = np.zeros(shift)
    points1 = np.concatenate([z, points[:-len(z)]])
    result = np.corrcoef(points, points1)[0, 1]
    #        X.append(edgecor)
    return result


def psd(cell, edgedata, T):
    '''return the power spectrum counts (ps) as function of spatial freq u
    '''
    X = []
    Y = []
    for j in edgedata[cell]:
        if len(j) >= T:
            #    Z = np.hstack(Y[cell])
            data = j - np.mean(j)
            ps = (np.abs(np.fft.rfft(data))**2)/len(j)
            time_step = 1
            u = np.fft.rfftfreq(data.size, time_step)
            idx = np.argsort(u)
            X.append(u[idx])
            Y.append(ps[idx])
    return X, Y
