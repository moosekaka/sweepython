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


def autocorout(cell, edgedata):
    """Make the autocorrelation coefficient for every edge in a cell

    Parameters
    ----------
    cell :
        cell/file Name
    edgeData :
        array of edges points/data from one cell
    """
    X = []
    for j in edgedata[cell]:
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



def plotauto2(autodata):
    """
    Returns the autocorrelation coefficient for a population
    with various minimum edge lengths given by minL and maxL and
    return as a default dictionary with index f=media type  and
    index b= minimum edge lengths

    Parameters
    ----------
    autoData :
        Dictionary containing list of outputs from autoCorOut:

        autoData[media_type][cell_name] = autoCorOut(cell,edgeData)
    """

    minl = 10
    maxl = 45
    celltemp = {}

    media = sorted(autodata.keys())
    celltemp = defaultdict(dict)

    for mem in media:

        for _, b in enumerate(range(minl, maxl, 10)):

            automedia = autodata[mem].values()

            X = np.ravel(
                [el for lis in automedia for el in lis[0]
                 if len(el) > b])
            celltemp[mem][b] = []
            temp = []

            for i in islice(izip_longest(*X, fillvalue=.0), 0, 15):
                temp.append(np.mean(i))

            celltemp[mem][b].append(temp)
    return celltemp


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
