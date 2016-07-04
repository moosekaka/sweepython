# -*- coding: utf-8 -*-
"""
Created on Wed Jul 08 12:25:50 2015
Plot Autocorrelation curves by population
@author: sweel
"""
import numpy as np
import pandas as pd
# pylint: disable=C0103
# pylint: disable=R0204

def acf(vec):
    """Statistical method of autocorrelation using numpy corrcoef
    https://en.wikipedia.org/wiki/Autocorrelation
    """
    L = len(vec)
    cor = [np.corrcoef(vec[:L-i], vec[i:L])[0, 1]
           for i in range(1, L)]
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
        edgecor = acf(j)  # estimated_autocorrelation(j)
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


def psd(edgedata, threshold_length):
    '''return the power spectrum counts (ps) as function of spatial freq u
    '''
    X = []
    Y = []
    for j in edgedata:
        if len(j) >= threshold_length:
            #    Z = np.hstack(Y[cell])
            data = j - np.mean(j)
            ps = (np.abs(np.fft.rfft(data))**2)/len(j)
            time_step = 1
            u = np.fft.rfftfreq(data.size, time_step)
            idx = np.argsort(u)
            X.append(u[idx])
            Y.append(ps[idx])
    return X, Y


def conv_to_pd(autodata):
    """
    Returns flattened dataframe of PSD values for each cell (columns)

    Parameters
    ----------

    autodata:
        Dictionary of PSD by cells (thresholded to a min. length)

    Returns
    -------

    pdx, pdy:
        Dataframe of the power spectrum power (y) at frequencies (x)

    """
    x1 = pd.concat([pd.Series(autodata[k][0][0], name=k) for k
                    in autodata.keys()], axis=1)
    y1 = pd.concat([pd.Series(autodata[k][0][1], name=k) for k
                    in autodata.keys()], axis=1)
    x_ls = x1.stack()
    y_ls = y1.stack()
    x_ls.reset_index(drop=True, inplace=True)
    y_ls.reset_index(drop=True, inplace=True)
    pdx = pd.DataFrame({i: pd.Series(x_ls.ix[i]) for i in x_ls.index})
    pdy = pd.DataFrame({i: pd.Series(y_ls.ix[i]) for i in y_ls.index})
    return pdx, pdy


def tidy_psd(pdx, pdy, bins, typelabel):
    """
    Compresses the flat dataframe into a longform tidydata format

    Parameters
    ----------

    pdx, pdy:
        Dataframe of the power spectrum power (y) at frequencies (x)
    bins:
        the histogram bins/ half open intervals of edgelens
    typelabel:
        categorical labels for each type

    Returns
    -------

    tidypd:
        Dataframe in longform with var psd and categorical var typelabel
        and id variables of frequencies (u)

    """
    temp = []
    tidypd = pd.DataFrame()
    labs = np.round(bins[1:], 2)
    for edge in pdx.columns:
        mean_y = pdy.ix[:, edge].groupby(pd.cut(pdx.ix[:, edge],
                                                bins,
                                                labels=labs)).mean()
        temp.append(mean_y)
    psdframe = pd.concat(temp, axis=1)
    psdframe['u'] = labs
    psdtidy = pd.melt(psdframe,
                      id_vars='u',
                      var_name='lineid',
                      value_name='psd')
    psdtidy['type'] = typelabel
    tidypd = tidypd.append(psdtidy, ignore_index=True)
    return tidypd


def binedges(coldata, interval):
    """
    Bins the edge lists according to length of edge array

    Parameters
    ----------

    coldata:
        column from dataframe, each row contains vector of autocor coeffs for
        and edge
    interval:
        the histogram bins/ half open intervals of edgelens
    """
    coldata = coldata.dropna()
    edgelen = coldata.apply(len)
    edgs = pd.DataFrame({'auto_cor': coldata, 'len': edgelen})
    edgs['cut'] = pd.cut(edgs.len,
                         interval,
                         labels=interval[:-1],
                         include_lowest=True,
                         right=False)
    edgs = edgs.dropna()  # drops edges < min len
    return edgs


def lagdist(edgebinned, thresh, label):
    """
    Returns an array of lag distances (k) with the autocor coef at that k

    Parameters
    ----------

    edgebinned:
        edges with autocorr binned by lengths
    thresh:
        cutoff (lower bound) threshold length
    labels:
        categorical labels for each type
    """
    binned_edges = edgebinned.loc[edgebinned.cut == thresh]
    dftemp = pd.DataFrame({i: pd.Series(j) for i, j
                           in binned_edges.ix[:, 'auto_cor'].iteritems()})

    dftemp = dftemp.stack().reset_index(0)
    dftemp.columns = ['lag', 'auto_cor']
    dftemp['thresh'] = thresh
    dftemp['type'] = label
    return dftemp


def calcdlag(df, lag):
    """Calc the lags in an edge or vector as a delta of two intensity points
    separated by a lag value and return the absolute value of that edge
    """
    delta = df.shift(lag) - df
    delta = delta.ix[lag:]
    delta = delta.dropna(how='all', axis=1)
    abso = delta.abs()
    return abso


def iterlagspd(data, labs, lagdist=None):
    """Take an edge data per cell and compute the lags in lagList using apply
    and calclagPD, and output in  wideform Dataframe with a 'cat' label=labs.
    This is done for all cell names in cols and the result of each iteration
    appended to the dat DataFrame for output.

    Parameters
    ----------
    data:
        DataFrame consisting of columns = cells, each row in the df
        is an edge list containing the values of interest
    labs:
        media type of the cells
    lagdist:
        list of lag distances, defaults [1, 5, 10, 15, 20]

    Returns
    -------
        Dataframe in wideform of lag values with column name = lag
    """
    if lagdist is None:
        lagdist = [1, 5, 10, 15, 20]  # default values
    dat = pd.DataFrame()
    data = np.squeeze(data)

    dfedge = pd.DataFrame(
        {row: pd.Series(edge) for row, edge in enumerate(data)})

    dfedge.dropna(how='all', inplace=True)
    listoflags = [calcdlag(dfedge, k).mean() for k in lagdist]

    out = pd.concat([l for l in listoflags], axis=1)
    dat = dat.append(out, ignore_index=True)
    dat.rename(columns=lambda x: lagdist[x], inplace=True)
    dat['cat'] = labs
    return dat
