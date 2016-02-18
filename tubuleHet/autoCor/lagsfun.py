# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 00:02:35 2015
Functions for calculating the lags of an edge
@author: sweel
"""
import pandas as pd


def calcdlag(df, lag):
    """Calc the lags in an edge or vector as a delta of two intensity points
    separated by a lag value and return the absolute value of that edge
    """
    delta = df.shift(lag) - df
    delta = delta.ix[lag:]
    delta = delta.dropna(how='all', axis=1)
#    rms = dv.apply(lambda x: np.sqrt(x**2))
    abso = delta.abs()
    return abso


def iterlagspd(data, cols, labs):
    """Take an edge data per cell and compute the lags in lagList using apply
    and calclagPD, and output in  wideform Dataframe with a 'cat' label=labs.
    This is done for all cell names in cols and the result of each iteration
    appended to the dat DataFrame for output.

    Parameters
    ----------
        data:
            DataFrame consisting of columns = cells, each row in the df
            is an edge list containing the values of interest
        cols:
            name of the cells
        labs:
            media type of the cells

    Returns
    -------
        Dataframe in wideform  of lag values with column name = lag
    """
    dat = pd.DataFrame()
    laglist = [1, 5, 10, 15, 20]
    for cell in sorted(cols):
        df = data[cell]  # Series of edge data per cell
        dfedge = pd.DataFrame(
            {i: pd.Series(df[i]) for i in df.index})\
            .dropna(axis=1, how='all')

        listoflags = [calcdlag(dfedge, k).mean() for k in laglist]

        out = pd.concat([l for l in listoflags], axis=1)
        dat = dat.append(out, ignore_index=True)
    dat.rename(columns=lambda x: laglist[x], inplace=True)
    dat['cat'] = labs
    return dat


#def calclagpd(col, lags=1):
#    """Function for calculating the lags in a vector as a delta of
#    two intensity points separated by a lag value and return the absolute
#    value of that edge as a Series of lag values. default kw lag value =1
#
#    """
#    shifted = col.shift(lags)
#    delta = shifted - col
#    delta = delta.ix[lags:]
#    delta = delta.dropna()
#    delta = delta.abs()
#    return delta


#def iterLags(df, cols):
#    """Iterate thru the col list of media type and calcualted the lags for
#    a range of lag values k
#
#    Returns
#    -------
#        Default dic of lag values indexed by cell name and lag value
#    """
#    out = defaultdict(dict)
#    for cell in cols:
#        df1 = df[cell]  # Series of edge data per cell
#        df2 = pd.DataFrame(
#            dict([(k, pd.Series(v)) for k, v
#                 in df1.iteritems()]))
#        for k in [1, 5, 10, 15, 20]:
#            out[cell][k] = []
#            abso, delta = calcDlag(df2, k)
#            out[cell][k].append(abso.mean())
#    return (out)
#
#
#def stackMeans(df, cellList, labs):
#    """Create a wide form dataframe with columns consisting of lags values
#    with a media type category columns
#
#    Parameters
#    ----------
#    df:
#        dataframe output of iterLags
#    cellList:
#        list of cell names
#    labs:
#        media type of population
#
#
#    """
#    dat = pd.DataFrame(index=np.arange(10000))
#    for k in [1, 5, 10, 15, 20]:
#        frames = [pd.Series(
#            df[cell][k][0]) for cell in cellList]
#        temp = pd.concat(frames)
#        temp = temp.reset_index(drop=True)
#        dat[k] = pd.Series(temp.iloc[:], index=dat.index)
#    dat = dat.dropna(how='all', subset=dat.columns)
#    dat['cat'] = labs
#    return(dat)
