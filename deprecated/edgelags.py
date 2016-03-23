# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 00:07:58 2015
Output edge lags data into pkl for mungedata dataset
@author: sweel
"""
import matplotlib.pyplot as plt
import os
import cPickle as pickle
import seaborn as sns
import pandas as pd
sns.set_context("talk")
plt.close('all')


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


def lagedge(data, cols):
    """Take an edge data per cell and compute the lag=1
    This is done for all cell names in cols and the mean lag for that edge
    appended to the dat DataFrame for output.

    Parameters
    ----------
        data:
            DataFrame consisting of columns = cells, each row in the df
            is an edge list containing the values of interest
        cols:
            name of the cells

    Returns
    -------
        Series data of lag=1 list with index = cell names
    """
    lage = {}
    srtcol = sorted(cols)
    for cell in sorted(cols):
        df = data[cell]  # Series of edge data per cell
        dfedge = pd.DataFrame(
            {i: pd.Series(df[i]) for i in df.index})\
            .dropna(axis=1, how='all')

        lage[cell] = list(calcdlag(dfedge, 1).mean())
    dat = pd.Series([lage[cell] for cell in srtcol], index=srtcol)
    return dat

# =============================================================================
#    Data input
# =============================================================================
# pylint: disable=C0103
dirlist = []
# pylint: disable=C0103
for root, dirs, files in os.walk(os.getcwd()):
    for f in dirs:
        if f.startswith('YP'):
            dirlist.append(
                os.path.join(root, f))

DYL = pd.DataFrame()
FRAMES = []
for media in dirlist[:]:
    labs = media[-3:]
    print'\nNow on %s' % labs + "\n" + "="*79
    #   make sure the pkl file below exists, run MakeInputForLags.py otherwise
    with open('%s_lags.pkl' % labs, 'rb') as inpt:
        (randNDY, randUDY, Norm, NormPermute) = pickle.load(inpt)

    framesDY = [pd.Series([edge for edge in Norm[cells]],
                          name=cells)
                for cells in Norm.keys()]

    dfDY = pd.concat(framesDY, axis=1)
    cols1 = dfDY.columns

    FRAMES.append(lagedge(dfDY, cols1))
    print'DY complete'

dat = pd.concat(FRAMES)
with open('lagedges.pkl', 'wb') as OUT:
    pickle.dump(dat, OUT)
