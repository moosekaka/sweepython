# -*- coding: utf-8 -*-
"""
@author: sweel
process O2 data
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
import pandas as pd
import fnmatch
from collections import defaultdict
from networkHet.mungedata import MungeDataFuncs as md
import seaborn as sns
import math
import cPickle as pickle
# pylint: disable=C0103
D = defaultdict(dict)


def bootO2(X, N):
    '''boostrap the CI for the O2slope at OD=0.5
    '''
    result = []
    result2 = []
    result3 = []
    result4 = []
    result5 = []
    lenB = len(X)
#    lenB1 = len(Y)
    for _ in range(N):
        rX = X.sample(n=lenB, replace=True)
#        rY = Y.sample(n=lenB1, replace=True)
        A, B, _, _, _ = sp.linregress(rX.ix[:, r'$OD_{600}$'],
                                      rX.ix[:, 'O2slope'])
        ocr = A * .5 + B
        idx = X.index.tolist()
        cell_num = X.ix[idx[0], 'count'] / 2
        cellmass = X.ix[idx[0], 'mass'] / 2
#        cell_num = rY['countOD1'].median() / 2
        # 1E-6 * (1E-6)^2 gives micron^3, divide by 1E18 to be in picom^3
        mitovol = X.ix[idx[0], 'mitovol'] * (.15**2) * math.pi
        ocrm = ocr / cell_num / mitovol / 10  # factor to make consistent
        ocrc = ocr / cell_num
        ocrw = ocr / cellmass
        result.append(A * .5 + B)
        result2.append(ocrm)
        result3.append(cell_num)
        result4.append(ocrc)
        result5.append(ocrw)
    return result, result2, result3, result4, result5

# =============================================================================
# plot the O2 rates data points
# =============================================================================
L = range(7, 26)
L1 = range(0, 4)
L2 = np.concatenate([L1, L], 0)
plt.close('all')

#dfcount = pd.read_csv('cellcount1.csv',
#                      names=['countOD1', 'type'],
#                      header=None)

df0 = pd.read_csv('RawData.csv')
df0[r'$OD_{600}$'] = 0
df0['O2slope'] = 0
for i in df0.Name.tolist():
    j = i.partition(r'$OD_{600}$')[0].partition('_')
    df0.ix[df0.Name == i, r'$OD_{600}$'] = float(j[2][:3]) / 1000

for root, dirs, files in os.walk(os.getcwd()):
    for i in files:
        if fnmatch.fnmatch(i, '[Yy]*.csv'):
            D[i] = os.path.join(root, i)

for filename in df0.Name.tolist():
    mark = df0[df0.Name == filename].index.tolist()
    df = pd.read_csv(D[filename],
                     header=None,
                     names=['time', 'O2', 'Marker'],
                     skiprows=L2)
    t0 = df[df.Marker == 'Start'].index.tolist()
    t1 = df[df.Marker == 'Stop'].index.tolist()
    df1 = df.ix[t0[0]:t1[0], :2]  # first filter to rest to start/stop time
    df1 = df1.astype(float)  # make sure the values are floats
    linear_t0 = df0.ix[mark[0], 'T0']
    linear_t1 = df0.ix[mark[0], 'T1']
    df2 = df1[(df1.time >= linear_t0) & (df1.time <= linear_t1)]

    slope, _, _, _, _ = sp.linregress(df2['time'],
                                      df2['O2'])
    df0.ix[df0.Name == filename, 'O2slope'] = -slope

with sns.plotting_context('talk', font_scale=1.5):
    g = sns.lmplot(x=r'$OD_{600}$',
                   y='O2slope',
                   col="type",
                   hue='type',
                   hue_order=['YPD', 'YPE', 'YPL', 'YPR'],
                   data=df0,
                   col_wrap=2,
                   size=5,
                   scatter_kws={"s": 25, "alpha": 1})
    for ax in g.axes:
        if ax.get_ylabel() != '':
            ax.set_ylabel('OCR (nmol/ml/s)')

# =============================================================================
# bootstrapped CI for the mito OCR
# =============================================================================
df = pd.DataFrame()
df['type'] = ''
F = {}
G = {}
H= {}
cn = {}
mvl = {}
ocr = {}
media = ['YPD', 'YPE', 'YPL', 'YPR']
for run in media:
    OCR, OCRmito, cellnum, OCRcell, OCRmass = bootO2(df0.ix[df0.type == run], 100)
#                                       dfcount.ix[dfcount.type == run],
#                                       1000)
    G[run] = OCRcell  # this is for OCR norm by cell number
    F[run] = OCRmito
    H[run] = OCRmass
    cn[run] = cellnum
    ocr[run] = OCR

A = pd.DataFrame({key: pd.Series(F[key]) for key in F.keys()})
B = A.stack().reset_index()
C = B.drop('level_0', 1)
C.columns = ['type', 'OCRmito']

D1 = pd.DataFrame({key: pd.Series(G[key]) for key in G.keys()})
E1 = D1.stack().reset_index()
F1 = E1.drop('level_0', 1)
F1.columns = ['type', 'OCRcell']

J = pd.DataFrame({key: pd.Series(H[key]) for key in H.keys()})
K = J.stack().reset_index()
M = K.drop('level_0', 1)
M.columns = ['type', 'OCRmass']

Q1 = C.groupby('type').quantile(0.025).values.flatten()
Q2 = C.groupby('type').quantile(0.5).values.flatten()
Q3 = C.groupby('type').quantile(0.975).values.flatten()
Q1n = F1.groupby('type').quantile(0.025).values.flatten()
Q2n = F1.groupby('type').quantile(0.5).values.flatten()
Q3n = F1.groupby('type').quantile(0.975).values.flatten()
Q1m = M.groupby('type').quantile(0.025).values.flatten()
Q2m = M.groupby('type').quantile(0.5).values.flatten()
Q3m = M.groupby('type').quantile(0.975).values.flatten()
with sns.plotting_context('talk', font_scale=1.5):
    fig, (ax1, ax2, ax3) = plt.subplots(3,1)
    fig.set_size_inches(8.5, 14, forward=True)
    g = sns.barplot(x='type',
                    y='OCRmito',
                    ax=ax1,
                    estimator=np.median,
                    ci=None,
                    ecolor=[.25, .25, .25],
                    data=C,
                    yerr=[Q2-Q1, Q3-Q2])

    g.set_ylabel(r'OCR per mito vol $\mu m^{3}$')
    ax1.set_xlabel('')
    ax1.set_yticks(np.arange(0,.015,.0025))
    for xlab in g.axes.xaxis.get_ticklabels():
        xlab.set_visible(False)

#with sns.plotting_context('talk', font_scale=1.5):
#    fig, ax2 = plt.subplots(1, 1)
    g = sns.barplot(x='type',
                    y='OCRcell',
                    ax=ax3,
                    estimator=np.median,
                    ci=None,
                    ecolor=[.25, .25, .25],
                    data=F1,
                    yerr=[Q2n-Q1n, Q3n-Q2n])
    g.set_ylabel(r'OCR per cell$ ( \times  10^{7})$ ')
    ax3.set_xlabel('')
    for xlab in g.axes.xaxis.get_ticklabels():
       xlab.set_visible(False)

#with sns.plotting_context('talk', font_scale=1.5):
#    fig, ax3 = plt.subplots(1, 1)
    g = sns.barplot(x='type',
                    y='OCRmass',
                    ax=ax2,
                    estimator=np.median,
                    ci=None,
                    ecolor=[.25, .25, .25],
                    data=M,
                    yerr=[Q2m-Q1m, Q3m-Q2m])
    g.set_ylabel(r'OCR per mg/ml ')
    ax2.set_xlabel('')
    for xlab in g.axes.xaxis.get_ticklabels():
       xlab.set_visible(False)


C.groupby('type').median()
#F.groupby('type').median()
res = md.multiple_test(C, 'OCRmito', 'type')
print res[2]

res = md.multiple_test(F1, 'OCRcell', 'type')
print res[2]

with open('o2data.pkl', 'wb') as outp:
    pickle.dump(C,outp)

#with sns.plotting_context('talk', font_scale=1.5):
#    fig, ax2 = plt.subplots(1, 1)
#    g = sns.barplot(x='type',
#                    y='countOD1',
#                    ax=ax2,
#                    estimator=np.median,
#                    ci=95,
#                    ecolor=[.25, .25, .25],
#                    data=dfcount)










# ONLY USED FOR PROCESSING RAW DATA
# ===========================================================================
#  plot out the raw O2 fitting to determine the linear region
# =============================================================================
#L = range(7, 26)
#L1 = range(0, 4)
#L2 = np.concatenate([L1, L], 0)
# plt.close('all')
#df0 = pd.read_csv('RawData.csv')
#
# for root, dirs, files in os.walk(os.getcwd()):
#    for f in files:
#        if fnmatch.fnmatch(f, 'YPR_520*.csv'):
#            print f
#            df = pd.read_csv(os.path.join(root, f),
#                             header=None,
#                             names=['time', 'O2', 'Marker'],
#                             skiprows=L2)
#
#            t0 = df[df.Marker == 'Start'].index.tolist()
#            t1 = df[df.Marker == 'Stop'].index.tolist()
#            df1 = df.ix[t0[0]:t1[0], :2]
#            df1 = df1.astype(float)
#            fig, ax1 =plt.subplots(1,1)
#            fig.canvas.set_window_title(os.path.join(root,f[:-4]))
#            df1.plot(x='time', y='O2',ax=ax1)
##
# ==============================================================================
# get start and end times linear region for prev marked files (T0,T1)
# ==============================================================================
# for root, dirs, files in os.walk(os.getcwd()):
#    for i in files:
#        if fnmatch.fnmatch(i, '*.csv'):
#            df=pd.read_csv(os.path.join(root, i),
#                           header=None,
#                           names=['time', 'O2', 'Marker'],
#                           skiprows=L2)
#            t0 = df[df.Marker == 'T0'].index.tolist()
#            t1 = df[df.Marker == 'T1'].index.tolist()
#            T0=float(df.ix[t0[0],'time'])
#            T1=float(df.ix[t1[0],'time'])
#            print '%s: %s,%s,%s' % (i,T0,T1,T1-T0)
