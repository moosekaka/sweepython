# -*- coding: utf-8 -*-
"""
       module to analyze mom bud asymmetry
"""

import matplotlib.pyplot as plt
import os
import numpy as np
from mayavi import mlab
import cPickle as pickle
import fnmatch
import pandas as pd
import seaborn as sns
from collections import defaultdict
import scipy.stats as sp
plt.rcParams['font.family'] = 'DejaVu Sans'

# pylint: disable=C0103
plt.close('all')
mlab.close(all=True)
vtkF = defaultdict(dict)
mombud = defaultdict(dict)

# =============================================================================
# filelist and graph list
# =============================================================================
with open('mombudtrans.pkl', 'rb') as inpt:
    dfmb = pickle.load(inpt)  # has columns base, neck, tip, media, bud, mom

for root, dirs, files in os.walk(os.getcwd()):
    for i in files:
        if fnmatch.fnmatch(i, '*.vtk'):
            media = root.rsplit('\\', 1)[1]
            vtkF[media][i[:-4]] = os.path.join(root, i)

filekeys = {item: vtkF[media][item] for media
            in sorted(vtkF.keys()) for item
            in sorted(vtkF[media].keys())}


# =============================================================================
# compute z positions dyRaw
# =============================================================================
cellall = pd.DataFrame(columns=['mom', 'bud'])
#binsbudprog = np.linspace(0.,1.,6)
binsbudprog = np.array([0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,2.])
binsaxis = np.linspace(0.,1.,6)
binsaxisbig=np.linspace(0,1.,11)
#bins = np.array([0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.])
binsvolbud = np.linspace(0,40,5)
binsvolmom = np.array([ 0, 30,40, 80.])
cellposmom = pd.DataFrame()
cellposbud = pd.DataFrame()

for _, key in enumerate(sorted(filekeys.keys())[:]):
    src = mlab.pipeline.open(filekeys[key],
                             figure=False)

    data = src.outputs[0]
    npx = data.points.to_array()  #is a column vec of R^3 (coordinates in the skel)
    xind = npx[:, 0].argsort()  # indices of npx that would sort npx according to the x-axis
    dy = data.point_data.get_array(u"DY_minmax").to_array()
    y = dy[xind]  # sorted Δψ values by x position
    #  individual skeletons xyz and Δψ
    cell = pd.DataFrame({'x': npx[:, 0][xind],
                         'DY': dy[xind]})
    xn, yn, zn = dfmb.ix[key, 'neck']
    xb, yb, zb = dfmb.ix[key, 'base']
    xt, yt, zt = dfmb.ix[key, 'tip']
#    P = (xt - xn) / (xn - xb)
    xn_scaled = (xn - cell.ix[0, 'x']) / (xt - xb)
    cell['posxcell'] = (cell.ix[:, 'x'] - cell.ix[0, 'x']) / (xt - xb)
#    cell['posx'] = (cell.ix[:, 'x'] - xn)

#    cell[u"DY neck"] = cell.ix[budneckreg].DY.mean()
    cell['type'] = ''
    cell.loc[cell.x > xn, ['type']] = 'bud'
    cell.loc[cell.x <= xn, ['type']] = 'mom'
    cell.ix[cell.type=='bud','posx'] =(cell.ix[:, 'x']-xn)/(xt-xn)
    cell.ix[cell.type=='mom','posx'] =(cell.ix[:, 'x']-xb)/(xn-xb)
    budneckreg = ((cell.posx >= xn_scaled - .05) &
                  (cell.posx <= xn_scaled + .05)).nonzero()
    cell.reset_index(drop=True, inplace=True)
    cell['pos'] = cell.index
    X = cell.groupby('type').median().DY
    X.name = key


    # scale Δψ to min-max of the GRADIENT within mom/bud
    binsmombud = pd.cut(cell.ix[:, 'posx'],
                        binsaxis,
                        labels=binsaxis[1:])
    cell['binposx'] = binsmombud
     # scale Δψ to min-max of the GRADIENT within cell
    binscellaxis = pd.cut(cell.ix[:, 'posxcell'],
                          binsaxisbig,
                          labels=binsaxisbig[1:])
    cell['binposxcell'] = binscellaxis

    Ymom = cell.ix[cell['type']=='mom'].groupby('binposx').DY.mean()
    Ybud = cell.ix[cell['type']=='bud'].groupby('binposx').DY.mean()
    Ycell = cell.groupby('binposxcell').DY.mean()
    Zmom = (Ymom-Ycell.min())/(Ycell.max()-Ycell.min())
    Zbud = (Ybud-Ycell.min())/(Ycell.max()-Ycell.min())
#    Z = (Y-W.min())/(W.max()-W.min())
    Zbud.name = key
    Zmom.name = key


    cellall = cellall.append(X)
    cellposbud = cellposbud.append(Zbud)
    cellposmom = cellposmom.append(Zmom)

cellall['budvol'] = dfmb.bud
cellall['momvol'] = dfmb.mom
cellall = cellall.reset_index()
cellall['type'] = cellall['index'].apply(lambda x: x[:3])
cellposbud = cellposbud.reset_index()
cellposmom = cellposmom.reset_index()
cellposbud = pd.concat([cellposbud, cellall.ix[:, ['type', 'neck']]], axis=1)
cellposmom = pd.concat([cellposmom, cellall.ix[:, ['type', 'neck']]], axis=1)

cellall['frac'] = cellall.ix[:, 'bud'] / cellall.ix[:, 'mom']
Q = cellall.groupby('type').quantile(.90)
cellall['q90'] = cellall.type.apply(lambda x: Q.ix[x].budvol)
gt90 = cellall[cellall['budvol']>cellall['q90']]
meangt90 = gt90.groupby('type').budvol.mean()
cellall['mean90']= cellall.type.apply(lambda x: meangt90.ix[x])
#  budvolratio is based on the mean of the largest 10% cells
#cellall['budvolratio'] = cellall.budvol / cellall.mean90
cellall['budvolratio'] = cellall.budvol / cellall.q90
#cellpos['budvolratio'] = cellall['budvolratio']
cellposbud['budvol'] = cellall['budvol']
cellposmom['momvol'] = cellall['momvol']


binbudvol = pd.cut(cellposbud.ix[:, 'budvol'], binsvolbud, labels=binsvolbud[1:])
binbudmom = pd.cut(cellposmom.ix[:, 'momvol'], binsvolmom, labels=binsvolmom[1:])
cellposbud['binvol'] = binbudvol
cellposmom['binvol'] = binbudmom

#cellall = cellall.dropna()
BIG = pd.melt(cellall,
              id_vars=['type'],
              value_vars=['frac'])
groups = BIG.groupby('type')
G = groups.apply(lambda x: x[x['value'] > 1])
a1 = G.groupby('type').value.count()
a2 = cellall.groupby('type').frac.count()
a1 / a2 * 1.

BIG2 = pd.melt(cellall,
               id_vars=['type'],
               value_vars=['mom', 'bud'])

A = pd.cut(cellall.ix[:, 'budvolratio'], binsbudprog, labels=binsbudprog[1:])
cellall['bin'] = A
cellall['binbudvol']=cellposbud['binvol']
rejectlist=cellposmom.ix[(np.asarray(cellposmom.momvol)>60) &  #reject super large cells
                         (cellposmom.type!='YPD'),'index']
cellall = cellall.ix[~cellall.ix[:,'index'].isin(rejectlist)]
cellposmom = cellposmom.ix[~cellposmom.ix[:,'index'].isin(rejectlist)]
cellposbud = cellposbud.ix[~cellposbud.ix[:,'index'].isin(rejectlist)]
# =============================================================================
# Progression of Δψ as move along the bud axis
# =============================================================================
sns.set_style('whitegrid')
bigbinsmom = pd.melt(cellposmom,
                  id_vars=['type','binvol'],
                  var_name='mom axis position',
                  value_name=r'$\Delta\Psi$ scaled gradient',
                  value_vars=binsaxis.tolist())
bigbinsmom = bigbinsmom.dropna()
bigbinsbud = pd.melt(cellposbud,
                  id_vars=['type','binvol'],
                  var_name='bud axis position',
                  value_name=r'$\Delta\Psi$ scaled gradient',
                  value_vars=binsaxis.tolist())
bigbinsbud = bigbinsbud.dropna()
#
with sns.plotting_context('talk', font_scale=1.1):
    plt.close('all')
    h = sns.FacetGrid(bigbinsmom,
                      col="type",
                      hue='type',
                      col_wrap=2)
    h = h.map(sns.pointplot,
              'mom axis position',
              r'$\Delta\Psi$ scaled gradient')
    h.set_xlabels('mom cell axis position')
    h.set(ylim=(0, 1.))
#
#    j = sns.FacetGrid(bigbinsbud,
#                      col="type",
#                      hue='type',
#                      col_wrap=2)
#    j = j.map(sns.pointplot,
#              'bud axis position',
#              r'$\Delta\Psi$ scaled gradient')
#    j.set_xlabels('bud axis position')
#    # Adjust the tick positions and labels
#    j.set(ylim=(0, 1.))
#
##    _, ax0 = plt.subplots(1, 1)
##    f = sns.pointplot(x='cell axis position',
##                      y=r'$\Delta\Psi$',
##                      hue='type',
##                      data=bigbins,
##                      ax=ax0)
##    f.get_legend().set_visible(False)
##    f.set_title(u" Δψ along cell axis")
##    f.set_xlabel("mom cell axis")
##    f.set_ylabel(u" gradien Δψ")
#
##    g = sns.FacetGrid(cellall, row="type",hue='type')
##    g = g.map(sns.distplot, "budvolratio")
#
#    e = sns.FacetGrid(cellall, row="type",hue='type')
#    e = e.map(plt.hist, "budvol")
#    g = sns.FacetGrid(cellall, row="type",hue='type')
#    g = g.map(plt.hist, "momvol")
##    for suba in g.axes:
##        for line in [ 30,40]:
##            suba[0].axvline(line)
#
#
#    #  facet the data into media and budvolumes
#    k = sns.FacetGrid(bigbinsmom,
#                      row="type",
#                      col="binvol",
#                      hue='type',
#                      col_order=binsvolmom[1:],
#                      )
#    k = k.map(sns.pointplot,
#              'mom axis position',
#              r'$\Delta\Psi$ scaled gradient')
#    k.set_xticklabels(fontsize=12)
#    k.set_xlabels("mom axis position")
#    k.set(ylim=(0, 1.))
#
    m = sns.FacetGrid(bigbinsbud,
                      row="type",
                      col="binvol",
                      hue='type',
                      col_order=binsvolbud[1:],
                      )
    m = m.map(sns.pointplot,
              'bud axis position',
              r'$\Delta\Psi$ scaled gradient')
    m.set_xticklabels(fontsize=14)
    m.set_xlabels("bud axis position")
    m.set(ylim=(0, 1.))

## =============================================================================
## frac Δψ as function of budratio
## =============================================================================
with sns.plotting_context('talk'):
    _, ax2 = plt.subplots(1, 1)
    h = sns.pointplot(x='bin',
                      y='frac',
                      hue='type',
                      data=cellall.dropna(),
                      ax=ax2)
    h.get_legend().set_visible(False)
    h.set_ylim(0, 3.)


    h.set_title(u"Δψ vs bud progression\n ")
    h.set_xlabel("bud progression")
    h.set_ylabel(u"Δψ bud/Δψ mom")

    p = sns.FacetGrid(cellall.dropna(), col="type", hue='type', col_wrap=2)
    p = p.map(sns.pointplot, 'bin', 'frac')

    q = sns.FacetGrid(cellall.dropna(),
                      col="binbudvol",
                      hue='type',
                      row='type')
    q = q.map(sns.pointplot,
              'bin',
              'frac')
    q.set_xlabels('bud progression')
## ============================================================================
## frac Δψ violinplots by media
## ============================================================================
#with sns.plotting_context('talk', font_scale=1.4):
#    _, ax1 = plt.subplots(1, 1)
#    g = sns.violinplot(x='type',
#                       y='value',
#                       hue='type',
#                       data=BIG,
#                       ax=ax1)
#    g.set_ylim(0, 3)
#    g.get_legend().set_visible(False)
#    g = sns.stripplot(x='type',
#                      y='value',
#                      hue='type',
#                      data=BIG,
#                      jitter=.15,
#                      size=5,
#                      ax=ax1)
#    g.set_ylim(0.4, 2.0)
#    g.get_legend().set_visible(False)
## =============================================================================
## violinplot mom vs bud Δψ scaled
## =============================================================================
#with sns.plotting_context('talk', font_scale=1.4):
#    _, ax3 = plt.subplots(1, 1)
#    h = sns.violinplot(x='type',
#                       y='value',
#                       hue='variable',
#                       data=BIG2,
#                       ax=ax3)
#    sns.stripplot(x='type',
#                  y='value',
#                  hue='variable',
#                  jitter=.15,
#                  size=4,
#                  data=BIG2,
#                  ax=ax3)
#    h.set_ylim(0, 1.)
#    h.get_legend().set_visible(False)
#
## =============================================================================
## frac Δψ as function of budvol
## =============================================================================
##with sns.plotting_context('talk', font_scale=1.4):
##    _, ax10 = plt.subplots(1, 1)
###    g = sns.FacetGrid(cellall.dropna(), col="type")
###    g = g.map(sns.regplot, "budvol", "frac")
##    datacell = cellall[cellall.bud <= bins2[-1]]
##    h = sns.pointplot(x='binbudratio',
##                      y='bud',
##                      hue='type',
##                      data=datacell.dropna(),
##                      ax=ax10)
##    h.get_legend().set_visible(False)
##
##    for i in ['YPD', 'YPE', 'YPL', 'YPR']:
##        data = cellall[(cellall.type == i) & (cellall.frac < 2)]
##        slope, _, r, p, _ = sp.linregress(data['budvol'],
##                                          data['frac'])
##        print 'slope= %6.4f r=%6.4f p=%6.4f' % (slope, r, p)
#
## =============================================================================
## Dy as budneckregion and budratio
## =============================================================================
##with sns.plotting_context('talk', font_scale=1.4):
##    _, ax2 = plt.subplots(1, 1)
##    h = sns.pointplot(x='bin',
##                      y='DYneck',
##                      hue='type',
##                      data=cellall.dropna(),
##                      ax=ax2)
##    h.get_legend().set_visible(False)
#
#
## with sns.plotting_context('talk', font_scale=1.4):
##    _, ax1 = plt.subplots(1, 1)
##    h = sns.pointplot(x='posx',
##                      y='DY',
##                      ci=None,
##                      markers='o',
##                      join=False,
##                      hue='type',
##                      data=cell,
##                      size=1,
##                      ax=ax1)
##    h.get_legend().set_visible(False)
##    h.set_xticks(np.linspace(cell.pos.min(), cell.pos.max(),11))
##    h.set_xticklabels(np.arange(0, 1.1 ,.1))
#
##==============================================================================
## budratio
##==============================================================================
##c2 = cellall.drop(cellall.index[[5, 15, 63, 46]])
##slope, _, r, p, std_err = sp.linregress(c2.ix[:, 'budratio'],
##                                        c2.ix[:, 'neck'])
## with sns.plotting_context('talk', font_scale=1.4):
##    _, ax5 = plt.subplots(1, 1)
##    h = sns.regplot(x='budratio',
##                      y='neck',
##                      data=c2[c2.neck>0.505],
##                      ax=ax5)
##    h.get_legend().set_visible(False)
