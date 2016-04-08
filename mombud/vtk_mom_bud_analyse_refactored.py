# -*- coding: utf-8 -*-
"""
module to analyze mom bud asymmetry
"""
import sys
import os
import os.path as op
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
from mayavi import mlab
import pandas as pd
import seaborn as sns
from mombud.vtk_viz import vtk_mbfuncs as vf
import wrappers as wr
# pylint: disable=C0103
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.close('all')
mlab.close(all=True)
datadir = op.join(os.getcwd(), 'data', 'transformedData')

# filelist and graph list
with open(op.join(datadir, 'mombudtrans.pkl'), 'rb') as inpt:
    dfmb = pickle.load(inpt)  # has columns base, neck, tip, media, bud, mom

try:
    filext = "*vtk"
    vtkF = wr.ddwalk(datadir, filext, stop=-4)
except LookupError:
    sys.exit("error filetypes %s not found in %s" % (filext, datadir))

filekeys = {item: vtkF[media][item] for media
            in sorted(vtkF.keys()) for item
            in sorted(vtkF[media].keys())}

# compute 'z' positions  dyRaw
cellall = pd.DataFrame(columns=['mom', 'bud'])
cellposmom = pd.DataFrame()
cellposbud = pd.DataFrame()
neckregion = pd.DataFrame()
binsbudprog = np.r_[np.arange(0, 1.1, .1), 2]
binsaxis = np.linspace(0., 1., 6)
binsaxisbig = np.linspace(0, 1., 11)
binsvolbud = np.linspace(0, 40, 5)
binsvolmom = np.array([0, 30, 40, 80.])

for key in sorted(filekeys)[:]:
    try:
        cell = vf.safecall(key, filekeys, dfmb)
    except LookupError:
        sys.exit("%s not in filekeys" % key)
    # pos along x-axis for inidivual mom or bud cell
    cell['binposx'] = vf.bincell(cell, 'posx', binsaxis)
    # scaled Δψ to min-max of the GRADIENT within mom/bud
    Xmom = cell.ix[cell['type'] == 'mom'].groupby('binposx').DY.mean()
    Xbud = cell.ix[cell['type'] == 'bud'].groupby('binposx').DY.mean()

    # pos along x-axis for the whole cell
    cell['binposxcell'] = vf.bincell(cell, 'posxcell', binsaxisbig)
    # scale Δψ to min-max of the GRADIENT within cell
    DYcell = cell.groupby('binposxcell').DY.mean()

    Xmom = vf.scaleminmax(Xmom, DYcell)
    Xbud = vf.scaleminmax(Xbud, DYcell)
    Xbud.name = key
    Xmom.name = key
    medianDY = cell.groupby('type').median().DY
    medianDY.name = key
    cellall = cellall.append(medianDY)
    cellposbud = cellposbud.append(Xbud)
    cellposmom = cellposmom.append(Xmom)
    # temp dict of mean Δψ at +- range of dist from budneck
    tempdic = {dist: vf.neckDY(cell, cell.neckpos, dist)
               for dist in [.15, .3, .5]}
    temp = pd.DataFrame({'bud': pd.Series({k: tempdic[k][0] for k in tempdic}),
                         'mom': pd.Series({k: tempdic[k][1] for k in tempdic}),
                         'cellname': key})
    temp['dist'] = temp.index
    temp.set_index('cellname', inplace=True)
    neckregion = neckregion.append(temp, ignore_index=False)


cellall['budvol'] = dfmb.bud
cellall['momvol'] = dfmb.mom
cellall = cellall.reset_index()
cellall['type'] = cellall['index'].apply(lambda x: x[:3])
cellposbud = cellposbud.reset_index()
cellposmom = cellposmom.reset_index()
cellposbud = pd.concat([cellposbud, cellall.ix[:, ['type', 'neck']]], axis=1)
cellposmom = pd.concat([cellposmom, cellall.ix[:, ['type', 'neck']]], axis=1)

cellall['frac'] = cellall.ix[:, 'bud'] / cellall.ix[:, 'mom']
Q = cellall.groupby('type').quantile(.90)  # 90th percentile of each cols
#  q90 = 90th percentile bud volume of each media type
cellall['q90'] = cellall.type.apply(lambda x: Q.ix[x].budvol)
gt90 = cellall[cellall['budvol'] > cellall['q90']]
meangt90 = gt90.groupby('type').budvol.mean()
cellall['mean90'] = cellall.type.apply(lambda x: meangt90.ix[x])
#  budvolratio is based on the largest 10% cells
cellall['budvolratio'] = cellall.budvol / cellall.q90
cellposbud['budvol'] = cellall['budvol']
cellposmom['momvol'] = cellall['momvol']

cellposbud['binvol'] = vf.bincell(cellposbud, 'budvol', binsvolbud)
cellposmom['binvol'] = vf.bincell(cellposmom, 'momvol', binsvolmom)

# =============================================================================
# cells binned by budding progression
# =============================================================================
BIG = pd.melt(cellall,
              id_vars=['type'],
              value_vars=['frac'])
groups = BIG.groupby('type')

BIG2 = pd.melt(cellall,
               id_vars=['type'],
               value_vars=['mom', 'bud'])

cellall['bin_budprog'] = vf.bincell(cellall, 'budvolratio', binsbudprog)


cellall['binbudvol'] = cellposbud['binvol']

# reject super large cells
rejectlist = cellposmom.ix[(np.asarray(cellposmom.momvol) > 60) &
                           (cellposmom.type != 'YPD'), 'index']
cellall = cellall.ix[~ cellall.ix[:, 'index'].isin(rejectlist)]
cellposmom = cellposmom.ix[~cellposmom.ix[:, 'index'].isin(rejectlist)]
cellposbud = cellposbud.ix[~cellposbud.ix[:, 'index'].isin(rejectlist)]

# =============================================================================
# Progression of Δψ as move along the bud axis
# =============================================================================
sns.set_style('whitegrid')
bigbinsmom = pd.melt(cellposmom,
                     id_vars=['type', 'binvol'],
                     var_name='mom axis position',
                     value_name=r'$\Delta\Psi$ scaled gradient',
                     value_vars=binsaxis.tolist())
bigbinsmom = bigbinsmom.dropna()
bigbinsbud = pd.melt(cellposbud,
                     id_vars=['type', 'binvol'],
                     var_name='bud axis position',
                     value_name=r'$\Delta\Psi$ scaled gradient',
                     value_vars=binsaxis.tolist())
bigbinsbud = bigbinsbud.dropna()
#
with sns.plotting_context('talk', font_scale=1.1):
    h = sns.FacetGrid(bigbinsmom,
                      col="type",
                      hue='type',
                      col_wrap=2)
    h = h.map(sns.pointplot,
              'mom axis position',
              r'$\Delta\Psi$ scaled gradient')
    h.set_xlabels('mom cell axis position')
    h.set(ylim=(0, 1.))

    m = sns.FacetGrid(bigbinsbud,
                      row="type",
                      col="binvol",
                      hue='type',
                      col_order=binsvolbud[1:])
    m = m.map(sns.pointplot,
              'bud axis position',
              r'$\Delta\Psi$ scaled gradient')
    m.set_xticklabels(fontsize=14)
    m.set_xlabels("bud axis position")
    m.set(ylim=(0, 1.))

# =============================================================================
# frac Δψ as function of budratio
# =============================================================================
with sns.plotting_context('talk'):
    _, ax2 = plt.subplots(1, 1)
    h = sns.pointplot(x='bin_budprog',
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
    p = p.map(sns.pointplot, 'bin_budprog', 'frac')

    q = sns.FacetGrid(cellall.dropna(),
                      col="binbudvol",
                      hue='type',
                      row='type')
    q = q.map(sns.pointplot,
              'bin_budprog',
              'frac')
    q.set_xlabels('bud progression')


# =============================================================================
#     Δψ at the bud neck region
# =============================================================================
with sns.plotting_context('talk'):
    _, ax3a = plt.subplots(1, 1)
    A = pd.melt(neckregion,
                id_vars=['dist'],
                value_vars=['bud', 'mom'])
A.dropna(inplace=True)
with sns.plotting_context('talk', font_scale=1.4):
    sns.barplot(x='dist', y='value', hue='variable', data=A, ax=ax3a)
#  ============================================================================
#  frac Δψ violinplots by media
#  ============================================================================
#with sns.plotting_context('talk'):
#    _, ax1 = plt.subplots(1, 1)
#    g = sns.violinplot(x='type',
#                       y='value',
#                       hue='type',
#                       data=BIG,
#                       ax=ax1)
#    g.set_ylim(0, 3)
#    g.get_legend().set_visible(False)
#    g = sns.stripplot(x='type',
#                      split=True,
#                      y='value',
#                      hue='type',
#                      data=BIG,
#                      jitter=.15,
#                      ax=ax1)
##    g.set_ylim(0.4, 2.0)
#    g.get_legend().set_visible(False)
# =============================================================================
# violinplot mom vs bud Δψ scaled
# =============================================================================
with sns.plotting_context('talk', font_scale=1.4):
    _, ax3 = plt.subplots(1, 1)
    h = sns.violinplot(x='type',
                       y='value',
                       hue='variable',
                       data=BIG2,
                       ax=ax3)
    sns.stripplot(x='type',
                  y='value',
                  hue='variable',
                  jitter=.15,
                  size=4,
                  data=BIG2,
                  ax=ax3)
    h.set_ylim(0, 1.)
    h.get_legend().set_visible(False)
#
# =============================================================================
# frac Δψ as function of budvol
# =============================================================================
# with sns.plotting_context('talk', font_scale=1.4):
##    _, ax10 = plt.subplots(1, 1)
###    g = sns.FacetGrid(cellall.dropna(), col="type")
###    g = g.map(sns.regplot, "budvol", "frac")
##    datacell = cellall[cellall.bud <= bins2[-1]]
# h = sns.pointplot(x='binbudratio',
# y='bud',
# hue='type',
# data=datacell.dropna(),
# ax=ax10)
# h.get_legend().set_visible(False)
##
# for i in ['YPD', 'YPE', 'YPL', 'YPR']:
##        data = cellall[(cellall.type == i) & (cellall.frac < 2)]
# slope, _, r, p, _ = sp.linregress(data['budvol'],
# data['frac'])
# print 'slope= %6.4f r=%6.4f p=%6.4f' % (slope, r, p)
#
# =============================================================================
# Dy as budneckregion and budratio
# =============================================================================
# with sns.plotting_context('talk', font_scale=1.4):
##    _, ax2 = plt.subplots(1, 1)
# h = sns.pointplot(x='bin_budprog',
# y='DYneck',
# hue='type',
# data=cellall.dropna(),
# ax=ax2)
# h.get_legend().set_visible(False)
#
#
# with sns.plotting_context('talk', font_scale=1.4):
##    _, ax1 = plt.subplots(1, 1)
# h = sns.pointplot(x='posx',
# y='DY',
# ci=None,
# markers='o',
# join=False,
# hue='type',
# data=cell,
# size=1,
# ax=ax1)
# h.get_legend().set_visible(False)
##    h.set_xticks(np.linspace(cell.pos.min(), cell.pos.max(),11))
##    h.set_xticklabels(np.arange(0, 1.1 ,.1))
#
# ==============================================================================
# budratio
# ==============================================================================
##c2 = cellall.drop(cellall.index[[5, 15, 63, 46]])
# slope, _, r, p, std_err = sp.linregress(c2.ix[:, 'budratio'],
# c2.ix[:, 'neck'])
# with sns.plotting_context('talk', font_scale=1.4):
##    _, ax5 = plt.subplots(1, 1)
# h = sns.regplot(x='budratio',
# y='neck',
# data=c2[c2.neck>0.505],
# ax=ax5)
# h.get_legend().set_visible(False)
