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
import pandas as pd
import seaborn as sns
from mombud.functions import vtk_mbfuncs as vf
import wrappers as wr
# pylint: disable=C0103

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.close('all')

#==============================================================================
# Data directory prep
#==============================================================================
datadir = op.join(os.getcwd(), 'mutants', 'transformedData', 'filtered')
datadir_old = op.join(os.getcwd(), 'data', 'transformedData')

with open(op.join(datadir, 'mombudtrans_new.pkl'), 'rb') as inpt:
    dfmb = pickle.load(inpt)  # has columns base, neck, tip, media, bud, mom

with open(op.join(datadir_old, 'mombudtrans.pkl'), 'rb') as inpt:
    dfmb_o = pickle.load(inpt)  # has columns base, neck, tip, media, bud, mom

rejectfold = op.join(datadir, os.pardir, 'reject')
reject = wr.swalk(rejectfold, '*png', stop=-4)

try:
    filext = "*vtk"
    vtkF = wr.ddwalk(datadir, filext, stop=-4)
except LookupError:
    sys.exit("error filetypes %s not found in %s" % (filext, datadir))

try:
    filext = "*vtk"
    vtkF_old = wr.swalk(datadir_old, filext, stop=-4)
except LookupError:
    sys.exit("error filetypes %s not found in %s" % (filext, datadir))


filekeys_old = {item: vtkF_old[item] for item
                in sorted(vtkF_old.keys()) if item.split('_')[0] != 'YPD'
                and item not in reject}

filekeys = {item: vtkF[media][item] for media
            in sorted(vtkF.keys()) for item
            in sorted(vtkF[media].keys())}

filekeys_f = {f: filekeys[f] for f in filekeys if f not in reject}
filekeys_f.update(filekeys_old)  # add YPE to dict
# dataframe of neck, mom and bud tip positions, bud and mom volumes
dfmb = dfmb.append(dfmb_o[dfmb_o.media != 'YPD'])

#=============================================================================
# Dataframe and binning setup
#=============================================================================
# bins for binning the bud progression ratio
binsbudprog = np.r_[np.arange(0, 1.1, .1), 2]
binsaxis = np.linspace(0., 1., 6)  # pos. along mom/bud cell
binsaxisbig = np.linspace(0, 1., 21)  # position along whole cell
binsvolbud = np.linspace(0, 40, 5)  # vol binning for bud
binsvolmom = np.array([0, 30, 40, 80.])  # vol binning for mom

# dataframe for budding progression and budratio, size distribution , fra Δψ
cellall = pd.DataFrame(columns=['mom', 'bud'])

# dataframe for Δψ distributino along mom/bud cell axis
cellposmom = pd.DataFrame()
cellposbud = pd.DataFrame()
cellposmom_fp = pd.DataFrame()
# dataframe for average Δψ around the neck region
neckregion = pd.DataFrame()

#=============================================================================
# compute Δψ distrbution along cellaxis for each ind cell and collect/append
# to the dataframes
#=============================================================================
for key in sorted(filekeys_f)[:]:

    # returns Dataframe of pos along x-axis for inidivual mom and bud cell
    cell = vf.cellpos(filekeys_f[key], dfmb)

    # bin the dataframe according to binxaxis
    cell['ind_cell_pos'] = vf.bincell(cell, 'indcell_xaxis', binsaxis)

    # Series of average Δψ (scaled to cell minmax values, DY_minmax)
    Xmom = cell.ix[cell['type'] == 'mom'].groupby('ind_cell_pos').DY.mean()
    Xbud = cell.ix[cell['type'] == 'bud'].groupby('ind_cell_pos').DY.mean()
    fp = cell[cell['type'] == 'bud'][:1]

    # pos along x-axis for the whole cell
    cell['whole_cell_pos'] = vf.bincell(cell, 'wholecell_xaxis', binsaxisbig)
    Xcell = cell.groupby('whole_cell_pos').DY.mean()
    Xcell = Xcell[Xcell < cell.neckpos_cellaxis.max()].reset_index()
    if len(fp.DY.values):
        Xcell = Xcell.append(
            pd.Series({'DY': fp.DY.values[0]}), ignore_index=True)
    else:
        Xcell = Xcell.append(pd.Series({'DY': 0}), ignore_index=True)
    g = Xcell.whole_cell_pos.cat.add_categories('fp')
    g.fillna('fp', inplace=True)
    Xcell['cellaxis_mom_budfp'] = g

    # Series of Δψ scaled to min-max of the MOM-BUD cell AXIS
    scaled_dy_wholecell = cell.groupby('whole_cell_pos').DY.mean()
    dy_wholecell_mean = cell.DY.mean()
    dy_wholecell_abs = cell.DY_abs.mean()
    dy_budmom_abs = cell.groupby('type').DY_abs.mean()
    dy_budmom_abs.rename({'bud': 'bud_abs_dy',
                          'mom': 'mom_abs_dy'}, inplace=True)
    dy_wholecell = cell.ix[:, ['DY', 'DY_abs']].mean()
    dy_wholecell.rename({'DY': 'whole_cell_mean',
                         'DY_abs': 'whole_cell_abs'}, inplace=True)

#    Xmom = vf.scaleminmax(Xmom, scaled_dy_wholecell)
#    Xbud = vf.scaleminmax(Xbud, scaled_dy_wholecell)

    Xmom = Xmom / dy_wholecell_mean
    Xbud = Xbud / dy_wholecell_mean

    Xbud.name = key
    Xmom.name = key
    DYseries = cell.groupby('type').median().DY
    DYseries = DYseries.append([dy_wholecell, dy_budmom_abs])
    DYseries.name = key
    cellall = cellall.append(DYseries)
    cellposbud = cellposbud.append(Xbud)
    cellposmom = cellposmom.append(Xmom)
    cellposmom_fp = cellposmom_fp.append(Xcell)  # momcell + first bud point
    # temp dict of mean Δψ at +- range of dist from budneck
    tempdic = {dist: vf.neckDY(cell, cell.neckpos, dist)
               for dist in [.15, .3, .5]}
    temp = pd.DataFrame({'bud': pd.Series({k: tempdic[k][0] for k in tempdic}),
                         'mom': pd.Series({k: tempdic[k][1] for k in tempdic}),
                         'cellname': key})
    temp['dist'] = temp.index
    temp.set_index('cellname', inplace=True)
    neckregion = neckregion.append(temp, ignore_index=False)

# cleanup and add. labels for dataframes
cellall['budvol'] = dfmb.bud
cellall['momvol'] = dfmb.mom
cellall = cellall.reset_index()
cellall['type'] = cellall['index'].apply(lambda x: x.split('_')[0])
cellposbud = cellposbud.reset_index()
cellposmom = cellposmom.reset_index()
cellposmom_fp = cellposmom_fp.reset_index()
cellposbud = pd.concat([cellposbud, cellall.ix[:, ['type']]], axis=1)
cellposmom = pd.concat([cellposmom, cellall.ix[:, ['type']]], axis=1)
cellposmom_fp = pd.concat([cellposmom_fp, cellall.ix[:, ['type']]], axis=1)

cellall['frac'] = cellall.ix[:, 'bud'] / cellall.ix[:, 'mom']
Q = cellall.groupby('type').quantile(.90)  # 90th percentile of each cols

#  q90 = 90th percentile bud volume of each media type
cellall['q90'] = cellall.type.apply(lambda x: Q.ix[x].budvol)
gt90 = cellall[cellall['budvol'] > cellall['q90']]
meangt90 = gt90.groupby('type').budvol.mean()
cellall['mean90'] = cellall.type.apply(lambda x: meangt90.ix[x])

#  budvolratio is based on the largest 10% cells
cellall['budvolratio'] = cellall.budvol / cellall.q90


cellall['date'] = cellall['index'].apply(lambda x: x.split('_')[1])
cellall['date'] = cellall.date.apply(lambda x: stripc(x))
YPE = cellall[(cellall.type == 'YPE') | (cellall.type == 'WT')].copy()
YPL = cellall[(cellall.type == 'YPL')].copy()
YPE = YPE.reset_index()

cellall = cellall.ix[~(((cellall.type == 'YPE') & (cellall.date == '052315')) |
                       ((cellall.type == 'WT') & (cellall.date == '032716')))]


cellposbud['budvol'] = cellall['budvol']
cellposmom['momvol'] = cellall['momvol']

cellposbud['binvol'] = vf.bincell(cellposbud, 'budvol', binsvolbud)
cellposmom['binvol'] = vf.bincell(cellposmom, 'momvol', binsvolmom)
#
# =============================================================================
# cells binned by budding progression
# =============================================================================
cellall['bin_budprog'] = vf.bincell(cellall, 'budvolratio', binsbudprog)
cellall['binbudvol'] = cellposbud['binvol']

# reject super large cells
rejectlist = cellposmom.ix[(np.asarray(cellposmom.momvol) > 100) &
                           (cellposmom.type != 'YPD'), 'index']
cellall = cellall.ix[~ cellall.ix[:, 'index'].isin(rejectlist)]
cellposmom = cellposmom.ix[~cellposmom.ix[:, 'index'].isin(rejectlist)]
cellposbud = cellposbud.ix[~cellposbud.ix[:, 'index'].isin(rejectlist)]

# =============================================================================
# Distribution of bud and mom volumes
# =============================================================================
budvol = cellall.ix[:, ['budvol', 'type', 'N']]
momvol = cellall.ix[:, ['momvol', 'type', 'N']]
budvol['N'] = budvol.groupby("type").transform('count')
momvol['N'] = budvol.groupby("type").transform('count')


def stripc(s):
    if s.startswith('c'):
        snew = s.replace(s[0], '0')
    else:
        snew = s
    return snew

Nype = YPE.groupby('date').index.count().to_dict()
Nbud = cellall.groupby(['type', 'binbudvol']).size()
N = cellall.groupby('type').index.count()
N = N.to_dict()  # dict to hold counts of each type
col_ord = ['MFB1', 'NUM1', 'YPT11', 'WT', 'YPE', 'YPL', 'YPR']
hue_ord = ['mom_abs_dy', 'bud_abs_dy', 'whole_cell_abs']


def label_n(handle, labeldic):
    """
    modifies title on facetgrid to include labeldic

    Parameters
    ----------

    handle : FacetGrid ref
        handle to FacetGrid obj

    labeldic : dict
        dictionary of text labels to be added to handle's title

    """
    if hasattr(handle.axes, 'flat'):
        for ax in handle.axes.flat:
            oldtitle = ax.get_title()
            print oldtitle
            if oldtitle.find('|') > -1:
                media = oldtitle.split('|')[0].split('=')[1].strip()
                budv = oldtitle.split('=')[-1].strip()
                newtitle = '{}, N = {}'.format(media,
                    labeldic.xs(media).get([float(budv)])[0])

                print "{} ,{}, {}".format(media, budv, newtitle)
                ax.set_title(newtitle)
            else:
                oldtitle = oldtitle.split('=')[1].strip()
                ax.set_title('%s, N=%d' % (oldtitle, labeldic[oldtitle]))
    else:

        labels = [xl.get_text().strip()
                  for xl in handle.axes.get_xticklabels()]
        new_labels = [
            '%s\n N=%d' %
            (old_lab, labeldic[old_lab]) for old_lab in labels]
        handle.axes.set_xticklabels(new_labels)

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


def plot(save=False):
    sns.set_style('whitegrid')
    with sns.plotting_context('talk', font_scale=1.1):
        g = sns.FacetGrid(budvol,
                          col="type",
                          col_wrap=4,
                          hue="type",
                          col_order=col_ord)
        g = (g.map(sns.distplot, "budvol")).set(xlim=(0.))
        label_n(g, N)
        if save:
            g.savefig(op.join(datadir, 'budsize_dist.png'))

        h = sns.FacetGrid(momvol,
                          col="type",
                          col_wrap=4,
                          hue="type",
                          col_order=col_ord)
        h = (h.map(sns.distplot, "momvol")).set(xlim=(0.))
        label_n(h, N)
        if save:
            h.savefig(op.join(datadir, 'momsize_dist.png'))
        # =================================================================
        # Progression of Δψ as move along the bud axis
        # =================================================================

        h = sns.FacetGrid(bigbinsmom,
                          col="type",
                          hue='type',
                          col_wrap=4,
                          sharex=True,
                          col_order=col_ord)
        h = h.map(sns.pointplot,
                  'mom axis position',
                  r'$\Delta\Psi$ scaled gradient').set(ylim=(0.7, 1.5))
        label_n(h, N)
        if save:
            h.savefig(op.join(datadir, 'mom_cell_dy.png'))

        m1 = sns.FacetGrid(bigbinsbud,
                           col="type",
                           hue='type',
                           col_wrap=4,
                           col_order=col_ord)

        m1 = m1.map(sns.pointplot,
                    'bud axis position',
                    r'$\Delta\Psi$ scaled gradient').set(ylim=(0.7, 1.5))
        label_n(m1, N)
        if save:
            m1.savefig(op.join(datadir, 'bud_cell_dy.png'))

    # with facetting by budvol
    with sns.plotting_context('talk', font_scale=.9):
        m0 = sns.FacetGrid(bigbinsbud,
                           row="type",
                           col="binvol",
                           hue='type',
                           row_order=col_ord,
                           col_order=binsvolbud[1:])

        m0 = m0.map(sns.pointplot,
                    'bud axis position', r'$\Delta\Psi$ scaled gradient'
                    ).set(yticks=np.arange(0.5, 1.9, 0.25),
                          ylim=(0.65, 2.))
        label_n(m0, Nbud)

        if save:
            m0.savefig(op.join(datadir, 'bud_cell_dy_facetted.png'))
    # =========================================================================
    # frac Δψ as function of budratio
    # =========================================================================
    with sns.plotting_context('talk'):
        _, ax2 = plt.subplots(1, 1)
        h = (sns.pointplot(x='bin_budprog',
                           y='frac',
                           hue='type',
                           data=cellall.dropna(),
                           ax=ax2))
        h.set(ylim=(0, 3),
              title=u"Δψ vs bud progression\n ",
              xlabel="bud progression",
              ylabel=u"Δψ bud/Δψ mom")
        leg = h.get_legend()
        plt.setp(leg, bbox_to_anchor=(0.85, 0.7, .3, .3))
        if save:
            plt.savefig(op.join(datadir, "DY vs bud progression.png"))

        p = sns.FacetGrid(cellall.dropna(),
                          col="type",
                          hue='type',
                          col_wrap=4,
                          col_order=col_ord)
        p = p.map(sns.pointplot, 'bin_budprog', 'frac')
        if save:
            p.savefig(op.join(datadir, "DY_bud_prog_facetted.png"))

    # ========================================================================
    #     Δψ at the bud neck region
    # ========================================================================
    with sns.plotting_context('talk'):
        A = pd.melt(neckregion,
                    id_vars=['dist'],
                    value_vars=['bud', 'mom'])
    A.dropna(inplace=True)
    with sns.plotting_context('talk', font_scale=1.4):
        _, ax1 = plt.subplots(1, 1)
        q1 = sns.barplot(x='dist',
                         y='value',
                         hue='variable',
                         data=A,
                         ax=ax1)
        leg = q1.get_legend()
        plt.setp(leg, bbox_to_anchor=(0.85, 0.7, .3, .3))
        if save:
            plt.savefig(op.join(datadir, "neckregionDY.png"))

    #  ========================================================================
    #  frac Δψ violinplots by media
    #  ========================================================================
    BIG = pd.melt(cellall,
                  id_vars=['type'],
                  value_vars=['frac'])

    with sns.plotting_context('talk'):
        _, ax4 = plt.subplots(1, 1)
        j = sns.violinplot(x='type',
                           y='value',
                           hue='type',
                           data=BIG.dropna(),
                           order=col_ord,
                           inner=None,
                           ax=ax4)
        j.set_ylim(0, 2.5)
        j.get_legend().set_visible(False)

        k = sns.boxplot(x='type',
                        y='value',
                        hue='type',
                        data=BIG.dropna(),
                        order=col_ord,
                        showmeans=True,
                        showbox=False,
                        showcaps=False,
                        showfliers=False,
                        medianprops={'linewidth': 0},
                        whiskerprops={'linewidth': 0},
                        meanprops={'marker': '_',
                                   'c': 'w',
                                   'ms': 5,
                                   'markeredgewidth': 2},
                        ax=ax4)
        k.get_legend().set_visible(False)
        label_n(j, N)
        if save:
            plt.savefig(op.join(datadir, "violin_fracDY.png"))

    # =====================================================================
    # violinplot mom vs bud Δψ scaled
    # =====================================================================
    BIG2 = pd.melt(cellall,
                   id_vars=['type'],
                   value_vars=['mom', 'bud'])

    with sns.plotting_context('talk', font_scale=1.):
        _, ax3 = plt.subplots(1, 1)
        h = sns.violinplot(x='type',
                           y='value',
                           hue='variable',
                           order=col_ord,
                           data=BIG2.dropna(),
                           ax=ax3)
        h.set_ylim(0, 1.)
        h.get_legend().set_visible(False)
        label_n(h, N)
        if save:
            plt.savefig(op.join(datadir, "Violin Mom_Bud_DY.png"))

    with sns.plotting_context('talk', font_scale=1.):
        BIG4 = pd.melt(cellall,
                       id_vars=['type'],
                       value_vars=['whole_cell_abs'])

        g = sns.FacetGrid(BIG4,
                          col="type",
                          col_wrap=4,
                          hue="type",
                          col_order=col_ord,
                          )
        g = (g.map(sns.distplot, "value")).set(xlim=(0.))
        label_n(g, N)

    # ========================================================================
    # violinplot mom vs bud Δψ raw by date for YPE
    # =====================================================================
    with sns.plotting_context('talk', font_scale=1.):
        BIG5 = pd.melt(YPE,
                       id_vars=['date'],
                       value_vars=['whole_cell_abs',
                                   'bud_abs_dy',
                                   'mom_abs_dy'])

        _, ax7 = plt.subplots(1, 1)
        g = sns.violinplot(x='date',
                           y='value',
                           hue='variable',
                           hue_order=hue_ord,
                           data=BIG5,
                           ax=ax7)
        leg = g.get_legend()
        plt.setp(leg,
                 bbox_to_anchor=(.75, .85, .1, .2))
        g.set_ylim(0, 4000)
        label_n(g, Nype)

        if save:
            plt.savefig(op.join(datadir, "Violin-DY_raw_ype_date.png"))

    # ========================================================================
    # violinplot raw by date for all types
    # =====================================================================
    with sns.plotting_context('talk', font_scale=1.):
        BIG6 = pd.melt(cellall,
                       id_vars=['type'],
                       value_vars=['whole_cell_abs',
                                   'bud_abs_dy',
                                   'mom_abs_dy'])

        _, ax8 = plt.subplots(1, 1)
        g = sns.violinplot(x='type',
                           y='value',
                           bw='scott',
                           hue='variable',
                           order=col_ord,
                           hue_order=hue_ord,
                           data=BIG6,
                           ax=ax8)
        leg = g.get_legend()
        plt.setp(leg,
                 bbox_to_anchor=(.75, .85, .1, .2))
        g.set_ylim(0, 2000)
        label_n(g, N)

        if save:
            plt.savefig(op.join(datadir, "Violin-DY_raw_ype_date.png"))

#==============================================================================
if __name__ == '__main__':
    pass
    plot(save=True)
