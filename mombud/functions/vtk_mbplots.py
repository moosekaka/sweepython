# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 16:10:16 2016
Module for plots of analysis of mother bud function in budding yeast
"""
import os.path as op
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from mombud.functions import vtk_mbfuncs as vf
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.close('all')

COL_ODR = ['MFB1', 'NUM1', 'YPT11', 'WT', 'YPE', 'YPL', 'YPR']
HUE_ODR = ['DY_abs_mean_mom', 'DY_abs_mean_bud', 'whole_cell_abs']

# pylint: disable=C0103
def getrval(df, x, y, labeldic):
    """
    return a subset DataFrame and R^2 values for columns x, y in original df
    """
    df = df.ix[:, [x, y, 'date']].reset_index(drop=True)
    df.rename(columns=labeldic, inplace=True)
    pear = df.groupby('date').corr().xs(labeldic[y],
                                        level=1).raw.to_dict()
    r_sqr = {key: value**2 for key, value in pear.iteritems()}

    return df, r_sqr


def label_n(handle, labeldic, Rsqr=False):
    """
    modifies title on facetgrid to include labeldic

    Parameters
    ----------

    handle : FacetGrid ref
        handle to FacetGrid obj

    labeldic : dict
        dictionary of text labels to be added to handle's title
    Rsqr : Bool
        if True, labels R^2 instead of N counts
    """

    if hasattr(handle.axes, 'flat'):
        for ax in handle.axes.flat:
            oldtitle = ax.get_title()

            if oldtitle.find('|') > -1:
                media = oldtitle.split('|')[0].split('=')[1].strip()
                budv = oldtitle.split('=')[-1].strip()
                newtitle = '{}, N = {}'.format(
                    media, labeldic.xs(media).get([float(budv)])[0])
                ax.set_title(newtitle)
            else:
                oldtitle = oldtitle.split('=')[1].strip()
                if not Rsqr:
                    ax.set_title('{}, N={}'
                                 .format(oldtitle,
                                         labeldic[oldtitle]))
                else:
                    ax.set_title('{}, R^2={:5.3f}'
                                 .format(oldtitle,
                                         labeldic[oldtitle]))
    else:

        labels = [xl.get_text().strip()
                  for xl in handle.axes.get_xticklabels()]
        new_labels = ['{}\n N={}'.format(
            old_lab, labeldic[old_lab]) for old_lab in labels]
        handle.axes.set_xticklabels(new_labels)


def plotDims(**kwargs):
    """
    Diameters of cells
    """
    df = kwargs.get('data')
    df2 = pd.DataFrame(x for x in df['cell_diameter'])
    df2.rename(columns={'bud': 'bud_diameter',
                        'mom': 'mom_diameter'},
               inplace=True)
    alldims = pd.concat([df, df2], axis=1)
    melt = pd.melt(alldims,
                   id_vars='type',
                   value_vars=['bud_diameter', 'mom_diameter'])

    with sns.plotting_context('talk'):
        g = sns.FacetGrid(melt,
                          col='type',
                          col_wrap=4,
                          hue="variable",
                          col_order=COL_ODR,
                          size=3,
                          aspect=1.5)
        g = (g.map(sns.stripplot,
                   "value", jitter=0.1)).set(xlim=(0.),
                                             xlabel='diameter/microns')


def plotSizeDist(**kwargs):
    """
    Distribution of bud and mom volumes
    """
    df = kwargs.get('data')
    datadir = kwargs.get('savefolder')
    N = kwargs.get('counts')
    save = kwargs.get('save', False)
    COL_ODR = kwargs.get('COL_ODR')

    budvol = df.ix[:, ['budvol', 'type']]
    momvol = df.ix[:, ['momvol', 'type']]
    budvol['N'] = budvol.groupby("type").transform('count')
    momvol['N'] = momvol.groupby("type").transform('count')
    sns.set_style('whitegrid')
    with sns.plotting_context('talk', font_scale=1.1):
        g = sns.FacetGrid(budvol,
                          col="type",
                          col_wrap=4,
                          hue="type",
                          col_order=COL_ODR)
        g = (g.map(sns.distplot, "budvol")).set(xlim=(0.))
        label_n(g, N)
        if save:
            g.savefig(op.join(datadir, 'budsize_dist.png'))

        h = sns.FacetGrid(momvol,
                          col="type",
                          col_wrap=4,
                          hue="type",
                          col_order=COL_ODR)
        h = (h.map(sns.distplot, "momvol")).set(xlim=(0.))
        label_n(h, N)
        if save:
            h.savefig(op.join(datadir, 'momsize_dist.png'))


def plotDyAxisDist(dfmom, dfbud, **kwargs):
    """
    Progression of Δψ as move along the bud axis
    """
    binsaxis = kwargs.get('mbax')
    binsvolbud = kwargs.get('binsvolbud')
    Nbud = kwargs.get('counts_buds')
    datadir = kwargs.get('savefolder')
    N = kwargs.get('counts')
    save = kwargs.get('save', False)
    COL_ODR = kwargs.get('COL_ODR')

    bigbinsmom = pd.melt(dfmom,
                         id_vars=['type', 'binvol'],
                         var_name='mom axis position',
                         value_name=r'$\Delta\Psi$ scaled gradient',
                         value_vars=binsaxis.tolist())
    bigbinsmom = bigbinsmom.dropna()
    bigbinsbud = pd.melt(dfbud,
                         id_vars=['type', 'binvol'],
                         var_name='bud axis position',
                         value_name=r'$\Delta\Psi$ scaled gradient',
                         value_vars=binsaxis.tolist())
    bigbinsbud = bigbinsbud.dropna()
    sns.set_style('whitegrid')
    with sns.plotting_context('talk', font_scale=1.):
        h = sns.FacetGrid(bigbinsmom,
                          col="type",
                          hue='type',
                          col_wrap=4,
                          sharex=True,
                          col_order=COL_ODR)
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
                           col_order=COL_ODR)

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
                           row_order=COL_ODR,
                           col_order=binsvolbud[1:])

        m0 = m0.map(sns.pointplot,
                    'bud axis position',
                    r'$\Delta\Psi$ scaled gradient').set(
                        yticks=np.arange(0.5, 1.9, 0.25), ylim=(0.65, 2.))
        label_n(m0, Nbud)

        if save:
            m0.savefig(op.join(datadir, 'bud_cell_dy_facetted.png'))


def plotBudProgr(**kwargs):
    """
    frac Δψ as function of budratio
    """
    df = kwargs.get('data')
    datadir = kwargs.get('savefolder')
    save = kwargs.get('save', False)
    COL_ODR = kwargs.get('COL_ODR')
    sns.set_style('whitegrid')

    with sns.plotting_context('talk'):
        _, ax2 = plt.subplots(1, 1)
        h = (sns.pointplot(x='bin_budprog',
                           y='frac',
                           hue='type',
                           data=df,
                           ax=ax2))
        h.set(ylim=(0, 3),
              title=u"Δψ vs bud progression\n ",
              xlabel="bud progression",
              ylabel=u"Δψ bud/Δψ mom")
        leg = h.get_legend()
        plt.setp(leg, bbox_to_anchor=(0.85, 0.7, .3, .3))
        if save:
            plt.savefig(op.join(datadir, "DY vs bud progression.png"))

        p = sns.FacetGrid(df,
                          col="type",
                          hue='type',
                          col_wrap=4,
                          col_order=COL_ODR)
        p = p.map(sns.pointplot, 'bin_budprog', 'frac')
        if save:
            p.savefig(op.join(datadir, "DY_bud_prog_facetted.png"))


def plotViolins(**kwargs):
    """
    Violinplots for frac Δψ, mom vs bud scaled and Δψ abs distr
    """
    df = kwargs.get('data')
    datadir = kwargs.get('savefolder')
    save = kwargs.get('save', False)
    N = kwargs.get('counts')
    COL_ODR = kwargs.get('COL_ODR')

    BIG = pd.melt(df,
                  id_vars=['type'],
                  value_vars=['frac'])

    BIG2 = pd.melt(df,
                   id_vars=['type'],
                   value_vars=['DY_median_mom', 'DY_median_bud'])

    sns.set_style('whitegrid')
    with sns.plotting_context('talk'):
        _, ax4 = plt.subplots(1, 1)
        j = sns.violinplot(x='type',
                           y='value',
                           hue='type',
                           data=BIG.dropna(),
                           order=COL_ODR,
                           inner=None,
                           ax=ax4)
        j.set_ylim(0, 2.5)
        j.get_legend().set_visible(False)

        k = sns.boxplot(x='type',
                        y='value',
                        hue='type',
                        data=BIG.dropna(),
                        order=COL_ODR,
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

        _, ax3 = plt.subplots(1, 1)
        h = sns.violinplot(x='type',
                           y='value',
                           hue='variable',
                           order=COL_ODR,
                           data=BIG2.dropna(),
                           ax=ax3)
        h.set_ylim(0, 1.)
        h.get_legend().set_visible(False)
        label_n(h, N)
        if save:
            plt.savefig(op.join(datadir, "Violin Mom_Bud_DY.png"))

        BIG4 = pd.melt(df,
                       id_vars=['type'],
                       value_vars=['whole_cell_abs'])

        g = sns.FacetGrid(BIG4,
                          col="type",
                          col_wrap=4,
                          hue="type",
                          col_order=COL_ODR,
                          size=3,
                          aspect=1.5)
        g = (g.map(sns.distplot, "value")).set(xlim=(0.))
        label_n(g, N)

        if save:
            plt.savefig(op.join(datadir, "DY_abs_dist.png"))


def plotGFP(**kwargs):
    """
    Violinplots for GFP uptake by date and carbon type
    """
    df = kwargs.get('data')
#    dfype = kwargs.get('data_ype')
    N = kwargs.get('counts')
    Ndate = kwargs.get('counts_date')
#    Nype = kwargs.get('counts_ype')
    datadir = kwargs.get('savefolder')
    save = kwargs.get('save', False)
    COL_ODR = kwargs.get('COL_ODR')
    HUE_ODR = kwargs.get('HUE_ODR')

    with sns.plotting_context('talk', font_scale=1.):
        BIG5 = pd.melt(df,
                       id_vars=['date'],
                       value_vars=['whole_cell_abs',
                                   'DY_abs_mean_bud',
                                   'DY_abs_mean_mom'])

        _, ax7 = plt.subplots(1, 1)
        g = sns.violinplot(x='date',
                           y='value',
                           hue='variable',
                           hue_order=HUE_ODR,
                           data=BIG5,
                           ax=ax7)
        leg = g.get_legend()
        plt.setp(leg,
                 bbox_to_anchor=(.75, .85, .1, .2))
        g.set_ylim(0, 4000)
        label_n(g, Ndate)

        if save:
            plt.savefig(op.join(datadir, "Violin-GFP_by_date.png"))

        BIG6 = pd.melt(df,
                       id_vars=['type'],
                       value_vars=['whole_cell_abs',
                                   'DY_abs_mean_bud',
                                   'DY_abs_mean_mom'])

        _, ax8 = plt.subplots(1, 1)
        g = sns.violinplot(x='type',
                           y='value',
                           bw='scott',
                           hue='variable',
                           order=COL_ODR,
                           hue_order=HUE_ODR,
                           data=BIG6,
                           ax=ax8)
        leg = g.get_legend()
        plt.setp(leg,
                 bbox_to_anchor=(.75, .85, .1, .2))
        g.set_ylim(0, 4000)
        label_n(g, N)

        if save:
            plt.savefig(op.join(datadir, "Violin-GFP_by_type.png"))


def plotRegr(**kwargs):
    """
    plot regression coeff of Δψ raw vs scaled
    """
    df = kwargs.get('data')
    datadir = kwargs.get('savefolder')
    save = kwargs.get('save', False)

    labeldic = {'whole_cell_abs': 'raw',
                'whole_cell_mean': 'scaled cell mean',
                'DY_median_mom': 'scaled mom mean',
                'DY_median_bud': 'scaled bud mean'}

    with sns.plotting_context('talk', font_scale=1.):
        for p in ['whole_cell_mean', 'DY_median_mom', 'DY_median_bud']:
            a, r2 = getrval(df, 'whole_cell_abs', p, labeldic)
            x, y, _ = a.columns
            sns.lmplot(x, y,
                       fit_reg=False,
                       legend_out=False,
                       data=a,
                       hue='date',
                       size=8,
                       aspect=1.5).set(xlim=(0, 2000), xlabel='raw GFP')
            if save:
                plt.savefig(op.join(datadir,
                                    "{} vs raw.png".format(labeldic[p])))

            g = sns.lmplot(x, y,
                           fit_reg=False,
                           data=a,
                           col='date',
                           hue='date',
                           col_wrap=3).set(xlim=(0, 2000),
                                           ylim=(0., 1.),
                                           xlabel='raw GFP')
            label_n(g, r2, Rsqr=True)

            if save:
                plt.savefig(op.join(datadir,
                                    "{} vs raw by date.png".format(
                                        labeldic[p])))


def plotmom_budfp(**kwargs):
    """
    plot Δψ of mom and first bud point
    """
    COL_ODR = kwargs.get('COL_ODR')
    df = kwargs.get('data_mfp')
    datadir = kwargs.get('savefolder')
    save = kwargs.get('save', False)
    sns.set_style('whitegrid')
    colnames = {'cellaxis_mom_budfp': 'mom position'}
#                'DY': u'Δψ scaled'}
    df = df.rename(columns=colnames)
    with sns.plotting_context('talk'):
        p = sns.FacetGrid(df,
                          col="type",
                          hue='type',
                          col_wrap=4,
                          col_order=COL_ODR,
                          size=3,
                          aspect=1.4)
        p = p.map(sns.pointplot, 'mom position', 'DY')
        s = df['mom position'].cat.categories.astype('str').tolist()
        s[:-1] = [str(float(i)/10.) for i in s[:-1]]
        p.set_xticklabels(labels=s)
        if save:
            p.savefig(op.join(datadir, "DY_mom_firstbud_facetted.png"))
