# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 16:10:16 2016
Module for plots of analysis of mother bud function in budding yeast
"""
import os
import os.path as op
import datetime
import inspect
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from wrappers import UsageError

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.close('all')
# pylint: disable=C0103


#set1 = dict(x='media', y='value', hue='variable', 'groupkey'='x',
#            data=gr.get_group('032016'), title='032016')
#set2 = dict(x='media', y='value', hue='variable',
#            data=gr.get_group('032716'), title='032716')


def get_group_counts(**kwargs):
    """
    Get counts of groups based on the first columns name or a groupkey
    """
    df = kwargs.get('data')
    groupkey = kwargs.pop('groupkey', df.columns[0])
    N = (df.groupby([groupkey, 'variable'])
         .size()
         .xs(df['variable'].iloc[0], level='variable'))
    return N


class plviol(object):
    """
    Wrapper class for generating violinplots
    """
    def __init__(self, plt_type='violinplot', **kwargs):
        self.plt_type = plt_type
        self.COL_ODR = kwargs.get('COL_ODR')
        self.q_lims = kwargs.get('q_lims')
        self.labeller = kwargs.get('labeller')
        self.ax = None

    def get_ylims(self, data):
        """
        get the ylims for the long form data set
        """
        if self.q_lims is not None:
            y_lims = data['value'].dropna().quantile(self.q_lims).tolist()
        else:
            y_lims = (data.value.min(), data.value.max())
        return y_lims

    def plt(self, ylims=None, **kwargs):
        """
        plot violinplots onto instance axes self.ax
        """
        data = kwargs.get('data')
        if ylims is None:  # if ylims is not specifice, get automatically
            ylims = self.get_ylims(data)

        N = get_group_counts(**kwargs)  # group counts

        # if COL_ODR is specified, get the subset of it
        if self.COL_ODR is not None:
            col_order = [i for i in self.COL_ODR if i in N]
        else:
            col_order = None

        _, self.ax = plt.subplots(1, 1)
        h = (getattr(sns, self.plt_type)(ax=self.ax,
                                         order=col_order, **kwargs))

        h.set_ylim(ylims[0], ylims[1])
        h.set_title(kwargs.get('title'))
        self.ax.legend(loc=0)
#        self.add_count_labels(h, N)
        if self.labeller is not None:
            self.labeller(h, N)


class plboxer(plviol):
    """
    subclass of plviol to have boxplot mean instead of violinplot medians
    """

    def bxplt(self, **kwargs):
        """
        extends plviol.plt to have boxplot mean markers
        """
        self.plt(**kwargs)
        # plt.boxplot doesnt have kwargs, hence need to check kws allowed
        allowable = inspect.getargspec(sns.boxplot).args
        kws = {k: v for k, v in kwargs.iteritems() if k in allowable}
        sns.boxplot(ax=self.ax,
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
                    **kws)


def ranger(lst, mult):
    """
    return array for y-ticks or x-ticks with multiples of `mult`
    from limits in `lst`
    """
    yt = [i - i % mult for i in lst]
    return np.arange(yt[0], yt[1], mult)


def getrval(df, x, y, labeldic):
    """
    return a subset DataFrame and R^2 values for columns x, y in original df
    """
    df = df[[x, y, 'date']].reset_index(drop=True)
    df.rename(columns=labeldic, inplace=True)
    pear = (df.groupby('date')
            .corr().xs(labeldic[y], level=1)
            ['raw'].to_dict())
    r_sqr = {key: value**2 for key, value in pear.iteritems()}

    return df, r_sqr


def labelhandler(htype):
    """
    modifies title on `handle` (an axes obj) to include additional
    labels in dict. `labeldic`
    """

    def fun1(handle, labeldic):
        """
        for normal plots
        """
        labels = [xl.get_text().strip()
                  for xl in handle.axes.get_xticklabels()]
        new_labels = ['{}\n N={}'
                      .format(old_lab,
                              labeldic[old_lab]) for old_lab in labels]
        handle.axes.set_xticklabels(new_labels)

    def fun2(handle, labeldic):
        """
        for single row facetted plots or rqr regressions
        """
        for ax in handle.axes.flat:
            oldtitle = ax.get_title()
            oldtitle = oldtitle.split('=')[1].strip()
            if not htype == 'rsqr':
                ax.set_title('{}, N={}'
                             .format(oldtitle,
                                     labeldic[oldtitle]))
            else:
                ax.set_title('{}, R^2={:5.3f}'
                             .format(oldtitle,
                                     labeldic[oldtitle]))

    def fun3(handle, labeldic):
        """
        for facetted grid plots with multi rows
        """
        for ax in handle.axes.flat:
            oldtitle = ax.get_title()

            media = (oldtitle
                     .split('|')[0]
                     .split('=')[1]
                     .strip())  # carbon/media type label
            budv = float(oldtitle
                         .split('=')[-1]
                         .strip())  # budvol label
            newtitle = ('{}, N = {}'
                        .format(media, labeldic.xs(media).get([budv])[0]))
            ax.set_title(newtitle)

    fundict = {'facet': fun2, 'normal': fun1, 'rowfacet': fun3, 'rsqr': fun2}

    return fundict[htype]

labelBudVol = labelhandler('rowfacet')
labelFacet = labelhandler('facet')
labelNormal = labelhandler('normal')
labelRsqr = labelhandler('rsqr')


def plotDims(**kwargs):
    """
    Diameters of cells
    """
    df = kwargs.get('data')
    save = kwargs.get('save', False)
    datadir = kwargs.get('savefolder')
    COL_ODR = kwargs.get('COL_ODR')
    diameters = pd.melt(df,
                        id_vars='media',
                        value_vars=['bud_diameter', 'mom_diameter'])

    with sns.plotting_context('talk'):
        g = sns.FacetGrid(diameters,
                          col='media',
                          col_wrap=4,
                          hue="variable",
                          col_order=COL_ODR,
                          size=3,
                          aspect=1.5)
        h = sns.FacetGrid(diameters,
                          col='media',
                          col_wrap=4,
                          hue="variable",
                          col_order=COL_ODR,
                          size=3,
                          aspect=1.5)
        g = (g.map(sns.stripplot, "value", jitter=0.1)
             .set(xlim=(0.), xlabel='diameter/microns'))

        h = (h.map(sns.violinplot, "value")
             .set(xlim=(0.), xlabel='diameter/microns'))

        if save:
            g.savefig(op.join(datadir, 'celldiameters_stripplot.png'))
            h.savefig(op.join(datadir, 'celldiameters_violin.png'))


def plotSizeDist(**kwargs):
    """
    Distribution of bud and mom volumes
    """
    df = kwargs.get('data')
    datadir = kwargs.get('savefolder')
    N = kwargs.get('counts')
    save = kwargs.get('save', False)
    COL_ODR = kwargs.get('COL_ODR')

    budvol = df.ix[:, ['budvol', 'media']]
    momvol = df.ix[:, ['momvol', 'media']]
    budvol['N'] = budvol.groupby('media').transform('count')
    momvol['N'] = momvol.groupby('media').transform('count')

    sns.set_style('whitegrid')
    with sns.plotting_context('talk', font_scale=1.1):
        g = sns.FacetGrid(budvol,
                          col='media',
                          col_wrap=4,
                          hue='media',
                          col_order=COL_ODR)
        g = (g.map(sns.distplot, "budvol")
             .set(xlim=(0.)))
        labelFacet(g, N)
        if save:
            g.savefig(op.join(datadir, 'budsize_dist.png'))

        h = sns.FacetGrid(momvol,
                          col='media',
                          col_wrap=4,
                          hue='media',
                          col_order=COL_ODR)
        h = (h.map(sns.distplot, "momvol")
             .set(xlim=(0.)))
        labelFacet(h, N)

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
                         id_vars=['media', 'binvol'],
                         var_name='mom axis position',
                         value_name=r'$\Delta\Psi$ scaled gradient',
                         value_vars=binsaxis.tolist())
    bigbinsmom = bigbinsmom.dropna()

    bigbinsbud = pd.melt(dfbud,
                         id_vars=['media', 'binvol'],
                         var_name='bud axis position',
                         value_name=r'$\Delta\Psi$ scaled gradient',
                         value_vars=binsaxis.tolist())
    bigbinsbud = bigbinsbud.dropna()

    ylims = (bigbinsbud[r'$\Delta\Psi$ scaled gradient']
             .quantile([0.1, 0.9]).tolist())

    sns.set_style('whitegrid')
    with sns.plotting_context('talk', font_scale=1.):
        h = sns.FacetGrid(bigbinsmom,
                          col='media',
                          hue='media',
                          col_wrap=4,
                          sharex=True,
                          col_order=COL_ODR)
        h = (h
             .map(sns.pointplot,
                  'mom axis position',
                  r'$\Delta\Psi$ scaled gradient')
             .set(ylim=tuple(ylims)))
        labelFacet(h, N)
        if save:
            h.savefig(op.join(datadir, 'mom_cell_dy.png'))

        m1 = sns.FacetGrid(bigbinsbud,
                           col='media',
                           hue='media',
                           col_wrap=4,
                           col_order=COL_ODR)

        m1 = (m1
              .map(sns.pointplot,
                   'bud axis position',
                   r'$\Delta\Psi$ scaled gradient')
              .set(ylim=tuple(ylims)))
        labelFacet(m1, N)
        if save:
            m1.savefig(op.join(datadir, 'bud_cell_dy.png'))

    # with facetting by budvol
    with sns.plotting_context('talk', font_scale=.9):
        ylims2 = (bigbinsbud[r'$\Delta\Psi$ scaled gradient']
                  .quantile([0.025, 0.975]).tolist())

        if (ylims[1] - ylims[0]) < 2:
            mult = 0.25
        else:
            mult = np.ceil((ylims2[1] - ylims2[0]) / 6)
        yt = ranger(ylims2, mult)

        m0 = sns.FacetGrid(bigbinsbud,
                           row='media',
                           col="binvol",
                           hue='media',
                           row_order=COL_ODR,
                           col_order=binsvolbud[1:])

        m0 = (m0.map(sns.pointplot,
                     'bud axis position',
                     r'$\Delta\Psi$ scaled gradient')
              .set(yticks=yt, ylim=tuple(ylims2)))
        labelBudVol(m0, Nbud)

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

    ylims = (df['frac'].dropna().quantile([0., 0.995]).tolist())

    with sns.plotting_context('talk'):
        _, ax2 = plt.subplots(1, 1)
        h = (sns.pointplot(x='bin_budprog',
                           y='frac',
                           hue='media',
                           data=df,
                           ax=ax2))
        h.set(ylim=tuple(ylims),
              title=u"Δψ vs bud progression\n ",
              xlabel="bud progression",
              ylabel=u"Δψ bud/Δψ mom")
        leg = h.get_legend()
        plt.setp(leg, bbox_to_anchor=(0.85, 0.7, .3, .3))
        if save:
            plt.savefig(op.join(datadir, "DY vs bud progression.png"))

        p = sns.FacetGrid(df,
                          col='media',
                          hue='media',
                          col_wrap=4,
                          col_order=COL_ODR)
        p = p.map(sns.pointplot, 'bin_budprog', 'frac').set(ylim=tuple(ylims))
        if save:
            p.savefig(op.join(datadir, "DY_bud_prog_facetted.png"))


def plotViolins(**kwargs):
    """
    Violinplots for frac Δψ, mom vs bud scaled and Δψ abs distr
    """
    df = kwargs.get('data')
    datadir = kwargs.get('savefolder')
    col_labels = kwargs.get('viol_plot_vars')
    save = kwargs.get('save', False)
    N = kwargs.get('counts')
    COL_ODR = kwargs.get('COL_ODR')

    frac = pd.melt(df,
                   id_vars=['media'],
                   value_vars=['frac'])

    mb_dy = pd.melt(df,
                    id_vars=['media'],
                    value_vars=col_labels[:-1])

    sns.set_style('whitegrid')
    with sns.plotting_context('talk'):
        _, ax4 = plt.subplots(1, 1)
        ylims = frac['value'].dropna().quantile([0.015, 0.985]).tolist()
        j = sns.violinplot(x='media',
                           y='value',
                           hue='media',
                           data=frac.dropna(),
                           order=COL_ODR,
                           inner=None,
                           ax=ax4)
        j.set_ylim(ylims[0], ylims[1])
        j.get_legend().set_visible(False)

        k = sns.boxplot(x='media',
                        y='value',
                        hue='media',
                        data=frac.dropna(),
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
        labelNormal(j, N)
        if save:
            plt.savefig(op.join(datadir, "violin_fracDY.png"))

        _, ax3 = plt.subplots(1, 1)
        ylims = mb_dy['value'].dropna().quantile([0., .9995]).tolist()
        h = sns.violinplot(x='media',
                           y='value',
                           hue='variable',
                           order=COL_ODR,
                           data=mb_dy.dropna(),
                           ax=ax3)
        h.set_ylim(ylims[0], ylims[1])
        h.get_legend().set_visible(False)
        labelNormal(h, N)
        if save:
            plt.savefig(op.join(datadir, "Violin Mom_Bud_DY.png"))


def plotGFP(**kwargs):
    """
    Violinplots for GFP uptake by date and carbon type
    """
    col_labels = kwargs.get('gfp_plot_vars')
    df = kwargs.get('data')
    N = kwargs.get('counts')
    Ndate = kwargs.get('counts_date')
    datadir = kwargs.get('savefolder')
    save = kwargs.get('save', False)
    COL_ODR = kwargs.get('COL_ODR')
    HUE_ODR = col_labels

    with sns.plotting_context('talk', font_scale=1.):
        gfp_date = pd.melt(df,
                           id_vars=['date'],
                           value_vars=col_labels)

        ylims = gfp_date['value'].dropna().quantile([0.005, 0.95]).tolist()

        _, ax7 = plt.subplots(1, 1)
        g = sns.violinplot(x='date',
                           y='value',
                           hue='variable',
                           hue_order=HUE_ODR,
                           data=gfp_date,
                           ax=ax7)
        leg = g.get_legend()
        plt.setp(leg,
                 bbox_to_anchor=(.75, .85, .1, .2))
        g.set_ylim(ylims[0], ylims[1])

        if Ndate is not None:
            labelNormal(g, Ndate)

        if save:
            plt.savefig(op.join(datadir, "Violin-GFP_by_date.png"))

        gfp_type = pd.melt(df,
                           id_vars=['media'],
                           value_vars=col_labels)

        _, ax8 = plt.subplots(1, 1)
        g = sns.violinplot(x='media',
                           y='value',
                           bw='scott',
                           hue='variable',
                           order=COL_ODR,
                           hue_order=HUE_ODR,
                           data=gfp_type,
                           ax=ax8)
        leg = g.get_legend()
        plt.setp(leg,
                 bbox_to_anchor=(.75, .85, .1, .2))
        g.set_ylim(ylims[0], ylims[1])

        if N is not None:
            labelNormal(g, N)

        if save:
            plt.savefig(op.join(datadir, "Violin-GFP_by_type.png"))


def plotRegr(**kwargs):
    """
    plot regression coeff of Δψ raw vs scaled
    """
    df = kwargs.get('data')
    datadir = kwargs.get('savefolder')
    save = kwargs.get('save', False)

    labeldic = {'DY_abs_cell_mean': 'raw',
                'DY_cell_mean': 'scaled cell mean',
                'DY_median_mom': 'scaled mom mean',
                'DY_median_bud': 'scaled bud mean'}

    with sns.plotting_context('talk', font_scale=1.):
        for p in ['DY_cell_mean', 'DY_median_mom', 'DY_median_bud']:
            a, r2 = getrval(df, 'DY_abs_cell_mean', p, labeldic)
            x, y, _ = a.columns
            (sns.lmplot(x, y,
                        fit_reg=False,
                        legend_out=False,
                        data=a,
                        hue='date',
                        size=8,
                        aspect=1.5)
             .set(xlim=(0, 2000), xlabel='raw GFP'))
            if save:
                plt.savefig(op.join(datadir,
                                    "{} vs raw.png".format(labeldic[p])))

            g = (sns.lmplot(x, y,
                            fit_reg=False,
                            data=a,
                            col='date',
                            hue='date',
                            col_wrap=3)
                 .set(xlim=(0, 2000), ylim=(0., 1.), xlabel='raw GFP'))

            labelRsqr(g, r2)

            if save:
                plt.savefig(op.join(datadir,
                                    "{} vs raw by date.png".format(
                                        labeldic[p])))


#def example_utils():
#    """
#    Just a few functions shared between all the examples. Ensures example plots are
#    all the same size, etc.
#
#    Usage
#    -----
#
#    example_utils.title(fig, 'fill/fill_between/stackplot: Filled polygons',
#                            y=0.95)
#    example_utils.label(ax, 'fill')
#    """
#    import matplotlib.pyplot as plt
#
#    def setup_axes():
#        fig, axes = plt.subplots(ncols=3, figsize=(6.5,3))
#        for ax in fig.axes:
#            ax.set(xticks=[], yticks=[])
#        fig.subplots_adjust(wspace=0, left=0, right=0.93)
#        return fig, axes
#
#    def title(fig, text, y=0.9):
#        fig.suptitle(text, size=14, y=y, weight='semibold', x=0.98, ha='right',
#                     bbox=dict(boxstyle='round', fc='floralwhite', ec='#8B7E66',
#                               lw=2))
#
#    def label(ax, text, y=0):
#        ax.annotate(text, xy=(0.5, 0.00), xycoords='axes fraction', ha='center',
#                    style='italic',
#                    bbox=dict(boxstyle='round', facecolor='floralwhite',
#                              ec='#8B7E66'))
#

