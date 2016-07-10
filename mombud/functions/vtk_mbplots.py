# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 16:10:16 2016
Module for plots of analysis of mother bud function in budding yeast
"""
import sys
import os
import os.path as op
import inspect
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from wrappers import UsageError
#from mombud.vtk_mom_bud_analyse_refactored import postprocess_df, _dySet
# pylint: disable=C0103

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.close('all')
error_text = ("labelhandler() called with type {}, which does not exist for"
              " this handler. Defaulting to unlabelled title for subplot")

class labelhandler(object):
    """ Callable object for labeling subplot titles"""
    def __init__(self, htype='normal'):
        self._htype = htype
        self._fundict = {'normal': 'fun1',
                         'facet': 'fun2', 'rowfacet': 'fun2', 'rsqr': 'fun2'}

    def __call__(self, handle_to_ax, label_dict):
        try:
            functype = self._fundict[self._htype]
            getattr(self, functype)(handle_to_ax, label_dict)
        except KeyError as e:
            print error_text.format(e)

    def fun1(self, handle, labeldic):
        """
        for normal plots, i.e. 'unfacetted'
        """
        assert self._htype == 'normal'  # sanity check here
        labels = [xtik.get_text().strip() for xtik in handle.axes.get_xticklabels()]
        new_labels = ['{}\n N={}'
                      .format(old_lab,
                              labeldic[old_lab]) for old_lab in labels]
        handle.axes.set_xticklabels(new_labels)

    def fun2(self, handle, labeldic):
        """
        for facetted plots
        """
        for ax in handle.axes.flat:
            oldtitle = ax.get_title()

            if self._htype == 'rowfacet':  # multi-row facetted plots"
                media = (oldtitle  # carbon/media type label
                         .split('|')[0]
                         .split('=')[1]
                         .strip())
                budv = float(oldtitle   # budvol label
                             .split('=')[-1]
                             .strip())
                newtitle = ('{}, N = {}'
                            .format(media,
                                    labeldic.xs(media).get([budv])[0]))
                ax.set_title(newtitle)

            else:
                oldtitle = oldtitle.split('=')[1].strip()
                if self._htype != 'rsqr':  # defaults for "facet"
                    ax.set_title('{}, N={}'
                                 .format(oldtitle,
                                         labeldic[oldtitle]))
                else:
                    ax.set_title('{}, R^2={:5.3f}'  # adds R^2 corr. labels"
                                 .format(oldtitle,
                                         labeldic[oldtitle]))


class plviol(object):
    """
    Wrapper class for generating violinplots
    """
    # class/wrapper level attr.
    def_init_kwargs = ('col_order', 'default_ylims', 'labeller')

    def __init__(self, plt_type='violinplot', nofig=False, **kwargs):
        self.pltobj = getattr(sns, plt_type)  # will be a plt axes obj
        if not nofig:  # toggle to disable subplot in case of FacetGrid
            self.fig, self.ax = plt.subplots(1, 1)
        for name in type(self).def_init_kwargs:
            try:
                setattr(self, name, kwargs[name])
            except KeyError:
                continue
        self.data = None
        self.x = None
        self.y = None
        self.n_counts = None
        self.ylim = None

    def get_ylims(self, df):
        """
        get the ylims for the long form data set
        """
        try:
            self.default_ylims
            y_lims = (df.iloc[:, -1]
                      .dropna()
                      .quantile(self.default_ylims)
                      .tolist())
        except AttributeError:
            print "default_ylims not specified, returning minmax instead"
            y_lims = (df.iloc[:, -1].min(), df.iloc[:, -1].max())
        return tuple(y_lims)

    def get_group_counts(self, **kwargs):
        """
        Get counts of groups based on the first columns name or a groupkey
        """
        groupkey = kwargs.get('group_key', self.data.columns[0])

        self.n_counts = (self.data.groupby([groupkey, 'variable'])
                         .size()
                         .xs(self.data['variable'].iloc[0], level='variable'))

    def call_sns_plotter(self, order, ylims, **kwargs):
            """
            Draws the axes obj using the seaborn method self.pltobj
            """
            h = self.pltobj(data=self.data, ax=self.ax,
                            order=order, **kwargs)
            h.set_ylim(ylims[0], ylims[1])
            h.set_title(kwargs.get('title'))

            try:
                self.labeller(h, self.n_counts)
            except AttributeError:
                "Skipping labels of categorical vars."
                pass

    def plt(self, data=None, ylim=None, **kwargs):
        """
        Parameters
        ----------
        data : DataFrame
            data for plot input
        ylim : [None | 'auto' | tuple]
            `ylim` is an axes level parameter to manually set the y-limits,\
            as the axes level seaborn plot functions do not have a `ylim`\
            parameter, unlike FacetGrid. If set to `auto`, ylim is calculated
            according to a default_ylims parameter during object instantiation
        """
        self.data = data
        try:
            assert self.data is not None
        except AssertionError:
            etype, val, tb = sys.exc_info()
            raise UsageError("{} on {}".format(etype.__name__, tb.tb_lineno))

        # get the y-axis limits
        if ylim is not None:
            if ylim == 'auto':
                self.ylim = self.get_ylims(self.data)
            else:
                self.ylim = ylim

        # get the counts for the categorical vars
        self.get_group_counts(**kwargs)

        # if a subset of data is chosen, make col_order a subset too
        try:
            col_order_sub = [i for i in self.col_order if i in self.n_counts]
        except AttributeError:
            col_order_sub = None

        # tries to plot using the relevant seaborn plot type, except if this
        # instance's  plt() is overriden
        try:
            self.call_sns_plotter(order=col_order_sub,
                                  ylims=self.ylim, **kwargs)
        except AttributeError:
            pass

    def save_figure(self, savepath=None):
        if savepath is None:
            savepath = op.join(os.getcwd(), '_fig.png')
        self.fig.savefig(savepath)

    def turn_off_legend(self):
        try:
            self.ax.legend_.remove()
        except AttributeError:
            pass


class plbox(plviol):
    """
    subclass of plviol to have boxplot mean instead of violinplot medians
    """

    def plt(self, **kwargs):
        """
        extends plviol.plt to have boxplot mean markers
        """
        # plt.boxplot doesnt have kwargs, hence need to check kws allowed
        allowable = inspect.getargspec(sns.boxplot).args
        kws = {k: v for k, v in kwargs.iteritems() if k in allowable}
        sns.boxplot(ax=self.ax,
                    showmeans=True,
                    showbox=False,
                    showcaps=False,
                    showfliers=False,
                    order=self.col_order,
                    medianprops={'linewidth': 0},
                    whiskerprops={'linewidth': 0},
                    meanprops={'marker': '_',
                               'c': 'w',
                               'ms': 5,
                               'markeredgewidth': 2}, **kws)
        # plt this after boxplot so that labelling macros are preserved
        super(plbox, self).plt(**kwargs)


class plfacet(plviol):
    def __init__(self, **kwargs):
        # nofig: ensure that the base class does not call plt.figure
        super(plfacet, self).__init__(nofig=True, **kwargs)
        self.facet_obj = sns.FacetGrid
        self.ylim = None

    def get_group_counts(self, **kwargs):
        """
        Get counts of groups based on the first columns name or a groupkey
        """
        groupkey = kwargs.get('group_key', self.data.columns[0])
        try:
            self.n_counts = (self.data.groupby([groupkey, 'variable'])
                             .size()
                             .xs(self.data['variable'].iloc[0],
                                 level='variable'))
        except KeyError:
            # reshape data using pivot table
            dfcopy = self.data.reset_index()
            c1, c2, c3, c4, c5 = dfcopy.columns
            self.n_counts = (pd.pivot_table
                             (dfcopy, c1, index=c2, columns=c4, aggfunc=len)
                             ).iloc[:, 0]

    def plt(self, data=None, mapargs=None, **kwargs):
        super(plfacet, self).plt(data=data, **kwargs)  # consumes ylim

        # add args for set function of FacetGrid
        setargs = kwargs.pop('setargs', None)
        allowable = inspect.getargspec(sns.FacetGrid.__init__).args
        kws = {k: v for k, v in kwargs.iteritems() if k in allowable}
        kws.update({'ylim': self.ylim})
        print "{}\n".format(kws)
        self.facet_obj = sns.FacetGrid(data, **kws)

        try:
            self.facet_obj.map(self.pltobj, *mapargs).set(**setargs)
        except TypeError:
            try:
                self.facet_obj.map(self.pltobj, *mapargs)
            except TypeError:
                raise UsageError("Missing kwargs `mapargs` on line {}"
                                 .format(sys.exc_info()[-1].tb_lineno))

        try:
            self.labeller(self.facet_obj, self.n_counts)
        except AttributeError as e:
            print "Error: {}, labelling skipped.".format(e)
            pass

    def save_figure(self, savepath=None):
        if savepath is None:
            savepath = op.join(os.getcwd(), '_fig.png')
        self.facet_obj.savefig(savepath)


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


#def main():
#    try:
#        def_args = {'regen': False,
#                    'save': False,  # toggle to save plots
#                    'inpdatpath': 'celldata.pkl',
#                    'mbax': np.linspace(0., 1., 6),  # pos. along mom/bud cell
#                    'cellax': np.linspace(0, 1., 11),  # position along whole cell
#                    'binsvolbud': np.linspace(0, 40, 5),  # vol binning for bud
#                    'binsvolmom': np.array([0, 30, 40, 80.]),
#                    'gfp_plot_vars': ['DY_abs_mean_mom',
#                                      'DY_abs_mean_bud',
#                                      'DY_abs_cell_mean'],
#                    'COL_ODR': ['MFB1', 'NUM1', 'YPT11',
#                                'WT', 'YPE', 'YPL', 'YPR'],
#                    'HUE_ODR': ['DY_abs_mean_mom',
#                                'DY_abs_mean_bud',
#                                'whole_cell_abs']}
#
#        dydict = _dySet('DY')
#        def_args.update(dydict)
#
#        outputargs = postprocess_df(**def_args)  # calls inputdata(), mungedata()
#
#        kwargs = def_args
#        df = outputargs.get('data')
#        datadir = kwargs.get('savefolder')
#        col_labels = kwargs.get('gfp_plot_vars')
#        save = kwargs.get('save', False)
#        N = kwargs.get('counts')
#        COL_ODR = kwargs.get('COL_ODR')
#        frac = pd.melt(df,
#                       id_vars=['media'],
#                       value_vars=['frac'])
#        mb_dy = pd.melt(df,
#                        id_vars=['media', 'date'],
#                        value_vars=col_labels[:])
#        diameters = pd.melt(df,
#                            id_vars='media',
#                            value_vars=['bud_diameter', 'mom_diameter'])
#
#
#    #        h = (h.map(sns.violinplot, "value")
#    #             .set(xlim=(0.), xlabel='diameter/microns'))
#
#        m = ['NUM1', 'MFB1', 'YPT11', 'WT']
#
#        mutants = mb_dy[mb_dy['media'].isin(m)].reset_index(drop=True)
#        gr = mutants.groupby('date')
#        outkws2 = dict(default_ylims=[0.05, 0.95],
#                       labeller=labelhandler('normal'))
#        set1 = dict(data=mutants,
#                    x='date', y='value', hue='media', no_legend=False, # save=True,
#                    group_key='date')
#
#        plv = plviol(**outkws2)
#        plv.plt(**set1)
#
#        set2 = dict(data=frac,
#                    x='media', y='value', hue='media',# save=True,
#                    group_key='media', inner=None, no_legend=True)
#        plv2 = plbox(**outkws2)
#        plv2.plt(**set2)
#
#        set3 = dict(data=diameters, col_wrap=4, mapargs=['value',],
#                    col='media', hue='variable',# save=True,
#                    size=3, aspect=1.5,
#                    setargs=dict(xlim=(0.), xlabel='diameter'))
#        plv3 = plfacet()
#        plv3.plt(**set3)
#
#        print "Finished plotting!"
#        return 0
#
#    except UsageError as e:
#        print e
#        return 1
#
#
#if __name__=='__main__':
#    sys.exit(main())


#    a = plfacet()
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

