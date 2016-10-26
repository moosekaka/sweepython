# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 16:10:16 2016
Module for plots of analysis of mother bud function in budding yeast
"""
import os
import os.path as op
import inspect
import traceback
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from wrappers import FalseException, UsageError

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
        labels = [xtik.get_text().strip()
                  for xtik in handle.axes.get_xticklabels()]
        new_labels = [u'{}\n N={}'
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
                newtitle = (u'{}, N = {}'
                            .format(media,
                                    labeldic.xs(media).get([budv])[0]))
                ax.set_title(newtitle)

            else:
                oldtitle = oldtitle.split('=')[1].strip()
                if self._htype != 'rsqr':  # defaults for "facet"
                    ax.set_title(u'{}, N={}'
                                 .format(oldtitle,
                                         labeldic[oldtitle]))
                else:
                    ax.set_title(u'{}, R^2={:5.3f}'  # adds R^2 corr. labels"
                                 .format(oldtitle,
                                         labeldic[oldtitle]))


class plviol(object):
    """
    Wrapper class for generating violinplots
    """
    # class/wrapper level attr.
    def_init_kwargs = ('col_order', 'default_ylims', 'labeller')

    def __init__(self, plt_type='violinplot', nofig=False, **kwargs):
        """
        Store a ref to a seaborn plot type method in `pltobj`.
        `nofig` toggle to disable subplot in case of FacetGrid.
        """
        self.pltobj = getattr(sns, plt_type)
        if not nofig:
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
            assert self.default_ylims
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
        try:
            groupkey = kwargs.get('group_key', self.data.columns[0])
            self.n_counts = (self.data
                             .groupby([groupkey, 'variable'])
                             .size()
                             .xs(self.data['variable']
                                 .iloc[0],
                                 level='variable'))
        except KeyError:
            try:
                dfcopy = self.data.reset_index()
                c1, c2, _, c4, _ = dfcopy.columns
                self.n_counts = (pd.pivot_table(
                    dfcopy, c1, index=c2,
                    columns=c4, aggfunc=len)).iloc[:, 0]
            except ValueError:
                traceback.print_stack(limit=4)
                print "Not enough columns to do group counts, skipping"

    def call_sns_plotter(self, order, ylims, **kwargs):
            """
            Draws the axes obj using the seaborn method self.pltobj
            """
            allowable = inspect.getargspec(self.pltobj).args
            kws = {k: v for k, v in kwargs.iteritems() if k in allowable}
            h = self.pltobj(data=self.data, ax=self.ax,
                            order=order, **kws)
            try:
                h.set_ylim(ylims[0], ylims[1])
                h.set_title(kwargs.get('title'))
            except TypeError:
                pass
            try:
                h.set(**kwargs['setargs'])
            except KeyError:
                pass

    def label_group_counts(self, obj):
        try:
            self.labeller(obj, self.n_counts)
        except (AttributeError, TypeError):
            print "Skipping labels of categorical vars \n"

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
            if self.data is None:
                raise FalseException
        except FalseException:
            traceback.print_stack(limit=4)
            raise UsageError("Must provide 'data' arg")

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
        except (AttributeError, TypeError):
            col_order_sub = None

        # tries to plot using the relevant seaborn plot type, except if this
        # instance's  plt() is overriden
        try:
            self.call_sns_plotter(order=col_order_sub,
                                  ylims=self.ylim, **kwargs)
        except AttributeError:
            pass
        else:
            self.label_group_counts(self.ax)

    def save_figure(self, savepath=None):
        """
        saves figure to path
        """
        if savepath is None:
            savepath = op.join(os.getcwd(), '_fig.png')
        self.fig.savefig(savepath)

    def turn_off_legend(self):
        """
        toggle to remove legend
        """
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
    """
    extends plviol to those with facet objects
    """
    def __init__(self, **kwargs):
        # nofig: ensure that the base class does not call plt.figure
        super(plfacet, self).__init__(nofig=True, **kwargs)
        self.facet_obj = sns.FacetGrid
        self.ylim = None

    def plt(self, data=None, mapargs=None, **kwargs):
        super(plfacet, self).plt(data=data, **kwargs)  # consumes ylim

        # add. args for FacetGrid.set() method
        setargs = kwargs.pop('setargs', None)
        allowable = inspect.getargspec(sns.FacetGrid.__init__).args
        kws = {k: v for k, v in kwargs.iteritems() if k in allowable}
        kws.update({'ylim': self.ylim})
        self.facet_obj = sns.FacetGrid(data, **kws)

        try:
            self.facet_obj.map(self.pltobj, *mapargs).set(**setargs)
        except TypeError:
            try:
                self.facet_obj.map(self.pltobj, *mapargs)
            except TypeError:
                traceback.print_stack(limit=4)
                raise UsageError("Missing kwargs `mapargs`")

        self.label_group_counts(self.facet_obj)

    def save_figure(self, savepath=None):
        if savepath is None:
            savepath = op.join(os.getcwd(), '_fig.png')
        self.facet_obj.savefig(savepath)
