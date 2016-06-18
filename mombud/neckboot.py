# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 00:06:42 2016
Generate bootstraps and plot bootstrap data along neck region of mom/bud
"""
import os
import os.path as op
import cPickle as pickle
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from mombud.functions import vtk_mbfuncs as vf
# pylint: disable=C0103


class wrapper(object):
    """
    wrapper class to autogenerate missing data
    """
    def __init__(self, **kwargs):
        self.celldfpath = kwargs.get('celldfdatapath')
        self.outdata = kwargs.get('datadict')
        self.vtkdf = vf.wrapper(inpdatpath=self.celldfpath)

    def gen_neck_boot(self, **kwargs):
        """
        generate bootstrap data
        """
        self.outdata['boot'] = bootNeck(self.vtkdf, save=True, **kwargs)

    def gen_neck_actual(self, save=True, **kwargs):
        """
        generate actual data
        """
        # note that vf.neckDY is modifying INPLACE the outdata dict!
        for key in self.vtkdf.keys():
            cell = self.vtkdf[key]['df']
            vf.neckDY(key, cell,
                      self.vtkdf[key]['neckpos'],
                      self.outdata['actual'])

        if save:
            with open(kwargs.get('neck_act_savepath',
                                 os.getcwd()), 'wb') as out:
                pickle.dump(self.outdata['actual'], out)


def bootNeck(vtkdf, dd=0.3, num_runs=10, save=False, **kwargs):
    """
    Parameters
    ----------
    vtkdf : DataFrame
        Ind cell data datafrmes
    dd : Float
        distance from neck
    num_runs : Int
        num of sample runs for bootstrap

    kwargs
    ------
    boot_savepath : Str
        filepath for saving bootstrap results, def. to curdir if not specified
    """

    merge = defaultdict(dict)
    for key in sorted(vtkdf.keys()):
        print "now on cell {}".format(key)
        neck_position = vtkdf[key]['neckpos']
        cell = vtkdf[key]['df']
        merge[key]['bud'] = {}
        merge[key]['mom'] = {}
        cleft = cell.reset_index(level='name')
        for i in range(num_runs):
            cright = cell.sample(n=cell.shape[0])
            cright = cright[['DY', 'DY_abs']].reset_index(level='name')
            c3 = cleft[['x', 'type']].merge(cright,
                                            left_index=True,
                                            right_index=True,
                                            indicator=True)
            merge[key]['bud'][i] = c3.loc[
                (c3.x >= neck_position) &
                (c3.x < (neck_position + dd))].DY.mean()
            merge[key]['mom'][i] = c3.loc[
                (c3.x < neck_position) &
                (c3.x >= (neck_position - dd))].DY.mean()

    if save:
        fpath = kwargs.get('boot_savepath', os.getcwd())
        with open(fpath, 'wb') as out:
            pickle.dump(merge, out)
    return merge


def plotNeck(**kwargs):
    """
    Plot Δψ at the bud neck region
    """
    df = kwargs.get('neckdata')
    datadir = kwargs.get('savefolder', os.getcwd())
    save = kwargs.get('save', False)

    df = df.reset_index()
    df['type'] = df.cellname.apply(lambda x: x.split('_')[0])
    actual = pd.melt(df[df['dist'] == 0.3],
                     id_vars=['type'],
                     value_vars=['bud', 'mom'])

    n_dist = pd.melt(df,
                     id_vars=['dist'],
                     value_vars=['bud', 'mom'])

    with sns.plotting_context('talk'):
        n_dist.dropna(inplace=True)
        _, ax1 = plt.subplots(1, 1)
        q1 = sns.barplot(x='dist',
                         y='value',
                         hue='variable',
                         data=n_dist,
                         ax=ax1)
        leg = q1.get_legend()
        plt.setp(leg, bbox_to_anchor=(0.85, 0.7, .3, .3))
        if save:
            plt.savefig(op.join(datadir, "neckregionDY.png"))

        _, ax2 = plt.subplots()
        sns.violinplot(x='type',
                       y='value',
                       hue='variable',
                       data=actual,
                       order=kwargs.get('COL_ODR'),
                       ax=ax2).set(title='Actual', ylim=(0, 0.95))


def plotBoot(gen_boot_data, **kwargs):
    """
    Plot bootstrap Δψ neckregion data
    """
    runs = []
    cellid = []

    for ids, d in gen_boot_data.iteritems():
        cellid.append(ids)
        runs.append(pd.DataFrame.from_dict(d))

    df = pd.concat(runs, keys=cellid)

    df2 = df[~df.bud.isnull()]
    df2.index.names = ['cellname', 'run_no']
    df2 = df2.reset_index()
    df2['type'] = df2.cellname.apply(lambda x: x.split('_')[0])
    df3 = df2.groupby(['type', 'cellname']).mean()
    df4 = df3.reset_index()
    boot = pd.melt(df4,
                   id_vars=['type'],
                   value_vars=['bud', 'mom'])

    with sns.plotting_context('talk'):
        _, ax1 = plt.subplots()
        sns.violinplot(x='type',
                       y='value',
                       hue='variable',
                       order=kwargs.get('COL_ODR'),
                       data=boot,
                       ax=ax1).set(title='Bootstrapped', ylim=(0, 0.95))


def main():
    """
    main
    """
    sns.set_style('whitegrid')
    data = defaultdict(dict)
    params = {'celldfdatapath': 'celldata.pkl',
              'num_runs': 100,
              'boot_savepath': 'neck_boot.pkl',
              'neck_act_savepath': 'neck_actual.pkl',
              'datadict': data,
              'COL_ODR': ['MFB1', 'NUM1',
                          'YPT11', 'WT',
                          'YPE', 'YPL', 'YPR']}
    datobj = wrapper(**params)

    for files in ['neck_boot', 'neck_actual']:
        keyname = files.split('_')[1]  # key for storing actual and boot data
        try:
            with open('{}.pkl'.format(files), 'rb') as inp:
                data[keyname] = pickle.load(inp)
        except IOError:
            # must pass **params in () !!
            getattr(datobj, 'gen_%s' % files)(**params)
            data[keyname] = getattr(datobj, 'outdata')[keyname]

    cell_id = []
    frames = []
    for ids, d in data['actual'].iteritems():
        cell_id.append(ids)
        frames.append(pd.DataFrame.from_dict(d, orient='columns'))
    pd.concat(frames, keys=cell_id)
    params['neckdata'] = pd.concat(frames, keys=cell_id)
    params['neckdata'].index.names = ['cellname', 'dist']
    params['neckdata'].reset_index(level='dist', inplace=True)

    plotNeck(**params)
    plotBoot(data['boot'], **params)

# _____________________________________________________________________________
if __name__ == '__main__':
    plt.close('all')
    main()
