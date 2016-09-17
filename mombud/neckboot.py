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


class GenData(object):
    """
    wrapper class to autogenerate missing data
    """

    def __init__(self, **kwargs):
        self.celldfpath = kwargs.get('celldfdatapath')
        self.outdata = kwargs.get('datadict')
        self.vtkdf = vf.gen_data(inpdatpath=self.celldfpath)

    def gen_neck_boot(self, save=True, **kwargs):
        """
        generate bootstrap data
        """
        self.outdata['boot'] = bootNeck(self.vtkdf, **kwargs)

        if save:
            fpath = kwargs.get('boot_savepath', os.getcwd())
            with open(fpath, 'wb') as out:
                pickle.dump(self.outdata['boot'], out)

    def gen_neck_actual(self, save=True, **kwargs):
        """
        generate actual data
        """
        # note that vf.neckDY is modifying INPLACE the outdata dict!
        for key in self.vtkdf.keys():
            cell = self.vtkdf[key]['df']
            vf.neckDY(key, cell,
                      self.vtkdf[key]['celldata']['neckpos'],
                      self.outdata['actual'])

        if save:
            with open(kwargs.get('neck_act_savepath',
                                 os.getcwd()), 'wb') as out:
                pickle.dump(self.outdata['actual'], out)


def bootNeck(vtkdf, dd=0.3, num_runs=10, save=False, **kwargs):
    """
    Bootstrap points along cellaxis for a given radius around a neckposition

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
    maskMom = vf.makeMask('mom')
    maskBud = vf.makeMask('bud')

    for key in sorted(vtkdf.keys()):  # cellname
        print "now on cell {}".format(key)
        neck_position = vtkdf[key]['celldata']['neckpos']
        cell = vtkdf[key]['df']

        for i in ['bud', 'mom', 'fk_bud', 'fk_mom']:
            merge[key][i] = {}

        cleft = cell[['x', 'type']].reset_index(level='name')
        for i in range(num_runs):
            # bootstrap the pixel positions along the cell axis
            cright = cell.sample(n=cell.shape[0])  # bootstrapped DY vector
            cright = cright['DY'].reset_index(level='name', drop=True)
            cleft['DY'] = cright
            merge[key]['bud'][i] = (cleft.loc[
                                    maskBud(cleft, neck_position, dd)]['DY']
                                    .mean())
            merge[key]['mom'][i] = (cleft.loc[
                                    maskMom(cleft, neck_position, dd)]['DY']
                                    .mean())

            # take random positions along the cell axis as the neck
            random_x = cell.sample(n=1)['x'].values[0]
            merge[key]['fk_bud'][i] = (cell
                                       .loc[maskBud(cell, random_x, dd)]['DY']
                                       .mean())
            merge[key]['fk_mom'][i] = (cell
                                       .loc[maskMom(cell, random_x, dd)]['DY']
                                       .mean())
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
#
        _, ax2 = plt.subplots()
        sns.violinplot(x='type',
                       y='value',
                       hue='variable',
                       data=actual,
                       order=kwargs.get('COL_ODR'),
                       ax=ax2).set(title='Actual', ylim=(0, 0.95))
    return actual


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
                   value_vars=['bud', 'mom', 'fk_bud', 'fk_mom'])

    with sns.plotting_context('talk'):
        _, ax1 = plt.subplots()
        sns.violinplot(x='type',
                       y='value',
                       hue='variable',
                       order=kwargs.get('COL_ODR'),
                       data=boot,
                       ax=ax1).set(title='Bootstrapped', ylim=(0, 0.95))
    return boot

# def main(**kwargs):
#    """
#    main
#    """
plt.close('all')

kwargs = dict(run_boot=False, num_runs=35)
sns.set_style('whitegrid')
data = defaultdict(dict)
params = {  # 'celldfdatapath': 'celldata.pkl',
    'celldfdatapath': 'filtered_neckbootdf.pkl',
    'num_runs': 1,
    'boot_savepath': 'neck_boot.pkl',
    'neck_act_savepath': 'neck_actual.pkl',
    'datadict': data,
    'COL_ODR': ['MFB1', 'NUM1',
                      'YPT11', 'WT',
                      'YPE', 'WT_COMBINED', 'YPL', 'YPR']}
params.update(kwargs)
datobj = GenData(**params)

run_boot = params.get('run_boot', False)

# if run_boot switch is False, will try to read in the pickle files below;
# except if the file is not found, will regenerate the files anyway
while not run_boot:
    for files in ['neck_boot', 'neck_actual']:
        keyname = files.split('_')[1]

        try:
            with open('{}.pkl'.format(files), 'rb') as inp:
                data[keyname] = pickle.load(inp)

        except IOError as e:
            print "{} not found, will regenerate".format(e.filename)
            run_boot = True
            break

    break

if run_boot:
    for files in ['neck_boot', 'neck_actual']:
        keyname = files.split('_')[1]
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

actual = plotNeck(**params)
actual = actual.replace(
    {'variable': {'bud': 'actual_bud', 'mom': 'actual_mom'}})
boot = plotBoot(data['boot'], **params)
combined = pd.concat([actual, boot])
wt = combined[combined.type.isin(['YPE', 'WT'])].copy()
wt.loc[:, 'type'] = 'WT_COMBINED'
combined = pd.concat([combined, wt])
c2 = combined[combined.variable.isin(['actual_mom', 'actual_bud'])]
with sns.plotting_context('talk', font_scale=1.35):
    _, ax3 = plt.subplots()
    sns.boxplot(x='type',
                y='value',
                hue='variable',
                data=combined,
                order=params.get('COL_ODR'),
                notch=True,
                ax=ax3).set(title='Neck DY')

    _, ax4 = plt.subplots()
    sns.boxplot(x='type',
                y='value',
                hue='variable',
                data=c2,
                order=params.get('COL_ODR'),
                notch=True,
                ax=ax4).set(title='Neck DY')
