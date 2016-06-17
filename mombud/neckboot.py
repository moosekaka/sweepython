# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 00:06:42 2016
Generate bootstraps and plot bootstrap data along neck region of mom/bud
"""
import sys
import cPickle as pickle
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import defaultdict
from mombud.functions import vtk_mbfuncs as vf
# pylint: disable=C0103

COL_ODR = ['MFB1', 'NUM1', 'YPT11', 'WT', 'YPE', 'YPL', 'YPR']
data = defaultdict(dict)
params = {'celldfdatapath': 'celldata.pkl',
          'boot_savepath': 'neck_boot.pkl',
          'neck_act_savepath': 'neck_actual',
          'datadict': data}


class wrapper:
    """
    wrapper class to autogenerate missing data
    """
    def __init__(self, **kwargs):
        self.celldfpath = kwargs.get('celldfdatapath')
        self.outdata = kwargs.get('datadict')
        self.vtkdf = vf.wrapper(inpdatpath=self.celldfpath)

    def gen_neck_boot(self, **kwargs):
        bootNeck(self.vtkdf, dd=0.3, num_runs=1, save=True, **kwargs)

    def gen_neck_actual(self, **kwargs):
#        note that vf.neckDY is modifying INPLACE the outdata dict!
        for key in self.vtkdf.keys():
            cell = self.vtkdf[key]['df']
            vf.neckDY(key, cell,
                      self.vtkdf[key]['neckpos'],
                      self.outdata['actual'])


def bootNeck(vtkdf, dd=0.3, num_runs=100, save=False, **kwargs):
    """
    Parameters
    ----------

    vtkdf : DataFrame
        Ind cell data datafrmes

    dd : Float
        distance from neck

    num_runs : Int
        num of sample runs for bootstrap

    Kwargs
    ------

    boot_savepath : Str
        filepath for saving bootstrap results

    """
    merge = defaultdict(dict)
    for k in sorted(vtkdf.keys()):
        print "now on cell {}".format(k)
        neck_position = vtkdf[k]['neckpos']
        cell = vtkdf[k]['df']
        merge[k]['bud'] = {}
        merge[k]['mom'] = {}
        cleft = cell.reset_index(level='name')
        for i in range(num_runs):
            cright = cell.sample(n=cell.shape[0])
            cright = cright[['DY', 'DY_abs']].reset_index(level='name')
            c3 = cleft[['x', 'type']].merge(cright,
                                            left_index=True,
                                            right_index=True,
                                            indicator=True)
            merge[k]['bud'][i] = c3.loc[
                (c3.x >= neck_position) &
                (c3.x < (neck_position + dd))].DY.mean()
            merge[k]['mom'][i] = c3.loc[
                (c3.x < neck_position) &
                (c3.x >= (neck_position - dd))].DY.mean()

    if save:
        fpath = kwargs.get('boot_savepath')
        with open(fpath, 'wb') as out:
            pickle.dump(merge, out)



DAT = wrapper(**params)

for k in ['neck_boot', 'neck_actual']:
    try:
        with open('{}.pkl'.format(k), 'rb') as inp:
            data[k] = pickle.load(inp)
    except IOError as e:
        getattr(DAT, 'gen_%s' % k )(**params)  # must pass **params in () !!



cell_id = []
frames = []
for ids, d in data['actual'].iteritems():
    cell_id.append(ids)
    frames.append(pd.DataFrame.from_dict(d, orient='columns'))
pd.concat(frames, keys=cell_id)
dicout['dfneck'] = pd.concat(frames, keys=cell_id)
dicout['dfneck'].index.names = ['cellname', 'dist']
dicout['dfneck'].reset_index(level='dist', inplace=True)
#with open('actual.pkl', 'rb') as inpt2:
#    actual = pickle.load(inpt2)

runs = []
cellid = []

for ids, d in data.iteritems():
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

df_actual = actual['neckres']
df_actual = df_actual.reset_index()
df_actual['type'] = df_actual.cellname.apply(lambda x: x.split('_')[0])
actual = pd.melt(df_actual[df_actual['dist'] == 0.3],
                 id_vars=['type'],
                 value_vars=['bud', 'mom'])

plt.close('all')
with sns.plotting_context('talk'):
    _, ax1 = plt.subplots()
    q1 = sns.violinplot(x='type',
                        y='value',
                        hue='variable',
                        order=COL_ODR,
                        data=boot,
                        ax=ax1).set(title='Bootstrapped', ylim=(0, 0.95))

    _, ax2 = plt.subplots()
    q2 = sns.violinplot(x='type',
                        y='value',
                        hue='variable',
                        data=actual,
                        order=COL_ODR,
                        ax=ax2).set(title='Actual', ylim=(0, 0.95))


#    with open('actual.pkl', 'wb') as out:
#        dicres = dict(zip(['cellres', 'momres', 'budres', 'neckres'],
#                          [cellall, cellposmom, cellposbud, neckregion]))
#        pickle.dump(dicres, out)

    # concatenate neckregion data to DataFrame
#    cell_id = []
#    frames = []
#    for ids, d in dicdf['neck'].iteritems():
#        cell_id.append(ids)
#        frames.append(pd.DataFrame.from_dict(d, orient='columns'))
#    pd.concat(frames, keys=cell_id)
#    dicout['dfneck'] = pd.concat(frames, keys=cell_id)
#    dicout['dfneck'].index.names = ['cellname', 'dist']
#    dicout['dfneck'].reset_index(level='dist', inplace=True)
