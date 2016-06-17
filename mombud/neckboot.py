# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 00:06:42 2016
Generate bootstraps and plot bootstrap data along neck region of mom/bud
"""
import cPickle as pickle
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

COL_ODR = ['MFB1', 'NUM1', 'YPT11', 'WT', 'YPE', 'YPL', 'YPR']

with open('boot.pkl', 'rb') as inp:
    data = pickle.load(inp)

with open('actual.pkl', 'rb') as inpt2:
    actual = pickle.load(inpt2)

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