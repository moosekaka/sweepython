# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 15:42:01 2015
grwoth curves
@author: sweel
"""
import matplotlib.pylab as plt
import pandas as pd
import seaborn as sns
import numpy as np
from scipy import stats as sp
sns.set_context("talk")
sns.set(style="whitegrid")

plt.close('all')

gc1 = pd.read_csv('gc1.csv',
                  parse_dates=[0])
gc2 = pd.read_csv('gc2.csv',
                  parse_dates=[0])
gc3 = pd.read_csv('gc3.csv',
                  parse_dates=[0])

colors = {'YPL': sns.xkcd_rgb["pale red"],
          'YPE': sns.xkcd_rgb["medium green"]}


def doubletime(i, gc, initial, end):
    gc['ETime'] = gc.ix[:, 'Time'].apply(lambda x: x - gc1.ix[0].Time)
    gc['ETime'] = gc['ETime'].astype('timedelta64[s]') / 3600.
    data = gc.ix[initial:, ['ETime', i]]
    gcy0 = data.ix[initial, i]
    data['ln'] = data.ix[initial:end, i].apply(lambda x: np.log(x))
    data['ln'] = data['ln'] - np.log(gcy0)
    slope, _, r, p, std_err = sp.linregress(data.ix[initial:end, 'ETime'],
                                            data.ix[initial:end, 'ln'])
    with sns.plotting_context('talk', font_scale=1.4):
        fig, ax = plt.subplots(1, 1)
        sns.regplot(x='ETime',
                    y='ln',
                    data=data.dropna(),
                    color=colors[i],
                    ax=ax,
                    scatter_kws={'s': 70})
        ax.text(0.4, 0.15,
                'slope = %6.4f /hr\ndoubling time = %4.3f hrs\nr =%7.4f' %
                (slope, np.log(2) / slope, r),
                fontsize=24, ha='left', va='center',
                transform=ax.transAxes)
        ax.set_title(i)

    return data, np.log(2) / slope, r


with sns.plotting_context('talk', font_scale=1.4):
    fig, ax = plt.subplots(1, 1)
    fig1, ax1 = plt.subplots(1, 1)
    fig2, ax2 = plt.subplots(1, 1)
    gc1.ix[:, :3].plot(marker='o',
                       ax=ax,
                       )
    gc2.ix[:, :3].plot(marker='o',
                       ax=ax1,
                       )
    gc3.ix[:, :3].plot(marker='o',
                       ax=ax2,
                       )
    ax.set(title='samp A')
    ax1.set(title='samp B')

print 'doubling time for YPL: %6.4f hrs (r=%6.4f)' % doubletime('YPL', gc1, 2, 7)[1:]
print 'doubling time for YPE: %6.4f hrs (r=%6.4f)' % doubletime('YPE', gc1, 2, 9)[1:]
print 'doubling time for YPL: %6.4f hrs (r=%6.4f)' % doubletime('YPL', gc2, 3, 8)[1:]
print 'doubling time for YPE: %6.4f hrs (r=%6.4f)' % doubletime('YPE', gc2, 3, 9)[1:]
