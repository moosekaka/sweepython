import pandas as pd
import matplotlib.pyplot as plt
import cPickle as pickle
import seaborn as sns
import numpy as np


def prolarea(mjr, mnr):
    e = np.sqrt(1-((mnr**2)/(mjr**2)))
    surfarea = 2 * np.pi * (mnr)**2 * (1 + mjr/(mnr*e) * np.arcsin(e))
    return surfarea


def prosimparea(mjr, mnr):
    A = (mjr*mnr)**1.6075
    B = (mnr*mnr)**1.6075
    srfarea = 4 * np.pi * ((A + A + B)/3)**(1/1.6075)
    return srfarea


Data = pd.read_table('Results.txt')

df = pd.DataFrame()
df['cell'] = Data.loc[:, 'Label'].apply(lambda x: x.partition(':')[2])
df['Minor'] = Data.Minor*.055*.5
df['Major'] = Data.Major*.055*.5
df['Vol'] = 4/3 * np.pi * df.Major * df.Minor**2
df['Surf'] = prolarea(df.Major, df.Minor)
df['Area2'] = prosimparea(df.Major, df.Minor)
df = df.drop(['Major', 'Minor'], axis=1)
dffull = df
df = df.groupby('cell').sum()
df = df.reset_index()
df['cat'] = df.loc[:, 'cell'].apply(lambda x: x[:3])

fig1, g = plt.subplots(1, 1)
sns.violinplot(x='cat', y='Surf', data=df, ax=g)
sns.stripplot(x='cat', y='Surf', data=df, jitter=.05, ax=g)
plt.show()

fig2, h = plt.subplots(1, 1)
sns.violinplot(x='cat', y='Vol',data=df,ax=h)
sns.stripplot(x='cat', y='Vol', data=df,jitter=.05, ax=h)
plt.show()

with open('cellVolume.pkl', 'wb') as output:
    pickle.dump(df, output)