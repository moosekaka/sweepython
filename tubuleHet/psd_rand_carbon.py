# -*- coding: utf-8 -*-
"""
Calculates and plots POWER SPEC DENSITY (psd) of Δψ for real and random dists.
@author: sweel
"""
import os
import os.path as op
import cPickle as pickle
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import wrappers as wr
from tubuleHet.autoCor.AutoPopFunc import psd, conv_to_pd, tidy_psd
# pylint: disable=C0103
# pylint: disable=R0204
sns.set_context("talk")
sns.set(style="whitegrid")
sns.set(rc={"legend.markerscale": 3})
colors = ["medium green",
          "greyish blue",
          "yellowy brown",
          "reddish grey"]

# =============================================================================
#           Data initialization
# =============================================================================
plt.close('all')
datadir = op.join(os.getcwd(), 'data')
rawdir = op.join(os.getcwd(), 'output')
vtkF = wr.ddwalk(op.join(rawdir, 'normalizedVTK'),
                 '*skeleton.vtk', start=5, stop=-13)

PSDY = defaultdict(dict)  # DY_scaled  autocors
PSDP = defaultdict(dict)  # Shuffled dist autocors
PSDN = defaultdict(dict)  # Normal dist autocors
PSDU = defaultdict(dict)  # uniform dist autocors

psd_type = {'actual YPE': PSDY['YPE'],
            'normal': PSDN['YPE'],
            'shuffled': PSDP['YPE'],
            'uniform': PSDU['YPE']}

bins = np.linspace(0, .5, 22)  # bins for the x-axis (freq spectrum)

# =============================================================================
# Load fitted and real data, calculate PSD
# =============================================================================
for mtype in sorted(vtkF.keys())[:]:
    for cell in vtkF[mtype].keys():
        PSDU[mtype][cell] = []
        PSDN[mtype][cell] = []
        PSDP[mtype][cell] = []
        PSDY[mtype][cell] = []
        with open(op.join(rawdir,
                          'fitted_data_scaled',
                          '%s.pkl' % cell), 'rb') as inpt:
            (lNorm, lNormP, randNDY, randUDY, llineId) = pickle.load(inpt)

        PSDY[mtype][cell].append(psd(np.squeeze(lNorm), 40))
        PSDU[mtype][cell].append(psd(np.squeeze(randUDY), 40))
        PSDN[mtype][cell].append(psd(np.squeeze(randNDY), 40))
        PSDP[mtype][cell].append(psd(np.squeeze(lNormP), 40))
        print "done psd for %s" % cell


# =============================================================================
# Power spectrum density actual (YPE) vs random
# =============================================================================
psd_tidydata = pd.DataFrame()
for dist_type in sorted(psd_type.keys()):
    realranddata = psd_type[dist_type]
    psdx, psdy = conv_to_pd(realranddata)

    psd_tidydata = psd_tidydata.append(
        tidy_psd(psdx, psdy, bins, dist_type),
        ignore_index=True)

# =============================================================================
# Power spectrum density actual by media
# =============================================================================
psd_tidydata2 = pd.DataFrame()
for carbon in sorted(PSDY.keys()):
    carbondata = PSDY[carbon]
    psdx, psdy = conv_to_pd(carbondata)

    psd_tidydata2 = psd_tidydata2.append(
        tidy_psd(psdx, psdy, bins, carbon),
        ignore_index=True)

# ============================================================================
# Plot
# ============================================================================
sns.set(style="whitegrid")
# vs random
with sns.plotting_context('talk', font_scale=1.4):
    _, ax1 = plt.subplots(1, 1)
    sns.pointplot(x='u',
                  y='psd',
                  hue='type',
                  palette=sns.xkcd_palette(colors),
                  scale=.95,
                  data=psd_tidydata,
                  ax=ax1)
    ax1.get_legend().set_visible(False)
# vs carbon type
with sns.plotting_context('talk', font_scale=1.4):
    _, ax2 = plt.subplots(1, 1)
    sns.pointplot(x='u',
                  y='psd',
                  hue='type',
                  scale=.95,
                  data=psd_tidydata2,
                  ax=ax2)
    ax2.get_legend().set_visible(False)
