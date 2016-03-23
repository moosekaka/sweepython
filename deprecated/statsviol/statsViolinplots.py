"""
    Plot a bunch of graphs based on stats calculated from statsDataset
    data imported from statsResults.pkl and functions are defined in svf
"""
import math
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import scipy.stats as sp
import seaborn as sns
import statsviol.StatViolinFunctions as svf
import config
sns.set(style="white")

config.counter = 0  # reset counter each time this script is run
# pylint: disable=C0103
media = []
mitoLength = {}
mitoNodes = {}
mitoEdgeNorm = {}
mitoStdNorm = {}
mitoEdges = {}
mitoNodesAll = {}
centralityBP = {}
closeCentralityBP = {}
rfpWidthCorCoef = {}
cellNormDY = {}
kNN3 = {}
ClstCoef = {}
charPathL = {}
AvgDeg = {}
phi = {}
pk3 = {}
with open('statsResults.pkl', 'rb') as inpt:
    (media, mitoEdges, mitoLength, rfpWidthCorCoef,
     mitoEdgeNorm, mitoStdNorm, cellNormDY, cellNormDY_raw,
     mitoEdgeNorm_unscaled, mitoStdNorm_unscaled,
     centralityBP, closeCentralityBP, ClstCoef, charPathL,
     phi, pk3, kNN3, AvgDeg,
     centralityBP_W, closeCentralityBP_W, ClstCoef_W, charPathL_W, kNN3_W) =\
     pickle.load(inpt)
plt.close('all')
# pylint: enable=C0103

# =============================================================================
#           Begin plots
# =============================================================================

#   Mito volume
X = [
    [j*math.pi*.15**2 for j in mitoLength[i]]
    for i in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(r'Total Mito Volume per cell ($\mu m^{3}$)',
                    'TotalMitoVolume', X, media, MLTCOMRES)

#    Average Edge Skew DY
X = [
    [sp.skew(i) for i in mitoEdgeNorm[j]]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(r'Average Edge $\Delta \Psi$ Skewness',
                    'EdgeSkew', X, media, MLTCOMRES,
                    pltLims=(0., 2.), cut=.1)

#   Average edge length
X = [
    [np.mean(i) for i in mitoEdges[j]]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(r'Average Edge Length ($\mu m$ / per cell)',
                    'EdgeLen', X, media, MLTCOMRES)

#   Average num edges
X = [
    [len(i) for i in mitoEdges[j]]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(r'Average Number Edges (per cell)',
                    'EdgeNums', X, media, MLTCOMRES)

#   Mean Edge Intensity DY
X = [
    [np.mean(i) for i in mitoEdgeNorm[j]]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(r'Mean Edge Intensity $\Delta \Psi$ (per cell)',
                    'EdgeDY', X, media, MLTCOMRES, cut=0.)

#   Edge Std Dev Intensity DY
X = [
    [np.mean(i) for i in mitoStdNorm[j]]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(r'Mean St Dev Edge Intensity $\Delta \Psi$ (per cell)',
                    'EdgeSTD', X, media, MLTCOMRES)

#   Mean Edge Intensity DY raw
X = [
    [np.mean(i) for i in mitoEdgeNorm_unscaled[j]]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(
    r'Mean Edge Intensity $\Delta \Psi$ Unscaled (per cell)',
    'EdgeDY_raw', X, media, MLTCOMRES, cut=0.)

#   Edge Std Dev Intensity DY raw
X = [
    [np.mean(i) for i in mitoStdNorm_unscaled[j]]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(
    r'Mean St Dev Edge Intensity $\Delta \Psi$ Unscaled (per cell)',
    'EdgeSTD_raw', X, media, MLTCOMRES)

#   Coef of variation for Edge DY
X = [
    [sp.variation(i) for i in mitoEdgeNorm[j]]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(
    r'Coef. of Var. Edge Intensity $\Delta \Psi$',
    'EdgeCoefVar', X, media, MLTCOMRES,
    pltLims=(0., 1.))

#   Coef of variation for Edge DY Raw
X = [
    [sp.variation(i) for i in mitoEdgeNorm_unscaled[j]]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(
    r'Coef. of Var. Edge Intensity $\Delta \Psi$ Unscaled',
    'EdgeCoefVar_raw', X, media, MLTCOMRES,
    pltLims=(0., 1.))

#   Average Degree
X = [AvgDeg[j][0] for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('Average Deg (2E/N)',
                    'AveDeg', X, media, MLTCOMRES)

#   Width Cor Coef with RFP inten
X = [
    [el[0] for el in rfpWidthCorCoef[j] if el[1] < 0.05]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('CorCoef RFP skel inten vs tubule Width',
                    'TubeWidthCorr', X, media, MLTCOMRES,
                    pltLims=(0., 1.))

#   Mean population DY
X = [
    [i for lis in mitoEdgeNorm[j] for i in lis]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(r'Mean Population $\Delta \Psi$ ',
                    'PopulationDY', X, media, MLTCOMRES, cut=0.)

#   Mean population DY Raw
X = [
    [i for lis in mitoEdgeNorm_unscaled[j] for i in lis]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(
    r'Mean Population $\Delta \Psi$ Unscaled',
    'PopulationDY_raw', X, media, MLTCOMRES, cut=0.)

#   Mean Cell DY
X = [cellNormDY[j] for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(r'Mean Whole Cell $\Delta \Psi$',
                    'CellDY', X, media, MLTCOMRES, cut=0.)

#   Mean Cell DY raw
X = [cellNormDY_raw[j] for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln(r'Mean Whole Cell $\Delta \Psi$ Unscaled',
                    'CellDY_raw', X, media, MLTCOMRES, cut=0.)

#   Phi
X = [phi[j] for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('Fraction of Edges of Largest Connected Component',
                    'Phi', X, media, MLTCOMRES,
                    pltLims=(0.45, 1.))

#   pk3
X = [pk3[j] for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('Fraction of Branchpoints (per cell)',
                    'Pk3', X, media, MLTCOMRES,
                    pltLims=(0., 1.))

#   Clustering Coef
X = [ClstCoef[j] for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('Clustering Coefficient',
                    'ClstCoef', X, media, MLTCOMRES, cut=0.)

#   Clustering Coef (weighted)
X = [ClstCoef_W[j] for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('Clustering Coefficient (weighted)',
                    'ClstCoefW', X, media, MLTCOMRES, cut=0.)

#   Char Path Lth
X = [charPathL[j] for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('Avg Shortest Path Length',
                    'CharPathL', X, media, MLTCOMRES)

#   Char Path Lth (weighted)
X = [charPathL_W[j] for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('Avg Shortest Path Length (weighted)',
                    'CharPathLW', X, media, MLTCOMRES)

#   Average Cell Betweeness Centrality
X = [
    [np.mean(el) for el in centralityBP[j]]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('Cell Mean Betweeness centrality',
                    'BtwCentral', X, media, MLTCOMRES, cut=0.)

#   Average Cell Betweeness Centrality (weighted)
X = [
    [np.mean(el) for el in centralityBP_W[j]]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('Cell Mean Betweeness centrality (weighted)',
                    'BtwCentralW', X, media, MLTCOMRES, cut=0.)

#   Average Cell Closeness Centrality
X = [
    [np.mean(el) for el in closeCentralityBP[j]]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('Cell Mean Closeness centrality',
                    'CloseCentral', X, media, MLTCOMRES)

#   Average Cell Closeness Centrality (weighted)
X = [
    [np.mean(el) for el in closeCentralityBP_W[j]]
    for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('Cell Mean Closeness centrality (weighted)',
                    'CloseCentralW', X, media, MLTCOMRES)

#   Bpts avg Deg Connectivity
X = [kNN3[j] for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('Branchpoints Avg Nearest Neighbr Connectivity',
                    'K_NearestNeighb', X, media, MLTCOMRES)

#   Bpts avg Deg Connectivity (weighted)
X = [kNN3_W[j] for j in media]
MLTCOMRES = svf.mltcomp(svf.mkdataset(X, media))
PLTOBJ = svf.pltvln('Branchpoints Avg Nearest Neighbr Connectivity (weighted)',
                    'K_NearestNeighbW', X, media, MLTCOMRES)
