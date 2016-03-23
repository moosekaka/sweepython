# =============================================================================
# To be able to read csv formated files, we will first have to import the
#   csv module.
# =============================================================================

import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as sp
import matplotlib.gridspec as gridspec

MassperMil = {}
CellCount = {}
O2Rate = {}
O2_500perMil = {}
OD2 = {}
CountperOD = {}
X = {}
Y = {}
W = {}
OD500 = {}
O2_500 = {}

media = [
    'ORI_YPD',
    'TUB1_YPD',
    'ORI_YPYE',
    'TUB1_YPYE']

colors = {
    'ORI_YPD': 'b',
    'ORI_YPYE': 'g',
    'TUB1_YPD': 'r',
    'TUB1_YPYE': 'm'}

cells = {}
for i in media:
    cells[i] = glob.glob(i+'*csv')

# linear regimes
L = {
 'ORI_YPD_252OD.csv': (196, 296),
 'ORI_YPD_450OD.csv': (200, 300),
 'ORI_YPD_499OD.csv': (168, 288),
 'ORI_YPD_544OD.csv': (140, 230),
 'ORI_YPD_578OD.csv': (210, 310),
 'ORI_YPD_870OD.csv': (148, 232),
 'ORI_YPD_874OD.csv': (147, 204),
 'ORI_YPYE_119OD.csv': (99, 175),
 'ORI_YPYE_170OD.csv': (150, 240),
 'ORI_YPYE_184OD.csv': (268, 355),
 'ORI_YPYE_210OD2.csv': (278, 364),
 'ORI_YPYE_300OD.csv': (445, 566),
 'ORI_YPYE_440OD.csv': (161, 227),
 'ORI_YPYE_562OD.csv': (296, 354),
 'TUB1_YPD_228OD.csv': (118, 188),
 'TUB1_YPD_268OD.csv': (165, 255),
 'TUB1_YPD_350OD.csv': (262, 375),
 'TUB1_YPD_399OD.csv': (205, 321),
 'TUB1_YPD_468OD.csv': (334, 456),
 'TUB1_YPD_478OD.csv': (199, 309),
 'TUB1_YPD_880OD.csv': (171, 247),
 'TUB1_YPD_947OD.csv': (137, 204),
 'TUB1_YPYE_190OD.csv': (215, 325),
 'TUB1_YPYE_200OD.csv': (168, 273),
 'TUB1_YPYE_320OD.csv': (254, 350),
 'TUB1_YPYE_336OD.csv': (261, 376),
 'TUB1_YPYE_451OD.csv': (125, 192),
 'TUB1_YPYE_520OD.csv': (142, 201)}

idx = 0   # counter
fig1 = plt.figure("MainFig", figsize=(8, 6))
gs = gridspec.GridSpec(2, 2)
plt.close('all')
mass = pd.read_csv('DryMass.csv')

# =============================================================================
#               Plot the O2 curve for each media type
# =============================================================================
for h, i in cells.iteritems():
    X.setdefault(h, [])
    Y.setdefault(h, [])
    plt.figure(h)

    for j in i:
        data = pd.read_csv(j, skiprows=26)
        x1 = L[j][0]
        x2 = L[j][1]

        time = data.Time[:-2].astype(float)
        O2 = data['Oxygen 1'][:-2].astype(float)
        ix1 = time[time == x1].index[0]
        ix2 = time[time == x2].index[0]

#   These are the actual O2 rate curves for each reading from Hansatech
        plt.plot(time[ix1:ix2],
                 O2[ix1:ix2])

        plt.text(time[ix2], O2[ix2], j)

        slope, intercept, r_value, p_value, std_err = (
            sp.linregress(time[ix1:ix2], O2[ix1:ix2]))

        Y[h].append(slope)
        X[h].append(float(j.rsplit('_', 1)[1][:3]))

        plt.show()

# =============================================================================
#    Plot the linear regression line from the O2 values
# =============================================================================
    plt.figure("MainFig")
    axes = plt.subplot(gs[idx])
    X[h] = np.array(X[h])
    Y[h] = np.array(Y[h])

    slope2, intercept2, r_value2, p_value2, std_err2 = (
        sp.linregress(X[h], Y[h]))

    line = slope2*X[h]+intercept2
    axes.plot(X[h], line, colors[h], X[h], Y[h], 'o', mfc=colors[h])
    plt.ylim(-1.0, -0.0)

    plt.title(j.rsplit('_', 1)[0])
    plt.text(250, -.5, "Slope= %7.5f"
             % slope2)
    plt.text(250, -.6, "Intercept= %5.3f"
             % intercept2)
    plt.text(250, -.7, "R_sq= %5.3f"
             % r_value2**2)
    plt.text(250, -.8, "O2@OD500= %5.3f"
             % (slope2*500+intercept2))

    O2_500[h] = slope2*500+intercept2
    idx += 1
    plt.show()

    W[h] = mass.Weight[
        mass.Label.str.contains(j.rsplit('_', 1)[0])]
    OD500[h] = mass.OD[
        mass.Label.str.contains(j.rsplit('_', 1)[0])]
    CellCount[h] = mass.cellCount[
        mass.Label.str.contains(j.rsplit('_', 1)[0])]
    OD2[h] = mass.OD2[
        mass.Label.str.contains(j.rsplit('_', 1)[0])]

    MassperMil[h] = [W[h]/40/OD500[h]*1000*.5]

    O2_500perMil[h] = [O2_500[h]/OD500[h]*.5]  # @OD0.5

    CountperOD[h] = [CellCount[h]/OD2[h]/2]  # @OD 0.5

    O2Rate[h] = [O2_500[h]/np.mean(MassperMil[h][0])]

plt.figure("DryWeight")
plt.boxplot(
    [MassperMil[k][0] for k in MassperMil],
    labels=media)
plt.ylim(0.0, 1.0)
plt.title('Dry Weight (mg/ml @ 0.5 OD)')
plt.show()

plt.figure("CellNumber")
plt.boxplot(
    [CountperOD[k][0] for k in CountperOD],
    labels=media)
plt.title('Number of Cells (10${^7}$/ml @ OD0.5)')
plt.show()
