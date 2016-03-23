import glob
import matplotlib.pyplot as plt
from itertools import cycle
import pandas as pd
import matplotlib.patches as mpatches
import seaborn as sns
sns.set_context("talk")
sns.set(style="darkgrid")

colors = ['r', 'b', 'grey', 'orange', 'y']
cycler = cycle(colors)


def pltSpectra(string, dic):
    """Function to plot emission spectra and return name of FP, legend handle

    Parameters
    ----------
    string :
        name of file

    dic :
        dictionary of ascii text files

    Returns
    -------
    string :
        name of file
    patch :
        patch data handle for making legend
    """
    data = pd.read_table(dic[string], header=None)
    fig = plt.gcf()
    ax = fig.gca()
    col = cycler.next()
    ax.fill_between(data[0], data[1], color=col, alpha=0.5)
    patch = mpatches.Patch(color=col, label=string, alpha=0.4)
    plt.show()
    return (string, patch)

plt.close('all')
files = glob.glob('*ascii.txt')

#   Plot z-quad filter bandwitdh
TRFilt = pd.read_table('TRFilt-ascii.txt', header=None)
fig = plt.figure(figsize=(11, 8.25))
ax = plt.subplot(1, 1, 1)
ax.plot(TRFilt[0], TRFilt[1], 'k', lw=.75)

#   List of available spectra data
Dict = {name.rsplit('-', 1)[0]: name for name in files}
count = 0
for key in sorted(Dict.keys()):
    count += 1
    print('%2g : %s' % (count, key))

#   Plot FP emission spectra
temp = []
temp.append(pltSpectra('dsRed-em', Dict))
temp.append(pltSpectra('dsRed-ex', Dict))

plt.legend(
    handles=[i[1] for i in temp],
    loc='right',
    fontsize='small')
plt.axvline(488, color='#00f7ff')
plt.axvline(561, color='#b3ff00')
plt.axvline(640, color='#ff2100')
plt.axvline(405, color='#8200c8')
ax.set_xlabel('Wavelength/nm')
ax.set_ylabel('Relative Intensity/%')
ax.set_title('Emission Spectra for %s' % (temp[0][0][:-3]))
