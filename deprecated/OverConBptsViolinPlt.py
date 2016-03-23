#Plot violin plots for overconnected Branchpoints

import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
import pandas as pd
import matplotlib.gridspec as gridspec
from statsmodels.graphics import boxplots

colors = ['b','g','r','m','c']
cycler =cycle(colors)
labs=['YPD','YPE','YPR','YPL','YPLatg32'] 
R=pd.read_csv('rejectlst.csv')

C=pd.read_csv('checklistCSV.csv',na_values=0)
Cfilt=C[~C['Cell Name'].isin(R['Label'])]
CellName=np.array(Cfilt['Cell Name'])
xBpts=np.array(Cfilt['OverConnected Node'].fillna(0))
Bpts=np.array(Cfilt['#bpts'])
#data=xBpts/Bpts
data=Bpts-xBpts

YPD=[i for h,i in enumerate(data) if CellName[h][:3]=='YPD' ]
YPE=[i for h,i in enumerate(data) if CellName[h][:3]=='YPE' ]
YPR=[i for h,i in enumerate(data) if CellName[h][:3]=='YPR' ]
YPL=[i for h,i in enumerate(data) if CellName[h][:4]=='YPL_' ]
YPLatg32=[i for h,i in enumerate(data) if CellName[h][:4]=='YPL-' ]
#A1RA1G=[i for h,i in enumerate(data) if CellName[h][:3]=='A1R' ]
N=[len(i) for i in [YPD,YPE,YPR,YPL,YPLatg32]]


plt.close('all')
fig=plt.figure(figsize=(11,8.25))
ax = fig.add_subplot(111)
for i in range(6):
    boxplots.violinplot([YPD,YPE,YPR,YPL,YPLatg32],ax=ax,labels=labs)
    
height=.92
ax.set_ylim(0)
#ax.set_title('Percentage of Overconnected Branchpoints')
ax.set_title('Average Number of Branchpoints per Cell')
#ax.text(0.1, height, 'N='+str(N[0]), transform=ax.transAxes, fontsize=14,horizontalalignment='left')
#ax.text(0.22 ,height, 'N='+str(N[1]), transform=ax.transAxes, fontsize=14,horizontalalignment='left')
#ax.text(0.38, height, 'N='+str(N[2]), transform=ax.transAxes, fontsize=14,horizontalalignment='left')
#ax.text(0.53, height, 'N='+str(N[3]), transform=ax.transAxes, fontsize=14,horizontalalignment='left')
#ax.text(0.72, height, 'N='+str(N[4]), transform=ax.transAxes, fontsize=14,horizontalalignment='left')
#ax.text(0.9, height, 'N='+str(N[5]), transform=ax.transAxes, fontsize=14,horizontalalignment='left')



#for patch, color in zip(bplot['boxes'], colors[:4]):
#    patch.set_facecolor(color)
#    patch.set_alpha(0.7)
#    patch.set_linewidth(1)


plt.show()

