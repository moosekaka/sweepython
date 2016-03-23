from mpl_toolkits.mplot3d import Axes3D
import glob,string,os
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
from tvtk.api import tvtk

parDir=os.path.dirname(os.getcwd())

r=pd.read_csv(parDir+'\\'+'BackgroundRFPGFP.csv')
g=pd.read_csv(parDir+'\\'+'BackgroundRFPGFP.csv')
backgroundGFP={i[0]:i[1] for i in zip(g['Label'],g['MeanGFP'])}
backgroundRFP={i[0]:i[1] for i in zip(r['Label'],r['MeanRFP'])}
with open(parDir+'\\fileMetas.pkl','rb') as inpt:
    s=pickle.load(inpt)

b=glob.glob(parDir+'\[0-9][0-9]*') # directory names
labs=[i.rsplit('\\',1)[1].split('_',1)[0][2:] for i in b]
meanFP={};dates=[];powersRFP=[];powersGFP=[]
bckgrndGFP=[];bckgrndRFP=[];

colors = ['orange','b','g','r','m','c']
colors={'YPD_':'b','YPE_':'g','YPR_':'r','YPL_':'m','YPL-':'c','A1RA':'orange'}
labs={'YPD_':'YPD','YPE_':'YPE','YPR_':'YPR','YPL_':'YPL','YPL-':'atg32','A1RA':'A1RA1G'}

fileLoc='*2.5_raw*vtk'    
files = glob.glob(fileLoc)      

for a in range(len(files)):

    fileKey=files[a][8:-13]
    reader=tvtk.PolyDataReader()
    reader.set(file_name=files[a])
    reader.update()
    data=reader.output
    temp=data.point_data
    rawGFP=np.ravel(temp.get_array('rGFP'))
    rawRFP=np.ravel(temp.get_array('rRFP'))
    
    
    fileKey=files[a][8:-13]
    dates.append(s[fileKey][0])
    powersRFP.append(s[fileKey][1])
    powersGFP.append(s[fileKey][2])
    meanFP.setdefault(fileKey[:4],[])
    
    
    temp=(np.mean(rawRFP),np.mean(rawGFP),fileKey,powersRFP[a],powersGFP[a])
    meanFP[fileKey[:4]].append(temp)


#fig = plt.figure()
#ax = fig.gca(projection='3d')
#condition='YPD_'
#x=[j[0] for j in meanFP[condition]]#inten
#z=[float(j[3]) for j in meanFP[condition]]
#
#y=[j[1] for j in meanFP[condition]]
#t=[j[2] for j in meanFP[condition]]
#ax.scatter(x, y, z)
#ax.set_ylabel('GFP')
#ax.set_xlabel('RFP')
#ax.set_title(condition)
#plt.show()

plt.close('all')
gs = gridspec.GridSpec(5, 2)
fig=plt.figure(figsize=(11,17))
for h,i in enumerate(sorted(colors.keys())[1:]):
    
    ax1=plt.subplot(gs[2*h])    
    ax2=plt.subplot(gs[2*h+1])  
    
    condition=i;z=[];zGFP=[]
    
    x=[j[0] for j in meanFP[condition]]
    for j in range(len(x)):
        try:
            z.append(float(meanFP[condition][j][3]))
            zGFP.append(float(meanFP[condition][j][4]))
        except:
            IndentationError
            z.append(0)
            zGFP.append(0)
    y=[j[1] for j in meanFP[condition]]
    t=[j[2] for j in meanFP[condition]]

    ax1.plot(z,x,'sr')
    ax2.plot(z,y,'sg')
    
    for k in range(len(t)): 
        if z[k]>1.7 and x[k]>5000:
            ax1.annotate(t[k].rsplit('_',1)[1],xy=(z[k],x[k]))

            ax2.annotate(t[k].rsplit('_',1)[1],xy=(z[k],y[k]))

    ax1.set_xlim(1.1,2.4)
    ax2.set_xlim(1.1,2.4)

    if h==0:
        ax1.set_title('RFP mean for '+condition)
        ax2.set_title('GFP mean for '+condition)
    ax1.set_xlabel('laser power RFP')
    ax2.set_xlabel('laser power RFP')
    ax1.text(0.85,0.5,labs[condition],transform=ax1.transAxes, fontsize=14)
    ax2.text(0.85,0.5,labs[condition],transform=ax2.transAxes, fontsize=14)

plt.show()

    
  

    #for i in range(len(t)): 
    #    if z[i]!='NA' and z[i]>2.1:
    #        plt.annotate(t[i],xy=(z[i],x[i]))
