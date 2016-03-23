#######  visualize the skel and surface vtk file ########################
#########################################################################
import matplotlib.pyplot as plt
import math,glob,string,os
import numpy as np
from skimage import io
import tifffile as tif
parDir=os.path.dirname(os.getcwd())

fileMatch=[\
'A1RA1G_038_RFPstack',\
'A1RA1G_044_RFPstack',\
'A1RA1G_046_RFPstack',\
'A1RA1G_051_RFPstack',\
'A1RA1G_053_RFPstack',\
'A1RA1G_055_RFPstack',\
'A1RA1G_061_RFPstack',\
'A1RA1G_067_RFPstack',\
'A1RA1G_071_RFPstack']

plt.close('all')


MaxProjRFP=tif.imread(os.getcwd()+'\RFP\MaxProjs.tif',multifile=1)
MaxProjGFP=tif.imread(os.getcwd()+'\GFP\MaxProjs.tif',multifile=1)

files=glob.glob('RFP\*tif')
Dict={}

for h,i in enumerate(files[:-1]):
    fileKey=i.rsplit('\\',1)[1][:-4]
    Dict[fileKey]=h

for k in fileMatch:
    image1=MaxProjRFP[Dict[k]]
    image2=MaxProjGFP[Dict[k]]
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(11,8.25))
    fig.canvas.set_window_title(k)
    ax1.imshow(image1,cmap='gray')
    ax2.imshow(image2,cmap='gray')
   
    ax1.text(0.5,0.9,'RFP',transform=ax1.transAxes,color='r')
    ax2.text(0.5,0.9,'GFP',transform=ax2.transAxes,color='g')
    fig.text(0.5,0.85,k,verticalalignment='top',\
    horizontalalignment='center')
    plt.savefig(parDir+'\\'+k+'.png')
    
    
      