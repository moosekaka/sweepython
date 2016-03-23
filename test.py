# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 15:27:07 2016

@author: sweel_rafelski
"""

#Macro to save FP and BF  files into a stack
from ij import IJ
from ij.gui import Roi, ShapeRoi
import os
import os.path as op
import glob
import string
import shutil as sh
import fnmatch as fn
import errno


#Power=getTag("DAC2_561-Volts");
#Gain=getTag("Hamamatsu_DCAM-EMGain");
#AcqTime=getTag("Exposure-ms");
#
#print(Power);
#print(Gain);
#print(AcqTime);

"""function getTag(tag) {
      info = getImageInfo();
      index1 = indexOf(info, tag);
      if (index1==-1) return "";
      index1 = indexOf(info, ":", index1);
      if (index1==-1) return "";
      index2 = indexOf(info, "\n", index1);
      value = substring(info, index1+1, index2);
      return value;
  }
"""
dr = IJ.getDirectory("Choose a Directory ");
for root, dirs, files in os.walk(dr):
    for f in files:
    	if fn.fnmatch(f, '*RFP*'):
        	print op.join(root, f)