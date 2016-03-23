import os,glob
from itertools import cycle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib.patches as mpatches

colors = ['r','b','grey','orange','y']
cycler =cycle(colors)

data=pd.read_csv(os.getcwd()+'\\'+'Results.csv')
A=data['mRuby2'];col=cycler.next()
B=data['mKO2'];col=cycler.next()
C=data['dsRed']
patch = mpatches.Patch(color=col, label=['mRuby2','mKO2','dsRed'],alpha=0.4)

coefficients = np.polyfit(range(15), A[:15],4)
polynomial = np.poly1d(coefficients)
xs = np.arange(0,15, 0.5)
ys = polynomial(xs)

coefficients = np.polyfit(range(19), B[:19],2)
polynomial = np.poly1d(coefficients)
xs1 = np.arange(0,19, 0.5)
ys1 = polynomial(xs1)


coefficients = np.polyfit(range(21), C,2)
polynomial = np.poly1d(coefficients)
xs2 = np.arange(0,21, 0.5)
ys2 = polynomial(xs2)


plt.figure()
plt.plot(xs,ys,'r',xs1,ys1,'g',xs2,ys2,'b',A,'rs',B,'gs',C,'bs')

#plt.plot(B,'gs')
plt.title('mRuby2 (2.09), mKO2 (1.05),dsRed (1.45) intensity vs timepoints',fontsize=12)

plt.show()
