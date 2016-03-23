#To be able to read csv formated files, we will first have to import the
#csv module.
import os,glob,string
import numpy as np
import matplotlib.pyplot as plt
import csv
from itertools import cycle
import scipy.stats as sp

fileLoc=os.getcwd()+'\\'+'a*csv'
files = glob.glob(fileLoc)
lines = ["--","-.",":"]
linecycler =cycle(lines)
Y=[];X=[]
def O2proc(fileNum):

#iterate thru the filelist, read csv data and store in list
    lists=[]
    with open(fileNum, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            lists.append(row)

#find the start and stop time and store the row num in lists3
    #lists2=['Start','Stop']
    lists2=['T0','T1']
    lists3=[]
    for h,i in enumerate(lists):
        for j in i:
            if j in lists2:
                lists3.append(h)

#actual time and O2 reading data stored  here
    data=lists[lists3[0]:lists3[1]]
    time=[float(i[0]) for i in data]
    O2=[float(i[1]) for i in data]
    slope, intercept, r_value, p_value, std_err = sp.linregress(time,O2)


#plt data
    l1=fileNum.split('\\')

    plt.plot(time,O2,label=l1[-1],linestyle=linecycler.next());
    plt.ylim(0,380)
    plt.show()
    temp=l1[-1] #file name desc media type and OD
    Y.append(slope) # respirate rate slope
    X.append(temp[1+temp.find('D'):temp.find('.')])

    return(slope,r_value,p_value,X,Y,temp)

plt.figure()
for a in range(len(files)):
    #plt.figure()
    O2proc(files[a])
    plt.legend(loc=0)

plt.show()

X=[float(i) for i in X]
X=np.array(X)
Y=np.array(Y)
slope, intercept, r_value, p_value, std_err = sp.linregress(X,Y)
line = slope*X+intercept
plt.figure()
plt.show()
plt.plot(X,line,'m-',X,Y,'o')
plt.ylim(-1.0,-0.0)
