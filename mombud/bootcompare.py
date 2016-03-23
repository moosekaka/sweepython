# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 01:24:48 2016

@author: sweel_rafelski
"""

def bootO2(X, N):
    '''boostrap the CI for the O2slope at OD=0.5
    '''
    result = []
    result2 = []
    result3 = []
    result4 = []
    result5 = []
    lenB = len(X)
#    lenB1 = len(Y)
    for _ in range(N):
        rX = X.sample(n=lenB, replace=True)
#        rY = Y.sample(n=lenB1, replace=True)
        A, B, _, _, _ = sp.linregress(rX.ix[:, r'$OD_{600}$'],
                                      rX.ix[:, 'O2slope'])
        ocr = A * .5 + B
        idx = X.index.tolist()
        cell_num = X.ix[idx[0], 'count'] / 2
        cellmass = X.ix[idx[0], 'mass'] / 2
#        cell_num = rY['countOD1'].median() / 2
        # 1E-6 * (1E-6)^2 gives micron^3, divide by 1E18 to be in picom^3
        mitovol = X.ix[idx[0], 'mitovol'] * (.15**2) * math.pi
        ocrm = ocr / cell_num / mitovol / 10  # factor to make consistent
        ocrc = ocr / cell_num
        ocrw = ocr / cellmass
        result.append(A * .5 + B)
        result2.append(ocrm)
        result3.append(cell_num)
        result4.append(ocrc)
        result5.append(ocrw)
    return result, result2, result3, result4, result5