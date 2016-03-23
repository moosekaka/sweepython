# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 15:23:43 2015

@author: sweel
"""

import pandas as pd

datelisttemp = pd.date_range('1/1/2014', periods=3, freq='D')
s = list(datelisttemp)*3
s.sort()
df = pd.DataFrame({'BORDER':['GERMANY','FRANCE','ITALY',
                             'GERMANY','FRANCE','ITALY',
                             'GERMANY','FRANCE','ITALY' ],
                   'HOUR1':[2 ,2 ,2 ,4 ,4 ,4 ,6 ,6, 6],
                   'HOUR2':[3 ,3 ,3, 5 ,5 ,5, 7, 7, 7],
                   'HOUR3':[8 ,8 ,8, 12 ,12 ,12, 99, 99, 99]},
                    index=s)

df = df.set_index(['BORDER'], append=True)
df.columns.name = 'HOUR'
df = df.unstack('BORDER')
df = df.stack('HOUR')
df = df.reset_index('HOUR')
df['HOUR'] = df['HOUR'].str.replace('HOUR', '').astype('int')
print(df)



df = pd.DataFrame({'BORDER':['GERMANY','FRANCE','ITALY',
                             'USA','CANADA','MEXICO',
                             'INDIA','CHINA','JAPAN' ],
                   'ASID':[21, 32, 99, 77,66,55,44,88,111],
                    'HOUR1':[2 ,2 ,2 ,4 ,4 ,4 ,6 ,6, 6],
                             'HOUR2':[3 ,3 ,3, 5 ,5 ,5, 7, 7, 7],
                             'HOUR3':[8 ,8 ,8, 12 ,12 ,12, 99, 99, 99],
                             'PRICE1':[2 ,2 ,2 ,4 ,4 ,4 ,6 ,6, 6],
                             'PRICE2':[2 ,2 ,2 ,4 ,4 ,4 ,6 ,6, 6],
                             'PRICE3':[2 ,2 ,2 ,4 ,4 ,4 ,6 ,6, 6] })

df = df[['ASID', 'BORDER', 'HOUR1', 'PRICE1', 'HOUR2', 'PRICE2', 'HOUR3', 'PRICE3']]

df.set_index(['ASID', 'BORDER'], inplace=True)
prices = df[['PRICE1','PRICE2', 'PRICE3']].stack()
print(prices)
prices.index = prices.index.droplevel(2)
print(prices)
hours = df[['HOUR1', 'HOUR2', 'HOUR3']].stack()
hours.index = hours.index.droplevel(2)

df_new = pd.concat([hours, prices], axis=1)
df_new.columns = ['HOUR', 'PRICE']