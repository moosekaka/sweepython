# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 12:29:53 2015

@author: sweel
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use('ggplot')
plt.close('all')

# =============================================================================
#   chase credit card
# =============================================================================
dataF = pd.read_csv('ChaseCredit.csv',
                    parse_dates=[1],
                    header=0,
                    names=['Type', 'Date', 'Desc', 'Amount'],
                    usecols=[0, 1, 3, 4])
sales = dataF.loc[dataF.Type == 'SALE']
returns = dataF.loc[dataF.Type == 'RETURN']

# =============================================================================
#    chase checking acct
# =============================================================================
dataJ = pd.read_csv('Checking.CSV',
                    parse_dates=[1],
                    header=0,
                    names=['Type', 'Date', 'Desc', 'Amount'],
                    usecols=[0, 1, 2, 3])
debit = dataJ[dataJ.Type == 'DEBIT']

#   exclude credit card purchases,bill me later/paypal and rent/utilities
pattern = 'Payment.*Chase|^CHASE|ALLY|PAYPAL|\
BILL.*LATER|UC IRVINE|SO CAL EDISON|CARDMEMBER'
debit = debit[~debit.Desc.str.contains(pattern)]
spend = pd.concat([sales,returns, debit]).reset_index(drop=True)
spend = spend.drop('Type', axis=1)
df1 = pd.pivot_table(spend,
                     values='Amount',
                     index='Date',
                     aggfunc='sum')
MthlySpd = df1.resample('M', how='sum')
MthlySpd = - MthlySpd.ix['2015-2':]

# =============================================================================
#       item breakdowns
# =============================================================================
spend['Category'] = ''
spend.loc[spend.Desc.str.contains('ATM'), ['Category']] = 'ATMcash'
pattern = r'COF|PEE|STAR|ARAMARK UCI'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = 'Coffee'
spend.loc[spend.Desc.str.contains('BCD|BOIL|Boil'), ['Category']] = 'BcdBoil'
pattern = 'JACK|^CHICK|MCD|^DEL|^IN\
|TACO BELL|SUBWAY|PANERA|^EL P|CHIPOTLE|BLAZE|HABIT|IKEA R|SEAFOOD C|DENNY|IHOP'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = 'FastFood'
pattern = 'J J BAKE|NORM|TENDER GR|UBRAN SEOUL|KAISEN|KULA|BLAKES\
|NANA SAN|TAPS|ORIGINAL FISH|DODDY|FOGO'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = 'Restaurants'
spend.loc[spend.Desc.str.contains('FLAME'), ['Category']] = 'FlameB'
spend.loc[spend.Desc.str.contains('BAJA|Baja'), ['Category']] = 'BajaF'
spend.loc[spend.Desc.str.contains('LAZY'), ['Category']] = 'Lazy'
pattern = 'YOGURT|MEET FRESH|BEER'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = 'DrinksDesert'
spend.loc[spend.Desc.str.contains('USA|A[Rr][Cc]'), ['Category']] = 'Petrol'
pattern = '99 RAN|PAVILIONS|VONS|ALBERT\
|RALP|TARG|^WAL|WM|^TRADER|MITSUWA MAR|WHOLESOME|H MART|SPROUTS|BED'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = 'Market'
pattern = 'AMC|FANDANGO|CINEMA|EDWARD|Kindle'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = 'Entertainment'
pattern = 'FARMERS IN|AM SOC|USPS 053711025|U-HAUL\
|UCI TRANSPO|VISTA DEL|IMAGE D|COINBASE|EVA|UMI|PRIME TIME'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = 'OneTime'

#   find cashback items, exclude cashback items from original set
grocash = spend.loc[spend.Desc.str.contains('Cash Back')]
spend = spend[~spend.index.isin(grocash.index)]
#   create two separate groups of actual cashback and actual grocery
values = grocash.Desc.str.rsplit(pat=None, n=1).str[1]
converted = [-1*float(i) for i in values.str.strip('$')]

cashback = pd.DataFrame(
    {'Date': grocash.Date,
     'Desc': grocash.Desc,
     'Amount': converted,
     'Category': 'ATMcash'},
    index=values.index)

netgro = pd.DataFrame(  # grocery items minus cashback
    {'Date': grocash.Date,
     'Desc': grocash.Desc,
     'Amount': grocash.Amount-cashback.Amount,  # net to take out cashback
     'Category': 'Market'},
    index=values.index)
#    concat back the excluded , new groups
total = pd.concat([spend, cashback, netgro], ignore_index=True)

# =============================================================================
#    Pivot to group by cat and months and plot
# =============================================================================
df = pd.pivot_table(total,
                    index='Date',
                    columns='Category',
                    values='Amount',
                    aggfunc='sum')
byMonth = df.resample('M', how='sum')
cols=[u'ATMcash', u'BajaF', u'BcdBoil', u'Coffee', u'DrinksDesert',
       u'Entertainment', u'FastFood', u'FlameB', u'Lazy', u'Market',
      u'Petrol', u'Restaurants']
Filt = -byMonth.ix['2015-2':, 1:]

#    Plot
fig = plt.plot(figsize=(11, 8.5))
ax = plt.subplot(111)
ax.set_ylim(0, MthlySpd.max()+200)
ax.set_axis_bgcolor('#C0C0C0')
ax.set_title('Stacked Monthly Expenditure')
Filt.ix[:,cols].plot(kind='Area',
          colormap=plt.cm.Paired,
          sort_columns=True,
          ax=ax)
MthlySpd.plot(ax=ax)

fig2 = plt.plot(figsize=(11, 8.5))
ax2 = plt.subplot(111)
ax2 = Filt.ix[:,cols].plot(kind='bar', colormap=plt.cm.Paired)
ax2.set_axis_bgcolor('#C0C0C0')
ax2.set_xticklabels(
    [dt.strftime('%b %y')
     for dt in Filt.index.to_pydatetime()])
plt.xticks(rotation=0)

#
x=total[(total.Date>'2015-05-31')]
y=x.groupby(['Date', 'Category']).sum()
ycount=x.groupby(['Date', 'Category']).count()
#
z=y.unstack(level=1).resample('M',how='sum')
zcount= ycount['Amount'].unstack(level=1).resample('M',how=sum)
ztot=z.sum(axis=1)
yf=y.query('Category == "Lazy"')
yf1=yf.reset_index(level=1)
print yf1.resample('M', how=len).Category
print z
print '%s\n\n\n%s\n' % (z, zcount)
uncat = total[(total.Date>'2015-09-30') &(total.Category =='')]
#uncat = total[(total.Date>'2015-10-30') &(total.Category =='ATMcash')]