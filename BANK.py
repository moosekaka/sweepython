# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 12:29:53 2015

@author: sweel
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from seaborn import xkcd_palette as scolor
# pylint: disable=C0103


def rgbcol(colorname):
    """
    returns RGB float values based on xkcd color name input string
    """
    return scolor([colorname])[0]

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
pattern = 'Payment.*Chase|^CHASE|ALLY|PAYMENT.*CHASE|\
BILL.*LATER|UC IRVINE|SO CAL EDISON|CARDMEMBER'
debit = debit[~debit.Desc.str.contains(pattern)]
spend = pd.concat([sales, returns, debit]).reset_index(drop=True)
spend = spend.drop('Type', axis=1)
df1 = pd.pivot_table(spend,
                     values='Amount',
                     index='Date',
                     aggfunc='sum')
MthlySpd = df1.resample('M').sum()

# =============================================================================
#       item breakdowns
# =============================================================================
spend['Category'] = ''

spend.loc[spend.Desc.str.contains('ATM'), ['Category']] = 'ATMcash'

pattern = r'COF|PEE|STARBU|ARAMARK UCI|BED'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = 'Coffee'

pattern = 'JACK|^CHICK|MCD|DEL TACO|IN-N-OUT|DAPHNE|RAISING\
|TACO BELL|SUBWAY|PANERA|EL POLLO LOCO|CHIPOTLE|BLAZE|HABIT|IKEA R\
|SEAFOOD C|DENNY|IHOP|POPEY|SLAPFISH|Flame|FLAME|KFC|Pizza 90\
|CHRONIC|PHO.*1|BRUXIE|MENDOCINO|TLT FOOD|YOGURT|MEET FRESH|MITSUWA-BEER'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = 'FastFood'

pattern = 'J J BAKE|NORM|TENDER GR|URBAN SEOUL|KAISEN|KULA|BLAKES\
|NANA SAN|TAPS|ORIGINAL FISH|DODDY|FOGO|RUTH|URBAN|BAJA|Baja|CA FISH\
|WOKCANO|BLACK BEAR DINER|CURRY.*AMERICA|I HEART PANCAKES|LUNA GRILL\
|DICKIE\'S|PHO.*PEARL|BCD|BOIL|Boil'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = 'Restaurants'

spend.loc[spend.Desc.str.contains('LAZY'), ['Category']] = 'Lazy'

spend.loc[spend.Desc.str.contains('USA|A[Rr][Cc]|76 - UNION\
|UNITED OIL'), ['Category']] = 'Gas'

pattern = '99 RAN|PAVILIONS|VONS|ALBERT\
|RALP|TARG|^WAL|WM|^TRADER|MITSUWA MAR|WHOLESOME|H MART|SPROUTS\
|MITSUWA MRKTPLACE|WHOLEFDS'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = 'Market'

pattern = 'Kindle|MICRO CENTER|FRY\'S\
|DSW METRO POINT|IRVINE-BOOKSTORE-HILL|T-MOBI|SAMY'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = '_Books, Electro.'

pattern = 'UCI CAMPUS RECREATION|AMC|FANDANGO|CINEMA|EDWARD'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = '_Gym, movies'

pattern = 'FARMERS IN|AM SOC|ASCB|U-HAUL|CA DMV\
|UCI TRANSPO|VISTA DEL|IMAGE D|COINBASE|EVA|UMI|PRIME TIME|VR|HARBOR J \
|BIOMEDICAL ENGINEERING|AUTOZONE|PERFORMANCE BIKE SHOP\
|STATE OF CALIF DMV|OPTOMETRIC|TOWING|AAA|CAPGOWN\
|H AND J AUTO REPAIR|COMFORT INN|CASA LOMA|PALA CHOICES\
|BIOPSYCHOSOCIAL|ENTERPRISE RENT|INSURE.*RENTAL|SEQUOIA KINGS|VILLAGE MARKET\
|PRINCE FOOD|DNPSKINGSCANYONGROCERY|SEQUOIA- LODGEPOLE'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = '_One-offs, Fees'
spend.loc[spend.Desc.str.contains('IKEA') &
         (spend.Amount < -12), ['Category']] = '_One-offs, Fees'
spend.loc[spend.Desc.str.contains('TARGET') &
         (spend.Date=='2016-07-08'), ['Category']] = '_One-offs, Fees'

pattern = 'UCI PARK|GOOGLE.*PRO|UCI TRANS'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = 'Phone, Parking'

pattern = 'AMZ  STORECARD|PAYPAL|AMAZON'
spend.loc[spend.Desc.str.contains(pattern), ['Category']] = '_Installments'


#   find cashback items, exclude cashback items from original set
grocash = spend.loc[spend.Desc.str.contains('Cash Back')]
spend = spend[~spend.index.isin(grocash.index)]
#   create two separate groups of actual cashback and actual grocery
values = grocash.Desc.str.rsplit(pat=None, n=1).str[1]
converted = [-1 * float(i) for i in values.str.strip('$')]

cashback = pd.DataFrame(
    {'Date': grocash.Date,
     'Desc': grocash.Desc,
     'Amount': converted,
     'Category': 'ATMcash'},
    index=values.index)

netgro = pd.DataFrame(  # grocery items minus cashback
    {'Date': grocash.Date,
     'Desc': grocash.Desc,
     'Amount': grocash.Amount - cashback.Amount,  # net to take out cashback
     'Category': 'Market'},
    index=values.index)

#    concat back the excluded , new groups
total = pd.concat([spend, cashback, netgro], ignore_index=True)
uncat = total[total.Category.apply(lambda x: x == '')]
total = total.set_value(uncat.index, 'Category', 'Uncat')

# =============================================================================
#    Pivot to group by cat and months and plot
# =============================================================================
date_range = '2016-01'
dollar_lim = 2200
mth_lim = 350
MthlySpd = - MthlySpd.ix[date_range:]
df = pd.pivot_table(total,
                    index='Date',
                    columns='Category',
                    values='Amount',
                    aggfunc='sum')
byMonth = df.resample('M').sum()
cols = sorted([c for c in pd.unique(total.Category.values) if c !='Uncat'])
Filt = -byMonth.ix[date_range:, :]

# =============================================================================
#     Plot
# =============================================================================
f1, ax = plt.subplots()
ax.set_ylim(0, MthlySpd.max() + dollar_lim)
ax.set_axis_bgcolor(rgbcol('grey'))
ax.set_title('Stacked Monthly Expenditure')
Filt.ix[:, cols].plot(kind='area',
                      colormap=plt.cm.Paired,
                      ax=ax)
MthlySpd.plot(ax=ax)
ax.set_ylim(0, dollar_lim)
ax.axhline(1000, c=rgbcol('hot pink'))

f2, ax2 = plt.subplots()
ax2 = Filt.ix[:, cols].plot(kind='bar',
                            colormap=plt.cm.Paired,
                            ax=ax2)
ax2.set_axis_bgcolor(rgbcol('grey'))
ax2.set_xticklabels(
    [dt.strftime('%b %y')
     for dt in Filt.index.to_pydatetime()])
ax2.set_ylim(0, mth_lim)
plt.xticks(rotation=0)

# =============================================================================
# Pivot table
# =============================================================================
x = total[(total.Date > '2016-01-31')].reset_index(drop=True)
y = x.groupby(['Date', 'Category']).sum()
y1 = y.reset_index()
ycount = x.groupby(['Date', 'Category']).count()
#
z = y.unstack(level=1).resample('M').sum()
z_count = y.unstack(level=1).resample('M').count()
z.columns = z.columns.droplevel(0)
z_count.columns = z_count.columns.droplevel(0)
z['Total'] = z.sum(axis=1)

# ============================================================================
# Number of visits to xxx
# ============================================================================
#print '\n{}\nNumber of occur. of \n{}\n{}\n'.format('='*79, '='*79, z_count)
print z_count
# ============================================================================
# Category spending
# ============================================================================
#print '\n{}\nSpending by category to \n{}\n{}\n'.format('='*79, '='*79, z)
print z
# =============================================================================
# Write out to file
# =============================================================================
with open('total.txt', 'w') as output:
    output.write('\n{}\nNumber of occur.\
    of \n{}\n{}\n'.format('='*79, '='*79, z_count))
    output.write('\n{}\nSpending by category\
    of \n{}\n{}\n'.format('='*79, '='*79, z))
