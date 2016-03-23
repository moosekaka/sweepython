# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 12:42:25 2015

@author: sweel
"""
import matplotlib.pyplot as plt
sns.set_context("talk")
sns.set(style="white")
plt.close('all')
plt.figure(figsize=(11,8.5))

annotation_string = r'$The\ vector\ F(k) \Rightarrow\ lagged\ intensity\ gradients:$'
annotation_string += '\n\n'
annotation_string += r'$F(k) \to \{ \frac{I_{n+k} - I_{n}}{k} \} \in \lbrack n = k .... N-k \rbrack$'
annotation_string += '\n\n'
annotation_string += r'$where\ N = \ number \ of \ points$'
annotation_string += '\n'
annotation_string += r'$I_{n} = intensity\ at\ point\ N$'
annotation_string += '\n'
annotation_string += r'$k = lag\ value$'

plt.annotate(annotation_string, xy=(0.15, 0.4), ha='left', fontsize=14)
plt.savefig('Texfile.png')
plt.show()
