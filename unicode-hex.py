# -*- coding: utf-8 -*-
"""
Created on Thu Apr 07 16:33:00 2016
@author: sweel_rafelski
Quick module to show unicode codepoints for greek and coptic and also proper
way to escape/decode the strings when using str.format() method
"""


hexrange = [i for i in range(0x00, 0xFF+1)]  # hex range are iterable

for h in hexrange:
    t = "{:02X}".format(h)          # upper case hex, pad leading zero 2 digits
    a = '\u03' + t                  # greek starts from 03XX
    b = a.decode('unicode_escape')  # decode \u escape to cast b to unicode
    print u"Unicode CodePoint : {} = {}".format(a, b)
