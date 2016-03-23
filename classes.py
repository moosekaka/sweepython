# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 17:48:48 2016

@author: sweel
"""

def gcd(m,n):
    while m%n != 0:
        oldm = m
        oldn = n

        m = oldn
        n = oldm%oldn
    return n

class Fraction:

    def __init__(self,top,bot):
        if (top % 1 ==0) and (bot % 1==0):
            com = gcd(top, bot)

            self.num=top/com
            self.den=bot/com

        else:
            raise RuntimeError ("must be integers")

    def __str__(self):
        return '%s/%s' % (self.num, self.den)


    def __add__(self,other):
        num1 = (self.num * other.den)+(other.num * self.den)
        den1 = self.den * other.den
        tot = Fraction(num1, den1)

        return tot

    def __sub__(self,other):
        num1 = (self.num * other.den)-(other.num*self.den)
        den1 = self.den*other.den
        tot = Fraction(num1, den1)


        return tot

    def __mul__(self,other):
        num1 = self.num * other.num
        den1 = self.den * other.den
        tot = Fraction(num1, den1)

        return tot


    def __div__(self,other):
        num1 = self.num * other.den
        den1 = self.den * other.num
        tot = Fraction(num1, den1)

        return tot

    def __lt__(self,other):
        com = gcd(self.den,other.den)
        return self.num*com < other.num*com


    def __gt__(self,other):
        com = gcd(self.den,other.den)
        return self.num*com > other.num*com





x = Fraction(3,4)
y= DFraction(1,2)
print x > y
