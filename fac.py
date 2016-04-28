# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 18:17:20 2016
Memoized versions for factorial and Fibonacci
@author: swee-i5-SSD
"""


def factorial(n):
    memo = {0: 1}

    def fac(N):
        print "calling Fac({})".format(N)
        try:
            return memo[N]
        except KeyError:
            memo[N] = N*fac(N-1)
            return memo[N]

    return fac(n), memo


def fib(n):
    memo = {0: 0, 1: 1}

    def fib_inner(N):
        print "calling Fib({})".format(N)
        try:
            return memo[N]

        except KeyError:
            memo[N] = fib_inner(N-1) + fib_inner(N-2)
            return memo[N]

    return fib_inner(n)
