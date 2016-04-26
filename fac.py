# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 18:17:20 2016

@author: swee-i5-SSD
"""


def retfactorial(n):
    results = []
    print "factorial(%u) = " % n,
    F = factorial(n, results)
    print "= {:,}".format(F)
    return F, results


def factorial(n, results):
    if n > 1:
        print "%u x" % n,
        fn = n*factorial(n-1, results)
        results.append(fn)
        return fn
    else:
        print "%u" % n,
        results.append(1)
        return 1

def fib(n):
#    m={0:0, 1:1}

#    def fib_inner(inner_n):
#        try:
#            return m[inner_n]
#
#        except KeyError:
#            print "%u not in m" % inner_n
#            m[inner_n] = fib(inner_n-1) + fib(inner_n-2)
#            return m[inner_n]
#
#        finally:
#            print "returning Fib({})={}".format(inner_n, m[inner_n])
    memo = {0:0, 1:1}
    def fibm(n):
        print "Calling fibm({})".format(n)
        if n not in memo:
            print "%u not in memo" % n
            memo[n] = fibm(n-1) + fibm(n-2)
            print "returning fibm({})={}".format(n, memo[n])
        return memo[n]

    return fibm(n)

