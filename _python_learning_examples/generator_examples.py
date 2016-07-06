# -*- coding: utf-8 -*-
"""
Created on Sun Apr 03 18:01:23 2016
@author: sweel_rafelski
"""
import glob
import os
import sys

def readfiles(filenames):
    for f in filenames:
        count=0
        for line in open(f):
            yield line, count, f

#            print "now on file %s, line %d " % (f, count)
#            count+=1

def grep(pattern, lines):
    return ((line, count, f) for line, count, f in lines if pattern in line)

def printlines(filtlines):
    for line, count, f in filtlines:
        print ("Matched in file %s\nLine %d\n%s" % (f, count, line))

# main(pattern, filenames):
if __name__ == '__main__':
    filenames = glob.glob('*.vtk')
    # this generator FUNCTION has all lines in all opened files
    genlines = readfiles(filenames)
    # this gen. EXPRESSION only lines containing pattern
    flines = grep("vtk", genlines)
    # this FUNCTION iterates over the gen expression
    printlines(flines)
    print lines
