# -*- coding: utf-8 -*-
"""
Created on Mon Jul 04 18:02:44 2016

@author: sweel_Rafelski
"""
import functools

class decoratorWithArguments(object):

    def __init__(self, arg1, arg2, arg3):
        """
        If there are decorator arguments, the function
        to be decorated is not passed to the constructor!
        """
        print "Inside __init__()"
        self.arg1 = arg1
        self.arg2 = arg2
        self.arg3 = arg3

    def __call__(self, f):
        functools.update_wrapper(self, f)
        """
        If there are decorator arguments, __call__() is only called
        once, as part of the decoration process! You can only give
        it a single argument, which is the function object.
        """
        print "Inside __call__()"

        def wrapped_f(*args):
            print "called {}".format(self.__name__)
            print "Decorator arguments:", self.arg1, self.arg2, self.arg3
            f(*args)
            print "After f(*args)"
        return wrapped_f

@decoratorWithArguments(1,2,3)
def sayHello(a1, a2, a3, a4):
    print 'sayHello arguments:', a1, a2, a3, a4

sayHello('a', 'ab' , 'c', 'd')