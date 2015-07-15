#!/usr/local/bin/python

import numpy
import math

def norm( x, p ):
   if p == '1':
      norm = abs(x).sum()
   elif p == '2':
      y = x**2
      norm = sqrt(y.sum())
   elif p == 'inf':
      norm = abs(x).max()
   else:
      raise "Norm kind p is wrong in btr.py!"
   return norm

def isarray(a):
    """
    Test for arrayobjects. Can also handle UserArray instances
    """
    try:
        sh = list(a.shape)
    except AttributeError:
        return 0
    try:
        sh[0] = sh[0]+1
        a.shape = sh
    except ValueError:
        return 1
    except IndexError:
        return 1 # ? this is a scalar array
    return 0
