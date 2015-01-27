# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 13:41:36 2015

@author: jaco_da
"""
from runtime import matlabarray
import numpy as np

#A = matlabarray([[1,2,3]]).T
#B = matlabarray([[10, 11, 12]]).T
#ilb = matlabarray([[True, False, True]]).T
#indfree = matlabarray([[0,1,2]]).T

#print "A\n", A
#print "B\n", B
#print "ilb\n", ilb
#print "indfree\n", indfree

#B[ilb] = A[indfree[ilb]]
#print "B after\n", B


lbounds = matlabarray([[-np.inf, -np.inf, -np.inf]]).T
ilb = matlabarray([[False, True, False]]).T
lb = matlabarray([[-0.5, 0.0, -np.inf]]).T
indfree = matlabarray([[1,2,3]]).T

print "lbounds\n", lbounds
print "ilb\n", ilb
print "lb\n", lb
print "indfree\n", indfree

lbounds[ilb]=lb[indfree[ilb]]

print "lboundsafter\n", lbounds