# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 14:08:04 2014

@author: jaco_da
"""
#import sys
import numpy as np
from runtime import *
#from bcdfo_build_QR_of_Y import *
#from bcdfo_solve_TR_MS_bc import *

X = matlabarray([1,2,3]).T

a = matlabarray([2.0, 3.0, 7.0])
b = matlabarray([4.0, 3.0, 3.0])

X[:,2] = a

#X[:,2] = b

#X[4,:] = matlabarray([7.0])


#X = np.resize(X, (3,2))

print "X", X
#print np.multiply(a, b)
#M = matlabarray([[1,2],[3,4]])
#print a.T*M*a


#print strcmp_(char('Hello'), char('Hello'))

#A = matlabarray([[1,     2,     3],
     #[1,     2,     3],
     #[1,     2,     3]])

#A = matlabarray([[1, 1, 1],
# [1, 2, 3], [0,0,0]])

#A = matlabarray([
#    [2,3,5],
#    [-4,2,3]
#    ])  
	
#print "np", np.linalg.solve(np.eye(3), np.zeros(3))	

#Z = null_(A)
#print "norm", norm_(Z)
#print "Z\n", Z
#print "Z^TZ\n", Z.T*Z

#gx = matlabarray([2, 3]).T
#print "gx", gx
##H = matlabarray(np.array([[4, 6],[6, 5]]))
#H = matlabarray([[4, 6],[6, 5]])
#print "H", H
##lb = np.array([-10, -10]).T
#lb = matlabarray([-10, -10]).T
#print "lb", lb
##ub = np.array([0, 0]).T
#ub = matlabarray([0, 0]).T
#print "ub", ub
##bcdfo_solve_TR_MS_bc( [2; 3], [4 6; 6 5], [-10; -10], [0; 0], 1.0, 0.001, 1 )
#ret2 = "fertig"
#ret1, ret2, _, _, _, _, _, _  = bcdfo_solve_TR_MS_bc( gx, H, lb, ub, matlabarray([1.0]), matlabarray([0.001]), matlabarray([1]) )
#print "ret1", ret1
#print "type ret1", type(ret1)
#print "ret2", ret2
#  should give
#    0
#   -0.6

