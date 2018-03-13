# -*- coding: utf-8 -*-
"""
#
#  Driver for the optimizer ECDFO.
#  ECDFO can solve a minimization problem of the form
#
#      minimize     f(x)
#      subject to   lx <=   x   <= ux
#                   li <= ci(x) <= ui
#                   ce(x) == 0,
#
#  where f: Rn -> R (hence x is the vector of n variables to optimize),
#  lx and ux are lower and upper bounds on x, ci: Rn -> Rmi, li and ui
#  are lower and upper bounds on ci(x), and ce: Rn -> Rme.
#
"""
from sqpdfo import sqpdfo
from func_f import func_f
from func_c import func_c
from numpy import array
import time

# starting time
tic = time.clock()

# initialize start values
x0 = array([[-1.2,1.0]])
lb = array([[-5.,-5.]]).T
ub = array([[10.,10.]]).T

# call sqpdfo
#x,f,ce = sqpdfo(func_f,x0,lb,ub)
x,f,ce = sqpdfo(func_f,x0,lb,ub,func_c)

print('')
print('x* = '+str(x))
print('f* = '+str(f))

# stop timer
toc = time.clock()
print('')
print("Elapsed time is " + str(toc - tic) + " seconds.")
