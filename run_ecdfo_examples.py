# -*- coding: utf-8 -*-
"""
#
#  Driver for the optimizer ECDFO
#  (Equality-Constrained Derivative-Free Optimization).
#
#  ECDFO can solve minimization problems of the form
#
#      minimize     f(x)
#      subject to   lx <=   x   <= ux
#                   ce(x) == 0,
#
#  where f: Rn -> R (hence x is the vector of n variables to optimize),
#  lx and ux are lower and upper bounds on x and ce: Rn -> Rme are
#  equality constraints.
#
"""
import helper
from runtime import *
from ecdfo_init_prob import ecdfo_init_prob_
from ecdfo_global_variables import set_prob, set_threshold,get_prob, set_check_condition
from ecdfo import ecdfo_
from evalfgh import evalfgh_
from numpy import array, zeros, arange


import time
tic = time.clock()

class options:
    pass

#---------------------------------------
# Initialize problem
#---------------------------------------

set_prob(5) #  definition of prob 1,...,5 in ecdfo_func(), extendable...
set_check_condition(1)
set_threshold(1e-08) # Threshold for violated bounds
prob=get_prob()
options.tol=zeros(3)

x,lb,ub,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(prob,nargout=13)

#---------------------------------------
# Set options
#---------------------------------------

options.hess_approx='model' # options: 'model' or 'bfgs'
options.bfgs_restart=0  # only taken into account if hess_approx is 'bfgs'
options.algo_descent='Powell' # options: 'Powell' or 'Wolfe'
options.tol[0]  = 1e-5       # tolerance on the gradient of the Lagrangian
options.tol[1]  = 1e-5       # tolerance on the feasibility
options.tol[2]  = 1e-5       # tolerance on the complementarity
options.miter   = 500        # max iterations
options.msimul  = 500        # max evaluations
options.verbose = 1          # verbosity level 0,...,3 (default: 1)

#------------------------------------
# Call ECDFO
#------------------------------------

lm=array([])
x,lm,info=ecdfo_(evalfgh_,x,lm,lb,ub,options,nargout=3)
print('x = '+str(x))
if info.ce.size > 0:
    print(' c = '+str(info.ce.T[0]))

toc = time.clock()
print("Elapsed time is " + str(toc - tic) + " seconds.")
