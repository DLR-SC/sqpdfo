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
import helper
from runtime import *
from ecdfo_init_prob import ecdfo_init_prob_
from ecdfo_global_variables import set_prob, set_threshold,get_prob, set_check_condition
from ecdfo import ecdfo_
from evalfgh import evalfgh_
from numpy import array, zeros, arange


import time
tic = time.clock()

#clear(char('all'))
#close(char('all'))
#_format(char('long'))
#global n,nb,mi,me,prob,threshold



#---------------------------------------
# Initialize problem
#---------------------------------------

set_prob(5) #  definition of prob 1,...,5 in ecdfo_func(), extendable...
set_check_condition(1)
prob=get_prob()
options = helper.dummyUnionStruct()
options.tol=zeros(3)

x,lx,ux,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(prob,nargout=13)
lb=zeros_(n,1)
ub=zeros_(n,1)
# bounds for variables and inequality constraints

lb[arange(0,n)]=lx
ub[arange(0,n)]=ux
if mi:
    lb[arange(n,n + mi)]=li
    ub[arange(n ,n + mi)]=ui
# Threshold for violated bounds 

set_threshold(1e-08)
#---------------------------------------
# Set options
#---------------------------------------

# default options
options.algo_method='quasi-Newton'
options.algo_globalization='trust regions'
options.hess_approx='model' # options: 'model' or 'bfgs'
options.bfgs_restart=0  # only taken into account if hess_approx is 'bfgs'
options.algo_descent='Powell' # options: 'Powell' or 'Wolfe'
if nb + mi + me == 0:
    options.algo_descent='Wolfe'

options.tol[0]  = 1e-5;       # tolerance on the gradient of the Lagrangian
options.tol[1]  = 1e-5;       # tolerance on the feasibility
options.tol[2]  = 1e-5;       # tolerance on the complementarity

options.dxmin   = dxmin;      # minimum size of a step

options.miter   = 500;        # max iterations
options.msimul  = 500;        # max evaluations

options.verbose = 2;          # verbosity level 0,...,3

#------------------------------------
# Call ECDFO
#------------------------------------

lm=array([])
x,lm,info=ecdfo_(evalfgh_,x,lm,lb,ub,options,nargout=3)
print(x)

toc = time.clock()
print("Elapsed time is " + str(toc - tic) + " seconds.")
