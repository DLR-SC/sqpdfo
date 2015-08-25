# -*- coding: utf-8 -*-
"""
%
%  Driver for the optimizer ECDFO.
%  ECDFO can solve a minimization problem of the form
%
%      minimize     f(x)
%      subject to   lx <=   x   <= ux
%                   li <= ci(x) <= ui
%                   ce(x) == 0,
%
%  where f: Rn -> R (hence x is the vector of n variables to optimize),
%  lx and ux are lower and upper bounds on x, ci: Rn -> Rmi, li and ui
%  are lower and upper bounds on ci(x), and ce: Rn -> Rme.
%
"""
#from __future__ import division
#try:
import helper
from runtime import *
from ecdfo_init_prob import ecdfo_init_prob_
from ecdfo_global_variables import set_prob, set_threshold,get_prob, set_check_condition
from ecdfo import ecdfo_
from evalfgh import evalfgh_
from numpy import array, zeros, arange
#except ImportError:
#from smop.runtime import *

import time
tic = time.clock()

#clear(char('all'))
#close(char('all'))
#_format(char('long'))
#global n,nb,mi,me,prob,threshold

set_prob(7) #  definition of prob 1,...,5 in ecdfo_func(), extendable...
set_check_condition(1)
prob=get_prob()
options = helper.dummyUnionStruct()
options.tol=zeros(3)

x,lx,ux,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(prob,nargout=13)
lb=zeros_(n,1)
ub=zeros_(n,1)
lb[arange(0,n)]=lx
ub[arange(0,n)]=ux
if mi:
    lb[arange(n,n + mi)]=li
    ub[arange(n ,n + mi)]=ui
set_threshold(1e-08)
options.algo_method='quasi-Newton'
options.algo_globalization='trust regions'
options.hess_approx='model'
options.bfgs_restart=0
options.algo_descent='Powell'
if nb + mi + me == 0:
    options.algo_descent='Wolfe'
options.tol[0]=1e-05
options.tol[1]=1e-05
options.tol[2]=1e-05
options.dxmin=dxmin
options.miter=500
options.msimul=500
options.verbose=2
lm=array([])
x,lm,info=ecdfo_(evalfgh_,x,lm,lb,ub,options,nargout=3)
print x

toc = time.clock()
print "Elapsed time is " + str(toc - tic) + " seconds."
