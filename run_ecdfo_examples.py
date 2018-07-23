# -*- coding: utf-8 -*-
import ecdfo_global_variables as glob
from ecdfo_init_prob import ecdfo_init_prob_
from ecdfo import ecdfo_
from numpy import array
import time

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

tic = time.clock()

#---------------------------------------
# Initialize problem
#---------------------------------------

glob.set_prob(3) #  definition of prob 1,...,5 in ecdfo_func_(), extendable...
x,lb,ub,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(glob.get_prob(),nargout=13)

#---------------------------------------
# Set options
#---------------------------------------

class optionsClass:
    def __init__(self):
        self.hess_approx='model'      # Type of Hessian approximation, options: 'model' or 'bfgs'
        self.bfgs_restart=0           # Restart of the BFGS Hessian approximation after how many iterations
                                      # only taken into account if hess_approx is 'bfgs'
        self.whichmodel = 'subbasis'  # options: 'Subbasis', 'Frobnorm', 'L2norm','Regression'
                                      # only taken into account if hess_approx is 'model'
        self.final_degree = 'quadratic' # Final degree of the model, options: 'linear', 'diagonal', 'quadratic'
        self.tol_grad  = 1e-5       # tolerance on the gradient of the Lagrangian
        self.tol_feas  = 1e-5       # tolerance on the feasibility
        self.tol_bnds  = 1e-5       # tolerance on the bounds
        self.miter   = 500        # max iterations
        self.msimul  = 500        # max evaluations
        self.verbose = 1          # verbosity level 0,...,3 (default: 1)

options = optionsClass()

#------------------------------------
# Call ECDFO
#------------------------------------

lm=array([])
x,lm,info=ecdfo_(x, lm, lb, ub, options)
print('x* = '+str(x))
print(' f* = '+str(info.f))
if info.ce.size > 0:
    print(' c* = '+str(info.ce.T[0]))

toc = time.clock()
print("Elapsed time is " + str(toc - tic) + " seconds.")
