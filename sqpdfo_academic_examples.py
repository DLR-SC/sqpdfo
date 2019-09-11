# -*- coding: utf-8 -*-
import ecdfo_global_variables as glob
from ecdfo_init_prob import ecdfo_init_prob_
from ecdfo import ecdfo_
from numpy import array, zeros, ones, concatenate
from runtime import *
import time

###############################################################################
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
###############################################################################

class optionsClass:
    def __init__(self):
        self.hess_approx='model'      # Type of Hessian approximation, options: 
                                      # 'model' or 'bfgs'
        self.bfgs_restart=0           # Restart of the BFGS Hessian approximation 
                                      # after how many iterations
                                      # only taken into account if hess_approx is 'bfgs'
        self.whichmodel = 'subbasis'  # options: 'Subbasis', 'Frobnorm', 'L2norm', 
                                      # 'Regression'
        self.final_degree = 'quadratic' # Final degree of the model, options: 
                                        # 'linear', 'diagonal', 'quadratic'
        self.tol_grad  = 1e-5       # tolerance on the gradient of the Lagrangian
        self.tol_feas  = 1e-5       # tolerance on the feasibility
        self.tol_bnds  = 1e-5       # tolerance on the bounds
        self.miter   = 200        # max iterations
        self.msimul  = 200        # max evaluations
        self.verbose = 1          # verbosity level 0,...,3 (default: 1)

def sqpdfo_academic_examples(prob=None):

    options = optionsClass()

    # Initialize problem
    glob.set_prob(prob) 
    x,lb,ub,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(prob,nargout=13)

    # Handling of inequalities if any
    if mi:
        glob.set_nbr_slacks(mi)
        glob.set_slacks(zeros((mi,1)))
        lb = concatenate((lb,zeros((mi,1))))
        ub = concatenate((ub,1e20*ones((mi,1))))
    
    # Call ECDFO
    lm=array([])
    x,lm,info=ecdfo_(x, lm, lb, ub, options)
    return x, info
    
if __name__ == '__main__':    
    
    # call sqpdfo with academic example problems:
    # definition of function and constraints in ecdfo_func_()
    # definition of starting values in ecdfo_init_prob_()
    
    prob = 16
    x, info = sqpdfo_academic_examples(prob)
    
    # final printout
    print('x* = '+str(x))
    print(' f* = '+str(info.f))
    if info.ce.size > 0:
        print(' c* = '+str(info.ce.T[0]))        
