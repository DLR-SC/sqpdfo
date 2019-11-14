# -*- coding: utf-8 -*-
from runtime import *
import sqpdfo_global_variables as glob
from sqpdfo import sqpdfo_

###############################################################################
#  Driver for the optimizer SQPDFO.
# (Sequential-Quadratic-Programming Derivative-Free Optimization).
#
#  SQPDFO can solve optimization problems of the form
#
#      minimize     f(x)
#      subject to   lx <= x <= ux
#      subject to   ce(x) == 0,
#                   ci(x) >=0
#
#  where f: Rn -> R (x is the vector of n variables to optimize),
#  lx and ux are lower and upper bounds on x, and ce: Rn -> Rme
#  are equality constraints and ci: Rn -> Rmi are inequality constraints.
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
        self.miter   = 500        # max iterations
        self.msimul  = 500        # max evaluations
        self.verbose = 1          # verbosity level 0,...,3 (default: 1)

if __name__ == '__main__':
    
    # call sqpdfo with academic example problems:
    # definition of function and constraints in sqpdfo_func_()
    # definition of starting values and bounds in sqpdfo_init_prob_()

    # set problem number (choose from 1...7 and 10...16)
    prob = 10
    glob.set_prob(prob)

    opts = optionsClass()
    
    x, lm, info = sqpdfo_(opts)
    
    # final printout
    print('x* = '+str(x))
    print(' f* = '+str(info.f))
    if info.ce.size > 0:
        print(' c* = '+str(info.ce.T[0]))        
