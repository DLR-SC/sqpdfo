# -*- coding: utf-8 -*-
from runtime import *
import ecdfo_global_variables as glob
from ecdfo import ecdfo_
from numpy import array, ones_like, infty
from ecdfo_init_prob import ecdfo_init_prob_
import pycutest

###############################################################################
#  Driver for the optimizer ECDFO to solve problems from the CUTEst library.
#
#  pycutest is downloadable from https://jfowkes.github.io/pycutest/
#
###############################################################################
#  ECDFO can solve a minimization problem of the form
#
#      minimize     f(x)
#      subject to   lx <= x <= ux
#      subject to   ce(x) == 0,
#
#  where f: Rn -> R (x is the vector of n variables to optimize),
#  lx and ux are lower and upper bounds on x, and ce: Rn -> Rme.
#
#  -------------------------------------------------------------
#  Example calls to ECDFO:
#      x,f    = run_ecdfo(func_f, x0)
#      x,f    = run_ecdfo(func_f, x0, lb, ub)
#      x,f,ce = run_ecdfo(func_f, x0, lb, ub, func_c)
###############################################################################


class optionsClass:
    def __init__(self):
        self.hess_approx = 'model'  # Type of Hessian approximation, options: 'model' or 'bfgs'
        self.bfgs_restart = 0  # Restart of the BFGS Hessian approximation after how many iterations
                               # only taken into account if hess_approx is 'bfgs'
        self.whichmodel = 'subbasis'  # options: 'Subbasis', 'Frobnorm', 'L2norm','Regression'
        self.final_degree = 'quadratic'  # Final degree of the model, options: 'linear', 'diagonal', 'quadratic'
        self.tol_grad = 1e-4  # tolerance on the gradient of the Lagrangian
        self.tol_feas = 1e-4  # tolerance on the feasibility
        self.tol_bnds = 1e-4  # tolerance on the bounds
        self.miter = 500  # max iterations
        self.msimul = 500  # max evaluations
        self.verbose = 1  # verbosity level 0,...,3 (default: 1)


def run_ecdfo_cutest(*args, **kwargs):

    # Set global variables
    glob.set_prob(1000)  # set problem number to 1000 (=problem from CUTEst library)
    glob.set_check_condition(0)
    #---------------------------------------
    # Initialize problem
    #---------------------------------------
    
    prob = glob.get_prob()
    x,lx,ux,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(prob)

    # initial Lagrange multipliers for the problem (not necessary)
    lm = array([])

    # Set options
    options = optionsClass()

    #------------------------------------
    # Call ECDFO
    #------------------------------------

    x,lm,info = ecdfo_(x,lm,lx,ux,options)
    return x,info

if __name__ == '__main__':

    #glob.set_prob_cuter(name,params)
    #glob.set_prob_cuter('AIRCRFTA',[])
    #glob.set_prob_cuter('ALJAZZAF', ['N', 3, 'N1', 2])
    #glob.set_prob_cuter('ARTIF', ['N', 10])
    glob.set_prob_cuter('BT13',[])
    
    # call run_ecdfo_cutest
    x,info = run_ecdfo_cutest()

    # final printout
    print('')
    print('x* = '+str(x))
    print('f* = '+str(info.f))
    if info.ce.size>0:
        print('ce* = '+str(info.ce.T[0]))
