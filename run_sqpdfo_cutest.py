# -*- coding: utf-8 -*-
from runtime import *
import sqpdfo_global_variables as glob
from sqpdfo import sqpdfo_
from numpy import array, ones_like, infty, zeros, ones, concatenate
from sqpdfo_init_prob import sqpdfo_init_prob_
import pycutest
import warnings

###############################################################################
#  Driver for the optimizer SQPDFO to solve problems from the CUTEst library.
#
#  pycutest is downloadable from https://jfowkes.github.io/pycutest/
#
###############################################################################
#  SQPDFO can solve a minimization problem of the form
#
#      minimize     f(x)
#      subject to   lx <= x <= ux
#      subject to   ce(x) == 0,
#
#  where f: Rn -> R (x is the vector of n variables to optimize),
#  lx and ux are lower and upper bounds on x, and ce: Rn -> Rme.
#
#  -------------------------------------------------------------
#  Example calls to SQPDFO:
#      x,f    = run_sqpdfo(func_f, x0, lb, ub, None, None, None, options)
#      x,f,c = run_sqpdfo(func_f, x0, lb, ub, me, mi, func_c, options)
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
        self.miter = 1500  # max iterations
        self.msimul = 5000  # max evaluations
        self.verbose = 1  # verbosity level 0,...,3 (default: 1)


def run_sqpdfo_cutest(*args, **kwargs):

    # Set global variables
    glob.set_prob(1000)  # set problem number to 1000 (=problem from CUTEst library)
    glob.set_check_condition(0)
    
    #---------------------------------------
    # Initialize problem
    #---------------------------------------
    
    prob = glob.get_prob()
    x,lx,ux,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=sqpdfo_init_prob_(prob)
    
    # initial Lagrange multipliers for the problem (not necessary)
    lm = array([])

    # Set options
    options = optionsClass()
    
    # Handling of inequalities if any
    if mi:
    
        # handling of two-sided inequalities in CUTEst
        for i in range(0,mi+me):
            if li[i] < ui[i]: 
                if li[i] > -1e7 and ui[i] < 1e7:
                    mi = mi+1
    
        glob.set_nbr_slacks(mi)
        glob.set_slacks(zeros((mi,1)))
        lx = concatenate((lx,zeros((mi,1))))
        ux = concatenate((ux,1e20*ones((mi,1))))
    

    #------------------------------------
    # Call SQPDFO
    #------------------------------------

    x,lm,info = sqpdfo_(x,lm,lx,ux,options)
    return x,info

if __name__ == '__main__':

    #glob.set_prob_cuter(name,params)
    glob.set_prob_cuter('TRUSPYR2',[])
    #glob.set_prob_cuter('ALJAZZAF', ['N', 3, 'N1', 2])
    #glob.set_prob_cuter('BDRY2', ['N', 3])
    #glob.set_prob_cuter('ALLINQP',['N',10])
    
    # call run_sqpdfo_cutest
    x,info = run_sqpdfo_cutest()

    # final printout
    print('')
    print('x* = '+str(x))
    print('f* = '+str(info.f))
    if info.ce.size>0:
        print('c* = '+str(info.ce.T[0]))
