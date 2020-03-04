# -*- coding: utf-8 -*-
#from runtime import *
import sqpdfo_global_variables as glob
import sqpdfo
import pycutest

###############################################################################
#  Driver for the optimizer SQPDFO to solve problems from the CUTEst library.
# (Sequential-Quadratic-Programming Derivative-Free Optimization)
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
#
#  pycutest can be downloaded from https://jfowkes.github.io/pycutest/
#
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
        self.miter = 5000  # max iterations
        self.msimul = 5000  # max evaluations
        self.verbose = 1  # verbosity level 0,...,3 (default: 1)

if __name__ == '__main__':

    # set problem number (1000 is problem from CUTEst library)
    glob.set_prob(1000)  # Do not modify this number!

    # set Cutest problem with name and parameters of the problem
    glob.set_prob_cuter('AVGASA',[])
    #glob.set_prob_cuter('ALJAZZAF', ['N', 3, 'N1', 2])
    #glob.set_prob_cuter('BDRY2', ['N', 3])
    #glob.set_prob_cuter('ALLINQP',['N',10])

    # Set options
    opts = optionsClass()
    
    # call run_sqpdfo_cutest
    x,lm,info = sqpdfo.optimize(opts)

    # final printout
    print('')
    print('x* = '+str(x))
    print('f* = '+str(info.f))
    if info.ce.size>0:
        print('c* = '+str(info.ce.T[0]))
