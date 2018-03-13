# -*- coding: utf-8 -*-
"""
#
#  Driver for the optimizer ECDFO.
#  ECDFO can solve a minimization problem of the form
#
#      minimize     f(x)
#      subject to   lx <=   x   <= ux
#                   ce(x) == 0,
#
#  where f: Rn -> R (hence x is the vector of n variables to optimize),
#  lx and ux are lower and upper bounds on x, and ce: Rn -> Rme.
#
"""
# python sqpdfo(func_f,[-1.2,1.0],[-5.,-5.],[10.,10.])
# python sqpdfo(func_f,[-1.2,1.0],[-5.,-5.],[10.,10.],func_c)

import helper
from runtime import *
import ecdfo_global_variables as glob
from ecdfo import ecdfo_
from evalfgh import evalfgh_
from numpy import array, zeros, arange, shape

def sqpdfo(func=None, x=None, lx=None, ux=None, cfunc=None, *args, **kwargs):

    # Set problem global variables
    glob.set_prob(100)                       # set problem number (100 is userdefined problem)
    if func == None:
        print('Error: Definition of the objective function (in func_f.py) is missing !')
    else:
        glob.set_filename_f(func)
    if cfunc == None:
        glob.set_filename_ce('')
    else:
        glob.set_filename_ce(cfunc)
    glob.set_check_condition(1)
    glob.set_fileoutput(1)
    glob.set_simul_not_initialized(1)
    glob.set_threshold(1e-08)                # Threshold for violated bounds

    #---------------------------------------
    # Initialize problem
    #---------------------------------------

    options = helper.dummyUnionStruct()
    options.tol = zeros(3)
    lm = array([])
    n = max(shape(x))

    # Check right format of x and bounds lx and ux
    if not isinstance(x, np.ndarray):
        x = array([x])
    if not isinstance(lx, np.ndarray):
        lx = array([lx]).T
    if not isinstance(ux, np.ndarray):
        ux = array([ux]).T

    # bounds for variables
    lb = zeros_(n, 1)
    ub = zeros_(n, 1)
    lb[arange(0,n)]=lx
    ub[arange(0,n)]=ux

    #---------------------------------------
    # Set options
    #---------------------------------------

    # default options
    options.hess_approx='model'  # options: 'model' or 'bfgs'
    options.bfgs_restart=0       # only taken into account if hess_approx is 'bfgs'
    options.algo_descent='Powell' # options: 'Powell' or 'Wolfe'
    options.tol[0]  = 1e-3       # tolerance on the gradient of the Lagrangian
    options.tol[1]  = 1e-2       # tolerance on the feasibility
    options.tol[2]  = 1e-5       # tolerance on the complementarity
    options.dxmin   = 1e-8       # minimum size of a step
    options.miter   = 1000        # max iterations
    options.msimul  = 1000        # max evaluations
    options.verbose = 2          # verbosity level 0,...,3

    #------------------------------------
    # Call ECDFO
    #------------------------------------

    x,lm,info=ecdfo_(evalfgh_,x,lm,lb,ub,options,nargout=3)

    # Return values
    f = info.f
    ce = info.ce
    return x, f, ce