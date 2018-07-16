# -*- coding: utf-8 -*-
###############################################################################
#  Driver for the optimizer ECDFO.
#  ECDFO can solve a minimization problem of the form
#
#      minimize     f(x)
#      subject to   lx <= x <= ux
#      subject to   ce(x) == 0,
#
#  where f: Rn -> R (x is the vector of n variables to optimize),
#  lx and ux are lower and upper bounds on x, and ce: Rn -> Rme.
#
#  ---------------------------------------------------------
#  See example call to ECDFO in '__main__' below !!
###############################################################################

from runtime import *
import ecdfo_global_variables as glob
from ecdfo import ecdfo_
from evalfgh import evalfgh_
from numpy import array, zeros, shape


def run_ecdfo(func=None, x=None, lx=None, ux=None, cfunc=None, *args, **kwargs):

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

    class options:
        pass
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

    #---------------------------------------
    # Set options
    #---------------------------------------

    # default options
    options.hess_approx='model'  # options: 'model' or 'bfgs'
    options.bfgs_restart=0       # only taken into account if hess_approx is 'bfgs'
    options.algo_descent='Powell' # options: 'Powell' or 'Wolfe'
    options.tol[0]  = 1e-5       # tolerance on the gradient of the Lagrangian
    options.tol[1]  = 1e-5       # tolerance on the feasibility of the equality constraints
    options.tol[2]  = 1e-5       # tolerance on the bounds
    options.dxmin   = 1e-8       # minimum size of a step
    options.miter   = 1000        # max iterations
    options.msimul  = 1000        # max evaluations
    options.verbose =2          # verbosity level 0,...,3

    #------------------------------------
    # Call ECDFO
    #------------------------------------

    x,lm,info = ecdfo_(evalfgh_,x,lm,lx,ux,options,nargout=3)

    # Return values
    f = info.f
    ce = info.ce
    return x, f, ce


def func_f(xvector):
    # 2D Rosenbrock function (constrained on the unitdisk if func_c() is considered)
    # x* = (1.0, 1.0)        (constrained on the unitdisk: x* = (0.7864, 0.6177))
    # f* = 0.0               (constrained on the unitdisk: f* = 0.045674824758137236)
    f =  (1-xvector[0])**2 + 100*(xvector[1]-xvector[0]**2)**2
    msg = 0
    return f,msg


def func_c(xvector):
    # solution on the unitdisk
    ce = np.zeros(1)
    ce[0] = xvector[0]**2 + xvector[1]**2 - 1
    msg = 0
    return ce, msg


if __name__ == '__main__':

    # initialize starting point and bounds
    x0 = array([[2.5,1.0]])
    lb = array([[-5.,-5.]]).T
    ub = array([[10.,10.]]).T

    # call run_ecdfo
    #x,f,ce = run_ecdfo(func_f,x0,lb,ub)
    x,f,ce = run_ecdfo(func_f,x0,lb,ub,func_c)

    # final printout
    print('')
    print('x* = '+str(x))
    print('f* = '+str(f))
    if ce.size>0:
        print('ce* = '+str(ce.T[0]))
