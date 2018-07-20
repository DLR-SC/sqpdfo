# -*- coding: utf-8 -*-
from runtime import *
import ecdfo_global_variables as glob
from ecdfo import ecdfo_
from numpy import array, zeros, shape

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


class optionsClass:
    def __init__(self):
        self.hess_approx = 'model'  # Type of Hessian approximation, options: 'model' or 'bfgs'
        self.bfgs_restart = 0  # Restart of the BFGS Hessian approximation after how many iterations
        # only taken into account if hess_approx is 'bfgs'
        self.whichmodel = 'subbasis'  # options: 'Subbasis', 'Frobnorm', 'L2norm','Regression'
        # only taken into account if hess_approx is 'model'
        self.final_degree = 'diagonal'  # Final degree of the model, options: 'linear', 'diagonal', 'quadratic'
        self.tol_grad = 1e-5  # tolerance on the gradient of the Lagrangian
        self.tol_feas = 1e-5  # tolerance on the feasibility
        self.tol_bnds = 1e-5  # tolerance on the bounds
        self.miter = 500  # max iterations
        self.msimul = 500  # max evaluations
        self.verbose = 1  # verbosity level 0,...,3 (default: 1)


def run_ecdfo(func=None, x=None, lx=None, ux=None, cfunc=None, *args, **kwargs):

    # Set problem global variables
    glob.set_prob(100)                       # set problem number to 100 (=userdefined problem which reads func_f() and func_c())
    if func is None:
        print('Error: Definition of the objective function (in func_f()) is missing !')
    else:
        glob.set_filename_f(func)
    if cfunc is None:
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

    options = optionsClass()

    #------------------------------------
    # Call ECDFO
    #------------------------------------

    x,lm,info = ecdfo_(x,lm,lx,ux,options)

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
