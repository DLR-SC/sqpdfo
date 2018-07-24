# -*- coding: utf-8 -*-
from runtime import *
import ecdfo_global_variables as glob
from ecdfo import ecdfo_
from numpy import array, ones_like, infty

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
        self.final_degree = 'diagonal'  # Final degree of the model, options: 'linear', 'diagonal', 'quadratic'
        self.tol_grad = 1e-5  # tolerance on the gradient of the Lagrangian
        self.tol_feas = 1e-5  # tolerance on the feasibility
        self.tol_bnds = 1e-5  # tolerance on the bounds
        self.miter = 500  # max iterations
        self.msimul = 500  # max evaluations
        self.verbose = 1  # verbosity level 0,...,3 (default: 1)


def run_ecdfo(func=None, x=None, lx=None, ux=None, cfunc=None, *args, **kwargs):

    # Set global variables
    glob.set_prob(100)  # set problem number to 100 (=userdefined problem which reads func() and cfunc())
    glob.set_filename_f(func)
    if cfunc is None:
        glob.set_filename_ce('')
    else:
        glob.set_filename_ce(cfunc)

    #---------------------------------------
    # Initialize problem
    #---------------------------------------

    # initial Lagrange multipliers for the problem (not necessary)
    lm = array([])

    # Check right format of x and bounds lx and ux
    if x is None:
        sys.exit('Error: No starting values of x given !')
    elif not isinstance(x, np.ndarray):
        x = array([x])
    if lx is None:
        lx = -1e20 * ones_like(x).T
    elif not isinstance(lx, np.ndarray):
        lx = array([lx]).T
    if ux is None:
        ux = 1e20 * ones_like(x).T
    elif not isinstance(ux, np.ndarray):
        ux = array([ux]).T

    # Set options
    options = optionsClass()

    #------------------------------------
    # Call ECDFO
    #------------------------------------

    x,lm,info = ecdfo_(x,lm,lx,ux,options)
    return x,info


def func_f(xvector):
    # 2D Rosenbrock function (constrained on the unitdisk if func_c() is considered)
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
    #x,info = run_ecdfo(func_f,x0)
    #x,info = run_ecdfo(func_f,x0,lb,ub)
    x,info = run_ecdfo(func_f,x0,lb,ub,func_c)

    # final printout
    print('')
    print('x* = '+str(x))
    print('f* = '+str(info.f))
    if info.ce.size>0:
        print('ce* = '+str(info.ce.T[0]))
