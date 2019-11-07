# -*- coding: utf-8 -*-
from runtime import *
import ecdfo_global_variables as glob
from ecdfo import ecdfo_
from numpy import array, zeros, ones, arange, shape, concatenate

###############################################################################
#  Driver for the optimizer ECDFO.
#  ECDFO can solve a minimization problem of the form
#
#      minimize     f(x)
#      subject to   lx <= x <= ux
#      subject to   ce(x) == 0,
#
#  where f: Rn -> R (hence x is the vector of n variables to optimize),
#  lx and ux are lower and upper bounds on x, and ce: Rn -> Rme.
###############################################################################

class optionsClass:
    def __init__(self):
        self.hess_approx = 'model'  # Type of Hessian approximation, options: 
                                    # 'model' or 'bfgs'
        self.bfgs_restart = 0  # Restart of the BFGS Hessian approximation after 
                               # how many iterations
                               # only taken into account if hess_approx is 'bfgs'
        self.whichmodel = 'subbasis'  # options: 'Subbasis', 'Frobnorm', 'L2norm', 
                                      # 'Regression'
        self.final_degree = 'quadratic'  # Final degree of the model, options: 
                                        # 'linear', 'diagonal', 'quadratic'
        self.tol_grad = 1e-4  # tolerance on the gradient of the Lagrangian
        self.tol_feas = 1e-4  # tolerance on the feasibility
        self.tol_bnds = 1e-4  # tolerance on the bounds
        self.miter = 500  # max iterations
        self.msimul = 500  # max evaluations
        self.verbose = 1  # verbosity level 0,...,3 (default: 1)

def sqpdfo(func=None, x=None, lx=None, ux=None, me=None, mi=None, cfunc=None, *args, **kwargs):

    # Set options
    options = optionsClass()

    # Set problem global variables
    glob.set_prob(100)               # set problem number (100 is userdefined problem)
    if func == None:
        print('Error: Definition of the objective function (in func_f.py) is missing !')
    else:
        glob.set_filename_f(func)
    if cfunc == None:
        glob.set_filename_cons('')
    else:
        glob.set_filename_cons(cfunc)

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

    # Handling of inequalities if any
    if mi:
        glob.set_nbr_slacks(mi)
        glob.set_slacks(zeros((mi,1)))
        lx = concatenate((lx,zeros((mi,1))))
        ux = concatenate((ux,1e20*ones((mi,1))))

    # Call ECDFO
    lm = array([])
    x,lm,info=ecdfo_(x,lm,lx,ux,options)

    # Return values
    f = info.f
    c = info.ce
    return x, f, c
    
def func_f(xvector):
    # 2D Rosenbrock function
    f =  (1-xvector[0])**2 + 100*(xvector[1]-xvector[0]**2)**2
    msg = 0
    return f,msg


def func_c(xvector):
    # equalities (ce == 0) first !!
    ce = np.zeros(2)
    ce[0] = xvector[0]**2 + xvector[1]**2 - 2
    ce[1] = - (xvector[0] - 1)**3 + xvector[1] - 1
    # inequalities (ci >= 0)
    ci = np.zeros(1)
    ci[0] = - xvector[0] - xvector[1] + 2    
    c = concatenate((ce,ci))  # equalities first !!!
    msg = 0
    return c, msg

if __name__ == '__main__':

    # initialize start values
    x0 = array([[-1.2,1.0]])
    lb = array([[-5.,-5.]]).T
    ub = array([[10.,10.]]).T
    
    me = 1; mi = 2

    # call sqpdfo
    #x,f,c = sqpdfo(func_f,x0,lb,ub)
    x,f,c = sqpdfo(func_f,x0,lb,ub,me,mi,func_c)

    # printout results
    print('')
    print('x* = '+str(x))
    print('f* = '+str(f))
    if c.size>0:
        print('c* = '+str(c.T[0]))
