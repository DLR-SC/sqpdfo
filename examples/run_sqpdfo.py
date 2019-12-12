# -*- coding: utf-8 -*-

from sqpdfo.sqpdfo import sqpdfo_
import numpy as np

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
        self.hess_approx = 'model'  # Type of Hessian approximation, options: 
                                    # 'model' or 'bfgs'
        self.bfgs_restart = 0  # Restart of the BFGS Hessian approximation after 
                               # how many iterations
                               # (only taken into account if hess_approx is 'bfgs')
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

def func_f(xvector):
    # 2D Rosenbrock function
    f =  (1-xvector[0])**2 + 100*(xvector[1]-xvector[0]**2)**2
    return f


def func_c(xvector):
    # equalities (ce == 0) first !!
    ce = np.zeros(2)
    ce[0] = xvector[0]**2 + xvector[1]**2 - 2
    ce[1] = - (xvector[0] - 1)**3 + xvector[1] - 1
    # inequalities (ci >= 0)
    ci = np.zeros(1)
    ci[0] = - xvector[0] - xvector[1] + 2    
    c = np.concatenate((ce,ci))  # equalities first !!!
    return c

if __name__ == '__main__':

    # initialize start values
    x0 = np.array([[-1.2,1.0]])
    lb = np.array([[-5.,-5.]]).T
    ub = np.array([[10.,10.]]).T
    
    me = 2
    mi = 1

    options = optionsClass()

    # call SQPDFO
    if me+mi > 0:
        x, lm, info = sqpdfo_(options, func_f, x0, lb, ub, me, mi, func_c)
    else:
        x, lm, info = sqpdfo_(options, func_f, x0, lb, ub)

    # printout results
    print('')
    print('x* = '+str(x))
    print('f* = '+str(info.f))
    if info.ce.size>0:
        print('c* = '+str(info.ce.T[0]))
