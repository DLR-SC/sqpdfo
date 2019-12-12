# -*- coding: utf-8 -*-

from sqpdfo.sqpdfo_check_convex import *
from sqpdfo.sqpdfo_check_cond import *
from sqpdfo.blls import *
from copy import copy
import sqpdfo.sqpdfo_global_variables as glob
from numpy import array, zeros, ones, concatenate


def sqpdfo_compute_multiplier_(x=None,lb=None,ub=None,info_=None,\
    options=None,values=None,*args,**kwargs):

    ###############################################################################
    # This procedure computes exact least-squares multipliers 'lm'. 
    # It solves in lm the quadratic optimization problem:
    #
    #   min || g+A'*lm ||^2
    #   subject to possible bounds on lm,
    #
    # where g = gradient, lm = lm(1:n+me) and A' = [ones(n) info.ae'] 
	# the Jacobian of the problem.
    #
    # A multiplier associated with an inequality constraint having
    # - infinite lower and upper bounds vanishes,
    # - infinite lower bound and finite upper bound is nonnegative,
    # - finite lower bound and infinite upper bound is nonpositive,
    # - finite lower bound and finite upper bound can have any sign.
    #
    # Output:
    #   lm: computed least-squares multiplier, more precisely
    #       lm(1:n): multiplier associated with the bounds on the variables
    #       lm(n+1:n+me): multiplier associated with the me equality
    #           constraints
    ###############################################################################

    info = copy(info_)
    me = size_(info.ae, 1)
    
    nbr_slacks = glob.get_nbr_slacks()
    sl = glob.get_slacks()
    
    glocal = copy(info.g)
    gconstraints = copy(info.ae)
    
    if nbr_slacks > 0:
       glocal = concatenate((glocal,zeros((nbr_slacks,1))))
       hh1 = concatenate((zeros((me-nbr_slacks,nbr_slacks)),-2*diag(sl.T[0])))
       gconstraints = concatenate((gconstraints,hh1),axis=1)

    nargin = 6 - [x is None, lb is None, ub is None, info is None, options is None,\
             values is None].count(True) + len(args)

    # If somewhere we have done set_check_condition (in the tests for instance we 
    # have set_check_condition(0)), then
    # we get this value, otherwise we take '1' by default.
    try:
        check_condition = glob.get_check_condition()
    except:
        check_condition = 1

    # Initialization

    lm = array([])
    info.flag = values.success
    badcond = 0

    # Dimensions

    n = length_(info.g)
    me = 0
    if (nargin >= 3):
        me = size_(info.ae, 1)

    # Check input arguments

    if (nargin < 4) or isempty_(lb):
        lb = - options.inf * ones_(n, 1)
    else:
        lb = lb[:]
        if any_(size_(lb) != [n, 1]):
            fprintf_('\n### sqpdfo_compute_multiplier: incorrect size of lb\n\n')
            info.flag = values.fail_unexpected
            return lm, info
    if (nargin < 5) or isempty_(ub):
        ub = options.inf * ones_(n, 1)
    else:
        ub = ub[:]
        if any_(size_(ub) != [n, 1]):
            fprintf_('\n### sqpdfo_compute_multiplier: incorrect size of ub\n\n')
            info.flag = values.fail_unexpected
            return lm, info

    # Form the matrix A

    A = concatenate_([eye_(n+nbr_slacks), gconstraints])

    # Compute the lower (lo) and upper (up) bounds on lm

    lo = - inf * ones_(n + nbr_slacks + me, 1)
    up = inf * ones_(n + nbr_slacks + me, 1)
    
    for i in arange(0, n):
        if (lb[i] <= - options.inf):
            lo[i] = 0
        if (ub[i] >= options.inf):
            up[i] = 0

        # at this point, if there are lower and upper bounds on [x,ci], 
        # the multiplier can take any sign;
        # take additional constraint in case [x,ci] is active at a bound

        if (lb[i] > - options.inf) and (abs(x[i] - lb[i]) < options.dxmin):
            up[i] = 0

        if (ub[i] < options.inf) and (abs(x[i] - ub[i]) < options.dxmin):
            lo[i] = 0
            
    for i in arange(n, n+nbr_slacks):
        lo[i] = 0
        up[i] = 0

    AA = np.dot(A,A.T)
    
    # check condition of matrix AA
    
    check_condition = 0

    if check_condition:
        cthreshold = 1e+17
        AA, badcond = sqpdfo_check_cond_(AA, cthreshold, options, nargout=2)


    #  check that matrix AA is convex (otherwise convexify)
    
    check_convex = 1

    if check_convex:
        AA = sqpdfo_check_convex_(AA, options)
        
    # compute A * g
    
    Ag = A.dot(glocal)
    #Ag = numpy.dot(A,glocal)
    
    AAn = copy(AA)
    Agn = copy(Ag)
    lon = copy(lo)
    upn = copy(up)
    ifree = ones_(size_(lo))
    k = 0
    
    for i in arange(0, length_(lo)):
        if lo[i] == up[i]:
            AAn = np.delete(AAn, k, 0)
            AAn = np.delete(AAn, k, 1)
            Agn = np.delete(Agn, k, 0)
            lon = np.delete(lon, k, 0)
            upn = np.delete(upn, k, 0)
            ifree[i] = 0
        else:
            k = k + 1
            
    if not isempty_(ifree[ifree > 0]):
        sn, rn, op, exitc = blls_(AAn, - Agn, lon, upn, nargout=4)
        I = eye_(length_(lo))
        lm = np.delete(I, find_(ifree <= 0), 1).dot(sn)
    else:
        lm = zeros_(size_(lo))
    
    if nbr_slacks:    
        lm = concatenate((lm[0:n],lm[n+nbr_slacks:n+nbr_slacks+me]))
    
    return lm, info
