# -*- coding: utf-8 -*-

import helper
import numpy as np
from ecdfo_main import ecdfo_main_
from ecdfo_global_variables import get_prob
from ecdfo_prelim import ecdfo_prelim_
from ecdfo_finish import ecdfo_finish_
from runtime import *
import time

from copy import copy
from numpy import array


def ecdfo_(func=None,x0=None,lm0=None,lb=None,ub=None,options_=None,*args,**kwargs):
    """
    # [x,lm,info] = ecdfo (func,x0);
# [x,lm,info] = ecdfo (func,x0,lm0);
# [x,lm,info] = ecdfo (func,x0,lm0,lb);
# [x,lm,info] = ecdfo (func,x0,lm0,lb,ub);
# [x,lm,info] = ecdfo (func,x0,lm0,lb,ub,options);
#
# The function computes the solution to a smooth nonlinear optimization
# problem possibly having bound constraints on the n variables x(1:n),
# nonlinear inequality constraints (ci) and nonlinear equality constraints
# (ce). The problem is supposed to be written in the form:
#
#   minimize    f(x)
#   subject to  lb(1:n)      <= x     <= ub(1:n)
#               lb(n+1:n+mi) <= ci(x) <= ub(n+1:n+mi)
#               ce(x) = 0
#
# The constraints functions are ci: Rn -> Rmi and ce: Rn -> Rme
# which means that there are mi (>=0) inequality constraints and 
# me (>=0) equality constraints. Lower and upper bounds can be infinite 
# (see the meaning of options.inf below). 
#
# On entry in 'ecdfo'
#   func: function handle of the user-supplied simulator; if this one
#     is named mysimul, use the handle @mysimul as the first argument
#   x0: n vector giving the initial guess of the solution, the number
#     of variables is assumed to be the length of the given x0
#   lm0 (optional): n+mi+me vector giving the initial guess of the
#     dual solution (Lagrange or KKT multiplier),
#     - lm0(1:n) is associated with the bound constraints on x
#     - lm0(n+1:n+mi) is associated with the inequality constraints
#     - lm0(n+mi+1:n+mi+me) is associated with the equality constraints
#     the default value is the least-squares multiplier; the dimensions
#     mi, and me are known after the first call to the simulator.
#   lb (optional): n+mi vector giving the lower bound on x (first n
#     components) and ci(x), it can have infinite components (default
#     -inf)
#   ub (optional): n+mi vector giving the upper bound on x (first n
#     components) and ci(x), it can have infinite components (default
#     inf)
#   options (optional): structure for tuning the behavior of 'ecdfo'
#     (when a string is required, both lower or uppercase letters can
#     be used, multiple white spaces are considered as a single white
#     space)
#     - options.algo_method = local algorithm to use:
#       . 'quasi-Newton' (default): requires a quasi-Newton algorithm;
#         only first order derivatives will be computed; the
#         full Hessian of the Lagrangian is approximated
#       . 'Newton' requires a Newton algorithm; second order
#         derivatives will be computed in the form of the
#         Hessian of the Lagrangian 
#     - options.algo_globalization specifies the type of globalization
#       technique to use:
#       . 'unit stepsize' prevents ecdfo from using a globalization
#         technique,
#       . 'linesearch' requires ecdfo to force convergence with
#         linesearch (default)
#       . 'trust-region' requires ecdfo to force convergence using a
#         composite step trust-region methodology
#     - options.dxmin = positive number that specifes the precision to
#       which the primal variables must be determined; if the solver
#       needs to make a step smaller than 'dxmin' in the infinity-norm
#       to progress to optimality, it will stop; a too small value will
#       force the solver to work for nothing at the very end when
#       rounding errors prevent making any progress (default 1.e-20)
#     - options.fout = FID for the printed output file (default 1,
#       i.e., the screen)
#     - options.inf = a lower bound lb(i) <= -options.inf will be
#       considered to be absent and an upper bound ub(i) >= options.inf
#       will be considered to be obsent (default = inf)
#     - options.miter = maximum number of iterations (default = 1000)
#     - options.tol = tolerance on optimality
#       (1) = tolerance on the Lagrangian gradient
#       (2) = tolerance on feasibility
#       (3) = tolerance on complementarity
#     - options.verbose = verbosity level for the outputs (default 1)
#        = 0: nothing is printed
#       >= 1: error messages (default)
#       >= 2: initial setting and final status
#       >= 3: one line per iteration
#       >= 4: details on the iterations
#       >= 5: details on the globalization
#       >= 6: some additional information, requiring expensive
#             operations, such as the computation of eigenvalues
#
# On return from 'ecdfo'
#   x: n vector giving the computed primal solution
#   lm: n+mi+me vector giving the dual solution (Lagrange or KKT
#     multiplier),
#     - lm(1:n) is associated with the bound constraints on x
#     - lm(n+1:n+mi) is associated with the inequality constraints
#     - lm(n+mi+1:n+mi+me) is associated with the equality constraints
#   info.flag: output status
#      0: solution found
#      1: an argument is wrong
#      2: unaccepted problem structure
#      3: error on a simulation
#      4: max iterations
#      5: max simulations
#      6: stop required by the simulator
#      7: stop on dxmin
#      8: infeasible QP
#      9: unbounded QP
#     10: nondescent direction (in linesearch)
#     99: should not have happened, call your guru
#   info.niter: number of realized iterations
#   info.nsimul(i): number of simulations with indic = i
#
# To have information on the problem, 'ecdfo' uses a simulator (direct
# communication), whose name is given by the argument 'func'. If
# 'mysimul' is the name of the simulator, the procedure is called by
# ecdfo by one of the following manner
#
#   [outdic] = mysimul (indic,x,lm);                 # indic=1
#   [outdic,f,ci,ce,cs,g,ai,ae] = mysimul (indic,x); # indic=2:4
#   [outdic,hl] = mysimul (indic,x,lm);              # indic=5:6
#   [outdic,hlv] = mysimul (indic,x,lm,v);           # indic=7
#   [outdic,mv] = mysimul (indic,x,v);               # indic=11:16
#
# where the input and output arguments of 'mysimul' have the following
# meaning:
#   indic is used by 'ecdfo' to drive the behavior of the simulator
#      1: the simulator will do whatever (the simulator is called with
#         indic = 1 at each iteration, so that these calls can be used
#         to count the iterations within the simulator)
#      2: the simulator will compute f(x), ci(x), ce(x), and cs(x)
#      3: the simulator will compute g(x), ai(x), and ae(x), and
#         prepare for subsequent evaluations of as(x)*v, asm(x)*v,
#         zsm(x)*v, etc
#      4: the simulator will compute f(x), ci(x), ce(x), cs(x), g(x),
#         ai(x), and ae(x), and prepare for subsequent evaluations of
#         as(x)*v, asm(x)*v, zsm(x)*v, etc
#      5: the simulator will compute hl(x,lm)
#     11: the simulator will compute as(x)*v, where as(x) is the
#         Jacobian of the state constraint cs at x and v is an n-vector
#     12: the simulator will compute as(x)'*v, where as(x) is the
#         Jacobian of the state constraint cs at x and v is an n-vector
#     13: the simulator will compute asm(x)*v, where asm(x) is a right
#         inverse of the Jacobian of the state constraint cs at x and v
#         is an ms-vector
#     14: the simulator will compute asm(x)'*v, where asm(x) is a right
#         inverse of the Jacobian of the state constraint cs at x and v
#         is an n-vector
#     15: the simulator will compute zsm(x)*v, where zsm(x) is a matrix
#         whose columns form a basis of the null space of as(x) and v
#         is an (n-ms)-vector
#     16: the simulator will compute zsm(x)'*v, where zsm(x) is a
#         matrix whose columns form a basis of the null space of as(x)
#         and v is an n-vector
#   x: n vector of variables at which the functions have to be evaluated
#   lm: KKT (or Lagrange) multiplier associated with the constraints,
#     at which the Hessian of the Lagrangian have to be computed,
#     - lm(1:n) is associated with the bound constraints on x
#     - lm(n+1:n+mi) is associated with the inequality constraints
#     - lm(n+mi+1:n+mi+me) is associated with the equality constraints
#   v: a vector that is aimed to multiply one of the matrices as(x),
#     asm(x), asm(x)', zsm(x), and zsm(x)'; its dimension depends on
#     the number of column of the matrix, hence on the value of indic
#     (see above)
#   outdic is supposed to describe the result of the simulation
#      0: the required computation has been done
#      1: x is out of an implicit domain
#      2: the simulator wants to stop
#      3: incorrect input parameters
#   f: cost-function value at x
#   ci: mi-vector giving the inequality constraint value at x
#   ce: me-vector giving the equality constraint value at x
#   cs: ms-vector giving the state constraint value at x
#   g: n vector giving the gradient of f at x
#   ai: matrix giving the Jacobian of the inequality constraints at x
#   ae: matrix giving the Jacobian of the equality constraints at x
#   hl: n x n matrix giving the Hessian of the Lagrangian at (x,lm)
#   mv: is one of the product as(x)*v, asm(x)*v, or zsm(x)*v, depending
#     on the value of indic
#
#-----------------------------------------------------------------------

# Authors: Jean Charles Gilbert, INRIA.
#      and Anke Troeltzsch, DLR.
#
# Copyright 2008, 2009, INRIA. 2013, DLR.
#
# EC-DFO is distributed under the terms of the Q Public License version
# 1.0.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the Q Public
# License version 1.0 for more details.
#
# You should have received a copy of the Q Public License version 1.0
# along with this program.  If not, see
# <http://doc.trolltech.com/3.0/license.html>.
#
    """
#    varargin = cellarray(args)

    options=copy(options_)
    nones = [func is None,x0 is None,lm0 is None,lb is None,ub is None,options is None].count(True)
    nargin = 6-nones+len(args)
#  Set some constants

    c = helper.dummyUnionStruct()
    info = helper.dummyUnionStruct()
    eps = 2.220446049250313e-16
    prob=get_prob()
    c.free=0
    c.fixed=1
    c.alwaysfixed=2
    c.in_=1
    c.out=0
    c.unused=0
    c.inY=1
    c.dummy=1
    c.nodummy=0
#  Initialize to zero
    x=copy(np.NaN)
    fx=copy(np.NaN)
    gx=copy(np.NaN)
    nit=0
    nitold=0
    neval=0
    errg=copy(np.Inf)
    X=array([])
    fX=array([])
    xstatus=array([])
    sstatus=array([])
    dstatus=array([])
    sspace_save=array([])
    xspace_save=array([])
    ndummyY=0
    info.flag=0

    #Reset the random generator : did not find something equivalent in python for rand_('seed',5) but this is commented in matlab
#    random.seed(np.pi / np.sqrt(2))
#    randn_('seed',5)
    #  Miscellaneous initializations
  
    if ( size_(x0,1) == 1 and size_(x0,2) > 1 ):    # make sure the starting point is 
       x0 = x0.T                                # a column vector
    
    n       = length_( x0 )                     # the dimension of the space
    #pquad   = int(( ( n + 1 ) * ( n + 2 ) ) / 2)    # the size of a fully quadratic model
    pquad   = int( 2* n )
    pdiag   = int(2 * n + 1)                        # the size of a diagonal model
    plin    = int(n + 1)                             # the size of a linear model
    
    msg              = 'Unexpected exit'       # no meaningful message at this point
    poisedness_known = 0                       # the poisedness is not known yet
    eps_rho          = 1.0e-14                 # stabilization constant for rho
    stallfact        = 10 * eps                # termination (stalled) when
                                                # ||s|| <= stallfact * norm( x )
    factor_Dmax      = 1e5                    # the ratio between Deltamax and Delta0
    factor_fmax      = 1e20                   # the ratio between the upper bound on the objective function
                                                # value and the absolute value of t(f(x0))
    CNTsin           = 0                       # variable to produce a quasi-random vector

    #  Set defaults
    
    Delta0       = 1                # the initial TR radius
    cur_degree   = copy(plin)             # the degree for the initial model
    rep_degree   = copy(plin)             # the minimum degree for the model after repair
    epsilon      = 1.0e-5           # gradient termination accuracy
    maxeval      = 200 * n          # maximum number of evaluations
    maxit        = copy(maxeval)          # maximum number of iterations
    verbose      = 1                # printout quantity
    show_errg    = 0                # display of the gradient error estimate
    initial_Y    = 'simplx'         # geometry of the initial interpolation set
    eta1         = 0.0001           # min rho value for successful iterations
    eta2         = 0.9              # min rho value for very successful iterations
    gamma1       = 0.01             # lower bound on Delta interval (unsuccessful its)
    gamma2       = 0.5              # upper bound on Delta interval (unsuccessful its)
    gamma3       = 2.0              # increase in Delta (successful its)
    interpol_TR  = 1                # use interpolation for the radius of the trust region
    factor_CV    = 100               # constant for termination (see above)
    Lambda_XN    = 1.0e-10          # poisedness for new iterates
    Lambda_CP    = 1.2              # poisedness for close points
    Lambda_FP    = 1.0e-10          # poisedness for far points
    factor_FPU   = 1                # multiple of TR radius defining far points
                                     # (for unsuccessful iterations)
    factor_FPR   = 10               # multiple of TR radius defining far points
                                     # (for reconstruction of a poised interpolation set
                                     # for large model gradient)
    criterion_S  = 'distance'       # selection of outgoing point at successful iterations:
    criterion_FP = 'distance'       # the same, but for far points at unsuccessful iterations
    criterion_CP = 'standard'       # the same, but for close points at unsuccessful iterations
    mu0          = 0                # initial reduction in gradient norm before repair
    mu           = 0                # subsequent reduction in gradient norm before repair
    theta        = 1                # ratio between gradient norm and radius after repair
    eps_TR       = 0.0001           # rel. accuracy on the trust-region constraint for steps
    eps_L        = 0.001            # rel. accuracy on the trust-region constraint for L max
    shift_Y      = 1                # shifting and scaling of the set Y is used
    lSolver      = 1                # local solver for minimizing the model. 1: 2-norm, 2: inf-norm.
    stratLam     = 1                # strategy to adjust lambda when solving the bc MS problem,
    kappa_ill    = 1e+15            # threshold to declare a system matrix as ill-conditioned
    kappa_th     = 2000             # threshold for a safely nondegenerate set of points
    eps_bnd      = epsilon/10       # epsilon to define a bound as nearly-active: |x - bnd|<eps_bnd
    whichmodel   = 2                # approach to build the local models
    hardcons     = 0                # apply hard bounds when maximizing the Lagrange polynomials
    noisy        = 0                # function supposed to be noisy
    scaleX       = 0                # scaling of variables is applied
    scalefacX    = ones_(1,n)        # scaling factors initialized to one
    shrink_Delta = 1                # shrink trust-region radius in every unsuccessful iteration
    
    #  Compute the maximal TR radius.
    
    Deltamax = factor_Dmax * Delta0
    
    # -------------------------------------------------------------------------
    # Check input arguments
    # -------------------------------------------------------------------------

    if nargin < 2:
        fprintf_('\n### EC-DFO: the first 2 arguments are required\n\n')
        x=array([])
        lm=array([])
        info.flag=1
        return x,lm,info
    if nargin < 3:
        lm0=array([])
    else:
        lm0=lm0[:]
    if nargin < 4:
        lb=array([])
    else:
        lb=lb[:]
    if nargin < 5:
        ub=array([])
    else:
        ub=ub[:]
    if nargin < 6:
        options.fout=1
        options.verbose=1
    # GTlab problem
    #if prob == 100:
        #Delta0=0.01
        #epsilon=0.001
        #scaleX = 1;
        #scalefacX = array([[100,0.01, 1, 0.1, 0.1]]).T;

    # -------------------------------------------------------------------------
    # Preliminaries:
    # -------------------------------------------------------------------------

    # - check bounds and the given options
    # - build initial poised interpolation set
    # - compute function + constraint values and initial multiplier (if not given)
    n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta,nfix,indfix,xfix,vstatus,xstatus,sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,poised_model,X,fX,Y,fY,ciX,ciY,ceX,ceY,poisedness_known,m,gx,normgx,fcmodel,ind_Y,i_xbest,cur_degree,rep_degree,plin,pdiag,pquad,indfree,info,options,values=ecdfo_prelim_(func,x0,lm0,Delta0,lb,ub,scaleX,scalefacX,cur_degree,rep_degree,plin,pdiag,pquad,c,initial_Y,kappa_ill,whichmodel,factor_FPR,Lambda_FP,Lambda_CP,eps_L,lSolver,hardcons,stratLam,xstatus,sstatus,dstatus,options,nargout=47)
    if info.flag:
        return x,lm,info
    # Initialize the current epsilon

    eps_current=max_(mu0 * normgx,epsilon)
    # Compute the maximal objective function value.

    fxmax=min_(1e+25,factor_fmax * abs(fx))
    # Initialize the Hessian of the Lagrangian or its approximation

    M=eye_(n)
    # printout to file convhist.m

    if (verbose):
        fid=fopen_('convhist.m','w')
        fprintf_(fid,'function A=history \n A=[ \n')
        fclose_(fid)
    # -------------------------------------------------------------------------
    # Call main optimization loop
    # -------------------------------------------------------------------------


    nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,Delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info=ecdfo_main_(func,n,nb,mi,me,lm,nitold,nit,i_xbest,lb,ub,m,X,fX,ciX,ceX,ind_Y,QZ,RZ,Delta,cur_degree,neval,maxeval,maxit,fcmodel,gx,normgx,show_errg,pquad,pdiag,plin,stallfact,eps_rho,Deltamax,rep_degree,epsilon,verbose,eta1,eta2,gamma1,gamma2,gamma3,interpol_TR,factor_CV,Lambda_XN,Lambda_CP,factor_FPU,factor_FPR,Lambda_FP,criterion_S,criterion_FP,criterion_CP,mu,theta,eps_TR,eps_L,lSolver,stratLam,eps_current,vstatus,xstatus,sstatus.T,dstatus,ndummyY,sspace_save,xspace_save,xfix,fxmax,poised_model,M,kappa_ill,kappa_th,eps_bnd,poised,Y_radius,c,'toplevel',whichmodel,hardcons,noisy,scaleX,scalefacX,CNTsin,shrink_Delta,scale,shift_Y,info,options,values,nargout=29)

    #  Append closing bracket in file convhist.m

    if (verbose):
        fid=fopen_('convhist.m','a')
        fprintf_(fid,'];')
        fclose_(fid)
    #  Re-assemble gradient at return

    if (nfix > 0):
        I=eye_(n + nfix)
        x=I[:,indfix] * xfix[indfix] + I[:,indfree].dot(x)
        gx=I[:,indfix] * zeros_(nfix,1) + I[:,indfree].dot(gx)
        Ilm=eye_(n + nfix + me + mi)
        indfree_lm=setdiff_(arange(0,n + nfix + me + mi),indfix)
        lm=Ilm[:,indfix] * zeros_(nfix,1) + Ilm[:,indfree_lm] * lm
        n=n + nfix
    #  Rescale x if necessary

    if (scaleX):
        x=x / scalefacX
    info_best=copy(info)
    info_best.f=fx
    x=X[:,i_xbest]
    ecdfo_finish_(nb,mi,me,info_best,options,values)
#    disp_('  Best result is function value number : ',str(i_xbest),' with fval = ', str(fX[i_xbest]))
#    if me:
#       fprintf_(1,'%s','  with nce ');
#       fprintf_(1,'%g ',ceX[:,i_xbest]);
#    fprintf_(1,'\n');
    
    #---------------------------------------
    # Recover the results
    #---------------------------------------
    if options.verbose > 2:
        if nb:
            fprintf_(options.fout,'VARIABLES:\n')
            fprintf_(options.fout,'i     lower bound          x            upper bound       multiplier\n')
            for i in range(0,min_(n,40)):
                fprintf_(options.fout,'%0i %+16.6e %+16.6e %+16.6e %+16.6e\n'%(i,lb[i],x[i],ub[i],lm[i]))
            if (n > 40):
                fprintf_(options.fout,'.....\n')
            else:
                fprintf_(options.fout,'\n')
        else:
            fprintf_(options.fout,'VARIABLES:\n')
            fprintf_(options.fout,'\n',x[0:min_(n,40)])
            if (n > 40):
                fprintf_(options.fout,'.....\n')
            else:
                fprintf_(options.fout,'\n')
        if mi:
            fprintf_(options.fout,'INEQUALITY CONSTRAINTS:\n')
            fprintf_(options.fout,'i     lower bound          ci           upper bound       multiplier\n')
            for i in range(0,min_(mi,40)):
                fprintf_(options.fout,'%0i %+16.6e %+16.6e %+16.6e %+16.6e\n'%(i,lb[n + i],ciX[i,i_xbest],ub[n + i],lm[n + i]))
            if (mi > 40):
                fprintf_(options.fout,'\n.....')
            else:
                fprintf_(options.fout,'\n')
        if me:
            fprintf_(options.fout,'EQUALITY CONSTRAINTS:\n')
            fprintf_(options.fout,'i         ce            multiplier\n')
            for i in range(0,min_(me,40)):
                fprintf_(options.fout,'%0i %+16.6e %+16.6e\n'%(i,ceX[i,i_xbest],lm[n + mi + i]))
            if (me > 40):
                fprintf_(options.fout,'.....\n')
            else:
                fprintf_(options.fout,'\n')
    return x,lm,info
				
