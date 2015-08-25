# -*- coding: utf-8 -*-
from __future__ import division
import helper
import numpy as np
from ecdfo_main import ecdfo_main_
from ecdfo_global_variables import get_prob
from ecdfo_prelim import ecdfo_prelim_
from ecdfo_finish import ecdfo_finish_
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from copy import copy
from numpy import array



def ecdfo_(func_=None,x0_=None,lm0_=None,lb_=None,ub_=None,options_=None,*args,**kwargs):
    """
    % [x,lm,info] = ecdfo (func,x0);
% [x,lm,info] = ecdfo (func,x0,lm0);
% [x,lm,info] = ecdfo (func,x0,lm0,lb);
% [x,lm,info] = ecdfo (func,x0,lm0,lb,ub);
% [x,lm,info] = ecdfo (func,x0,lm0,lb,ub,options);
%
% The function computes the solution to a smooth nonlinear optimization
% problem possibly having bound constraints on the n variables x(1:n),
% nonlinear inequality constraints (ci) and nonlinear equality constraints
% (ce). The problem is supposed to be written in the form:
%
%   minimize    f(x)
%   subject to  lb(1:n)      <= x     <= ub(1:n)
%               lb(n+1:n+mi) <= ci(x) <= ub(n+1:n+mi)
%               ce(x) = 0
%
% The constraints functions are ci: Rn -> Rmi and ce: Rn -> Rme
% which means that there are mi (>=0) inequality constraints and 
% me (>=0) equality constraints. Lower and upper bounds can be infinite 
% (see the meaning of options.inf below). 
%
% On entry in 'ecdfo'
%   func: function handle of the user-supplied simulator; if this one
%     is named mysimul, use the handle @mysimul as the first argument
%   x0: n vector giving the initial guess of the solution, the number
%     of variables is assumed to be the length of the given x0
%   lm0 (optional): n+mi+me vector giving the initial guess of the
%     dual solution (Lagrange or KKT multiplier),
%     - lm0(1:n) is associated with the bound constraints on x
%     - lm0(n+1:n+mi) is associated with the inequality constraints
%     - lm0(n+mi+1:n+mi+me) is associated with the equality constraints
%     the default value is the least-squares multiplier; the dimensions
%     mi, and me are known after the first call to the simulator.
%   lb (optional): n+mi vector giving the lower bound on x (first n
%     components) and ci(x), it can have infinite components (default
%     -inf)
%   ub (optional): n+mi vector giving the upper bound on x (first n
%     components) and ci(x), it can have infinite components (default
%     inf)
%   options (optional): structure for tuning the behavior of 'ecdfo'
%     (when a string is required, both lower or uppercase letters can
%     be used, multiple white spaces are considered as a single white
%     space)
%     - options.algo_method = local algorithm to use:
%       . 'quasi-Newton' (default): requires a quasi-Newton algorithm;
%         only first order derivatives will be computed; the
%         full Hessian of the Lagrangian is approximated
%       . 'Newton' requires a Newton algorithm; second order
%         derivatives will be computed in the form of the
%         Hessian of the Lagrangian 
%     - options.algo_globalization specifies the type of globalization
%       technique to use:
%       . 'unit stepsize' prevents ecdfo from using a globalization
%         technique,
%       . 'linesearch' requires ecdfo to force convergence with
%         linesearch (default)
%       . 'trust-region' requires ecdfo to force convergence using a
%         composite step trust-region methodology
%     - options.dxmin = positive number that specifes the precision to
%       which the primal variables must be determined; if the solver
%       needs to make a step smaller than 'dxmin' in the infinity-norm
%       to progress to optimality, it will stop; a too small value will
%       force the solver to work for nothing at the very end when
%       rounding errors prevent making any progress (default 1.e-20)
%     - options.fout = FID for the printed output file (default 1,
%       i.e., the screen)
%     - options.inf = a lower bound lb(i) <= -options.inf will be
%       considered to be absent and an upper bound ub(i) >= options.inf
%       will be considered to be obsent (default = inf)
%     - options.miter = maximum number of iterations (default = 1000)
%     - options.tol = tolerance on optimality
%       (1) = tolerance on the Lagrangian gradient
%       (2) = tolerance on feasibility
%       (3) = tolerance on complementarity
%     - options.verbose = verbosity level for the outputs (default 1)
%        = 0: nothing is printed
%       >= 1: error messages (default)
%       >= 2: initial setting and final status
%       >= 3: one line per iteration
%       >= 4: details on the iterations
%       >= 5: details on the globalization
%       >= 6: some additional information, requiring expensive
%             operations, such as the computation of eigenvalues
%
% On return from 'ecdfo'
%   x: n vector giving the computed primal solution
%   lm: n+mi+me vector giving the dual solution (Lagrange or KKT
%     multiplier),
%     - lm(1:n) is associated with the bound constraints on x
%     - lm(n+1:n+mi) is associated with the inequality constraints
%     - lm(n+mi+1:n+mi+me) is associated with the equality constraints
%   info.flag: output status
%      0: solution found
%      1: an argument is wrong
%      2: unaccepted problem structure
%      3: error on a simulation
%      4: max iterations
%      5: max simulations
%      6: stop required by the simulator
%      7: stop on dxmin
%      8: infeasible QP
%      9: unbounded QP
%     10: nondescent direction (in linesearch)
%     99: should not have happened, call your guru
%   info.niter: number of realized iterations
%   info.nsimul(i): number of simulations with indic = i
%
% To have information on the problem, 'ecdfo' uses a simulator (direct
% communication), whose name is given by the argument 'func'. If
% 'mysimul' is the name of the simulator, the procedure is called by
% ecdfo by one of the following manner
%
%   [outdic] = mysimul (indic,x,lm);                 % indic=1
%   [outdic,f,ci,ce,cs,g,ai,ae] = mysimul (indic,x); % indic=2:4
%   [outdic,hl] = mysimul (indic,x,lm);              % indic=5:6
%   [outdic,hlv] = mysimul (indic,x,lm,v);           % indic=7
%   [outdic,mv] = mysimul (indic,x,v);               % indic=11:16
%
% where the input and output arguments of 'mysimul' have the following
% meaning:
%   indic is used by 'ecdfo' to drive the behavior of the simulator
%      1: the simulator will do whatever (the simulator is called with
%         indic = 1 at each iteration, so that these calls can be used
%         to count the iterations within the simulator)
%      2: the simulator will compute f(x), ci(x), ce(x), and cs(x)
%      3: the simulator will compute g(x), ai(x), and ae(x), and
%         prepare for subsequent evaluations of as(x)*v, asm(x)*v,
%         zsm(x)*v, etc
%      4: the simulator will compute f(x), ci(x), ce(x), cs(x), g(x),
%         ai(x), and ae(x), and prepare for subsequent evaluations of
%         as(x)*v, asm(x)*v, zsm(x)*v, etc
%      5: the simulator will compute hl(x,lm)
%     11: the simulator will compute as(x)*v, where as(x) is the
%         Jacobian of the state constraint cs at x and v is an n-vector
%     12: the simulator will compute as(x)'*v, where as(x) is the
%         Jacobian of the state constraint cs at x and v is an n-vector
%     13: the simulator will compute asm(x)*v, where asm(x) is a right
%         inverse of the Jacobian of the state constraint cs at x and v
%         is an ms-vector
%     14: the simulator will compute asm(x)'*v, where asm(x) is a right
%         inverse of the Jacobian of the state constraint cs at x and v
%         is an n-vector
%     15: the simulator will compute zsm(x)*v, where zsm(x) is a matrix
%         whose columns form a basis of the null space of as(x) and v
%         is an (n-ms)-vector
%     16: the simulator will compute zsm(x)'*v, where zsm(x) is a
%         matrix whose columns form a basis of the null space of as(x)
%         and v is an n-vector
%   x: n vector of variables at which the functions have to be evaluated
%   lm: KKT (or Lagrange) multiplier associated with the constraints,
%     at which the Hessian of the Lagrangian have to be computed,
%     - lm(1:n) is associated with the bound constraints on x
%     - lm(n+1:n+mi) is associated with the inequality constraints
%     - lm(n+mi+1:n+mi+me) is associated with the equality constraints
%   v: a vector that is aimed to multiply one of the matrices as(x),
%     asm(x), asm(x)', zsm(x), and zsm(x)'; its dimension depends on
%     the number of column of the matrix, hence on the value of indic
%     (see above)
%   outdic is supposed to describe the result of the simulation
%      0: the required computation has been done
%      1: x is out of an implicit domain
%      2: the simulator wants to stop
%      3: incorrect input parameters
%   f: cost-function value at x
%   ci: mi-vector giving the inequality constraint value at x
%   ce: me-vector giving the equality constraint value at x
%   cs: ms-vector giving the state constraint value at x
%   g: n vector giving the gradient of f at x
%   ai: matrix giving the Jacobian of the inequality constraints at x
%   ae: matrix giving the Jacobian of the equality constraints at x
%   hl: n x n matrix giving the Hessian of the Lagrangian at (x,lm)
%   mv: is one of the product as(x)*v, asm(x)*v, or zsm(x)*v, depending
%     on the value of indic
%
%-----------------------------------------------------------------------

% Authors: Jean Charles Gilbert, INRIA.
%      and Anke Troeltzsch, DLR.
%
% Copyright 2008, 2009, INRIA. 2013, DLR.
%
% EC-DFO is distributed under the terms of the Q Public License version
% 1.0.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the Q Public
% License version 1.0 for more details.
%
% You should have received a copy of the Q Public License version 1.0
% along with this program.  If not, see
% <http://doc.trolltech.com/3.0/license.html>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    """
#    varargin = cellarray(args)

    func=copy(func_)
    x0=copy(x0_)
    lm0=copy(lm0_)
    lb=copy(lb_)
    ub=copy(ub_)
    options=copy(options_)

    nargin = 6-[func,x0,lm0,lb,ub,options].count(None)+len(args)

    c = helper.dummyUnionStruct()
    info = helper.dummyUnionStruct()
    eps = 2.220446049250313e-16
    prob=get_prob
    c.free=0
    c.fixed=1
    c.alwaysfixed=2
    c.in_=1
    c.out=0
    c.unused=0
    c.inY=1
    c.dummy=1
    c.nodummy=0
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
#    rand_('seed',np.pi / sqrt_(2))
#    randn_('seed',5)
    if (size_(x0,1) == 1 and size_(x0,2) > 1):
        x0=x0.T
    n=length_(x0)
    pquad=((n + 1) * (n + 2)) / 2
    pdiag=2 * n + 1
    plin=n + 1
    msg='Unexpected exit'
    poisedness_known=0
    eps_rho=1e-14
    stallfact=10 * eps
    factor_Dmax=100000.0
    factor_fmax=1e+20
    CNTsin=0
    Delta0=1
    cur_degree=copy(plin)
    rep_degree=copy(plin)
    epsilon=1e-05
    maxeval=200 * n
    maxit=copy(maxeval)
    verbose=options.verbose
    show_errg=0
    initial_Y='simplx'
    eta1=0.0001
    eta2=0.9
    gamma1=0.01
    gamma2=0.5
    gamma3=2.0
    interpol_TR=1
    factor_CV=100
    Lambda_XN=1e-10
    Lambda_CP=1.2
    Lambda_FP=1e-10
    factor_FPU=1
    factor_FPR=10
    criterion_S='distance'
    criterion_FP='distance'
    criterion_CP='standard'
    mu0=0
    mu=0
    theta=1
    eps_TR=0.0001
    eps_L=0.001
    shift_Y=1
    lSolver=1
    stratLam=1
    kappa_ill=1e+15
    kappa_th=2000
    eps_bnd=epsilon / 10
    whichmodel=0
    hardcons=0
    noisy=0
    scaleX=0
    scalefacX=ones_(1,n)
    shrink_Delta=1
    Deltamax=factor_Dmax * Delta0
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
    if prob == 100:
        Delta0=0.01
        epsilon=0.001
        if exist_('fvalues_ecdfo_karmanogive.dat','file') == 2: #Those 2 lines (exist_ and delete_) have not been dealt with, but only concerns prob==100
            delete_('fvalues_ecdfo_karmanogive.dat')
    print 'x0 = ',x0
    n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta,nfix,indfix,xfix,vstatus,xstatus,sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,poised_model,X,fX,Y,fY,ciX,ciY,ceX,ceY,poisedness_known,m,gx,normgx,fcmodel,ind_Y,i_xbest,cur_degree,rep_degree,plin,pdiag,pquad,indfree,info,options,values=ecdfo_prelim_(func,x0,lm0,Delta0,lb,ub,scaleX,scalefacX,cur_degree,rep_degree,plin,pdiag,pquad,c,initial_Y,kappa_ill,whichmodel,factor_FPR,Lambda_FP,Lambda_CP,eps_L,lSolver,hardcons,stratLam,xstatus,sstatus,dstatus,options,nargout=47)
    if info.flag:
        return x,lm,info
    x0=copy(x)
    eps_current=max_(mu0 * normgx,epsilon)
    fxmax=min_(1e+25,factor_fmax * abs(fx))
    M=eye_(n)
    if (verbose):
        fid=fopen_('convhist.m','w')
        fprintf_(fid,'function A=history \n A=[ \n')
        fprintf_(fid,'%6d  %+.14e %.2e \n'%(neval,fx,normgx))
        fclose_(fid)
    nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,Delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info=ecdfo_main_(func,n,nb,mi,me,lm,nitold,nit,i_xbest,lb,ub,m,X,fX,ciX,ceX,ind_Y,QZ,RZ,Delta,cur_degree,neval,maxeval,maxit,fcmodel,gx,normgx,show_errg,pquad,pdiag,plin,stallfact,eps_rho,Deltamax,rep_degree,epsilon,verbose,eta1,eta2,gamma1,gamma2,gamma3,interpol_TR,factor_CV,Lambda_XN,Lambda_CP,factor_FPU,factor_FPR,Lambda_FP,criterion_S,criterion_FP,criterion_CP,mu,theta,eps_TR,eps_L,lSolver,stratLam,eps_current,vstatus,xstatus,sstatus.T,dstatus,ndummyY,sspace_save,xspace_save,xfix,fxmax,poised_model,M,kappa_ill,kappa_th,eps_bnd,poised,Y_radius,c,'toplevel',whichmodel,hardcons,noisy,scaleX,scalefacX,CNTsin,shrink_Delta,scale,shift_Y,info,options,values,nargout=29)
    if (verbose):
        fid=fopen_('convhist.m','a')
        fprintf_(fid,'];')
        fclose_(fid)
    if (nfix > 0):
        I=eye_(n + nfix)
        x=I[:,indfix] * xfix[indfix] + I[:,indfree].dot(x)
        gx=I[:,indfix] * zeros_(nfix,1) + I[:,indfree].dot(gx)
        Ilm=eye_(n + nfix + me + mi)
        indfree_lm=setdiff_(arange(0,n + nfix + me + mi),indfix)
        lm=Ilm[:,indfix] * zeros_(nfix,1) + Ilm[:,indfree_lm] * lm
        n=n + nfix
    if (scaleX):
        x=x / scalefacX
    info_best=copy(info)
    info_best.f=fx
    x=X[:,i_xbest]
    ecdfo_finish_(nb,mi,me,info_best,options,values)
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
				
