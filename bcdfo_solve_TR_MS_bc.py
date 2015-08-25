# -*- coding: utf-8 -*-
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from copy import copy
from numpy import *

from bcdfo_solve_TR_MS import bcdfo_solve_TR_MS_

def bcdfo_solve_TR_MS_bc_(g_=None,H=None,lb=None,ub=None,Delta_=None,eps_D=None,stratLam=None,*args,**kwargs):
    """
#
#  An implementation of exact trust-region minimization based on the
#  More-Sorensen algorithm subject to bound constraints.
#
#  INPUT:
#
#  g        : the model's gradient
#  H        : the model's Hessian
#  lb       : lower bounds on the step
#  ub       : upper bounds on the step
#  Delta    : the trust-region's radius
#  eps_D    : the accuracy required on the equation ||s|| = Delta for a
#             boundary solution
#  stratLam : the strategy to adjust lamb to find an active bound
#             (1 - Newtons method, 0 - bisection method)
#
#  OUTPUT:
#
#  s        : the trust-region step
#  lamb   : the Lagrange multiplier corresponding to the trust-region constraint
#  norms    : the norm of s
#  value    : the value of the model at the optimal solution
#  gplus    : the value of the model's gradient at the optimal solution
#  nfact    : the number of Cholesky factorization required
#  neigd    : the number of eifenvalue decompositions required
#  msg      : an information message
#
#  DEPENDENCIES: bcdfo_solve_TR_MS
#
#  PROGRAMMING: A. Troeltzsch, S. Gratton, July 2009.
#              ( This version 14 I 2010 )
#
#  TEST:
#
#  bcdfo_solve_TR_MS_bc( [2; 3], [4 6; 6 5], [-10; -10], [10; 10], 1.0, 0.001, 1 )
#  should give
#    0.5153
#   -0.8575
#
#  bcdfo_solve_TR_MS_bc( [2; 3], [4 6; 6 5], [-0.1; -0.1], [10; 10], 1.0, 0.001, 1 )
#  should give
#   -0.1
#   -0.1
#
#  bcdfo_solve_TR_MS_bc( [2; 3], [4 6; 6 5], [-10; -10], [0; 0], 1.0, 0.001, 1 )
#  should give
#    0
#   -0.6
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
    """    
#    varargin = cellarray(args)
#    nargin = 7-[g,H,lb,ub,Delta,eps_D,stratLam].count(None)+len(args)

    Delta=copy(Delta_)    
    g=copy(g_)
    
    lambda_save = []
#    norms_b_save = []
    verbose=0
    theta=1e-13
    eps_bound=1e-05
    msg='no free variables'
    _lambda=0
    value=0
    nfact=0
    neigd=0
    Delta0=copy(Delta)
    g0=copy(g)
    gplus=copy(g)
    s=zeros_(size_(g))
    norms=0
    n=length_(g)
    I=eye_(n)
    ind_active=array([])
    ind_free=arange(0,n)
    nfree=copy(n)
    if (verbose):
        disp_('bcdfo_solve_TR_MS_bc: enter')
    if (not isempty_(find_(isnan(H)))):
        if (verbose):
            disp_('Error in bcdfo_solve_TR_MS_bc: H contains NaNs!')
        msg='error1'
        if (verbose):
            disp_('bcdfo_solve_TR_MS_bc: exit!')
        return s,_lambda,norms,value,gplus,nfact,neigd,msg
    if (not isempty_(find_(~ isreal(H)))):
        if (verbose):
            disp_('Error in bcdfo_solve_TR_MS_bc: H contains imaginary parts!')
        msg='error2'
        if (verbose):
            disp_('bcdfo_solve_TR_MS_bc: exit!')
        return s,_lambda,norms,value,gplus,nfact,neigd,msg
    if (not isempty_(find_(isinf(H)))):
        if (verbose):
            disp_('Error in bcdfo_solve_TR_MS_bc: H contains infinite elements!')
        msg='error3'
        if (verbose):
            disp_('bcdfo_solve_TR_MS_bc: exit!')
        return s,_lambda,norms,value,gplus,nfact,neigd,msg
#    ind_g_crit=find_((abs(lb) <= 1e-10 and g > 0) or (ub <= 1e-10 and g < 0))
    ind_g_crit=find_(logical_or(logical_and(g>0, abs(lb) <= 1e-10), logical_and(g < 0, ub <= 1e-10)))

    if (not isempty_(ind_g_crit)):
        ind_active=ind_free[ind_g_crit].T
        ind_free=setdiff_(ind_free,ind_active)
        nfree=length_(ind_free)
    j=0
    while nfree > 0:

        new_call_to_MS=1
        while (new_call_to_MS == 1):

            j=j + 1
            if (verbose >= 1):
                disp_('(',num2str_(j),') ---- minimizing in the (sub)space of ',num2str_(length_(ind_free)),' variable(s)')
            g_reduced=g[ind_free]
            H_reduced=H[ix_(ind_free,ind_free)]
            s_deltaMS,_lambda,norms_deltaMS,value_red,gplus_red,nfact_r,neigd_r,msg,hardcase=bcdfo_solve_TR_MS_(g_reduced,H_reduced,Delta,eps_D,nargout=9)
            nfact=nfact + nfact_r
            neigd=neigd + neigd_r
            gplus[ind_free]=gplus[ind_free] + gplus_red
            s_after_reduced_ms=s + I[:,ind_free] .dot( s_deltaMS )
#            ind_u_crit=find_((ub[ind_free] - s_after_reduced_ms[ind_free]) <= eps_bound and ub[ind_free] <= 1e-10)
            ind_u_crit=find_(logical_and(ub[ind_free] - s_after_reduced_ms[ind_free] <= eps_bound, ub[ind_free] <= 1e-10))
#            ind_l_crit=find_((s_after_reduced_ms[ind_free] - lb[ind_free]) <= eps_bound and lb[ind_free] >= - 1e-10)
            ind_l_crit=find_(logical_and(s_after_reduced_ms[ind_free] - lb[ind_free] <= eps_bound, lb[ind_free] >= -1e-10))
            if (length_(ind_u_crit) + length_(ind_l_crit) != 0):
#                ind_active=matlabarray([ind_active,ind_free[ind_u_crit],ind_free[ind_l_crit]])
                ind_active=concatenate_([ind_active,ind_free[ind_u_crit].T,ind_free[ind_l_crit].T],axis=1)
                ind_free=setdiff_(range(0,n),ind_active)
                nfree=length_(ind_free)
                if (verbose):
                    disp_('fixed one or more variables')
                if (nfree == 0):
                    norms=norm_(s)
                    value=0.5 * s.T.dot( H.dot( s )) + s.T .dot( g0 )
                    if (verbose):
                        disp_('no inactive variables anymore - return')
                        disp_('bcdfo_solve_TR_MS_bc: exit!')
                    return s,_lambda,norms,value,gplus,nfact,neigd,msg
            else:
                new_call_to_MS=0

        if (verbose == 2):
            disp_('check if step inside bounds')
        out_of_ubound=find_((ub[ind_free] - s_after_reduced_ms[ind_free]) < 0.0)
        out_of_lbound=find_((s_after_reduced_ms[ind_free] - lb[ind_free]) < 0.0)
        out_of_ubound_init=copy(out_of_ubound)
        out_of_lbound_init=copy(out_of_lbound)
        if (length_(out_of_ubound) + length_(out_of_lbound) != 0):
            back_inside=0
            lambda0=copy(_lambda)
            if (verbose == 2):
                disp_('step outside bounds!')
                out_of_ubound
                out_of_lbound
                disp_('lambda_0=',num2str_(lambda0))
            lower=copy(_lambda)
            if (stratLam == 0):
                _lambda=max_(2.0,2 * _lambda)
            gnorm=norm_(g)
            if (length_(out_of_ubound) > 0):
                delta_b=min_(abs(ub[ind_free[out_of_ubound]] - s[ind_free[out_of_ubound]]))
            if (length_(out_of_lbound) > 0):
                delta_b=min_(abs(lb[ind_free[out_of_lbound]] - s[ind_free[out_of_lbound]]))
                if (length_(out_of_ubound) > 0):
                    delta_b=min_(min_(abs(ub[ind_free[out_of_ubound]] - s[ind_free[out_of_ubound]])),delta_b)
            goverD=gnorm / delta_b
            Hnorminf=norm_(H,inf)
            if (Hnorminf > 0):
                HnormF=norm_(H,'fro')
            else:
                HnormF=0
            upper=max_(0,goverD + min_(Hnorminf,HnormF))
            ind_u_active=find_(abs(ub[ind_free] - s_after_reduced_ms[ind_free]) <= eps_bound)
            ind_l_active=find_(abs(s_after_reduced_ms[ind_free] - lb[ind_free]) <= eps_bound)
            i=0
            while (((length_(ind_u_active) + length_(ind_l_active)) == 0) or (length_(out_of_lbound) + length_(out_of_ubound) != 0)):

                i=i + 1
                old_lambda=copy(_lambda)
                new_lambda=- 1
                if (verbose):
                    disp_(' bcdfo_solve_TR_MS_bc (',int2str_(i),'): lower = ',num2str_(lower),' lambda = ',num2str_(_lambda),' upper = ',num2str_(upper))
                R,p=chol_(H[ix_(ind_free,ind_free)] + _lambda * eye_(nfree),nargout=2)
                if (not isempty_(find_(isnan(R)))):
                    H
                    _lambda
                    norm_(g)
                    R
                    p
                    disp_('Error in bcdfo_solve_TR_MS_bc: NaNs in Cholesky factorization')
                    msg='error4'
                    if (verbose):
                        disp_('bcdfo_solve_TR_MS_bc: exit!')
                    return s,_lambda,norms,value,gplus,nfact,neigd,msg
                nfact=nfact + 1
                if (p == 0 and hardcase == 0):
                    s_deltaH=numpy.linalg.solve(- R,(numpy.linalg.solve(R.T,g[ind_free])))
                    s_duringH=s + I[:,ind_free].dot( s_deltaH )
#                    ind_u_crit=find_((ub[ind_free] - s_duringH[ind_free]) <= eps_bound and ub[ind_free] <= 1e-10)
                    ind_u_crit=find_(logical_and(ub[ind_free] - s_duringH[ind_free] <= eps_bound, ub[ind_free] <= 1e-10))
#                    ind_l_crit=find_((s_duringH[ind_free] - lb[ind_free]) <= eps_bound and lb[ind_free] >= - 1e-10)
                    ind_l_crit=find_(logical_and(s_duringH[ind_free] - lb[ind_free] <= eps_bound, lb[ind_free] >= -1e-10))

                    if (length_(ind_u_crit) != 0):
                        s_deltaH[ind_u_crit]=0.0
                        s_duringH[ind_free[ind_u_crit]]=0.0
                    if (length_(ind_l_crit) != 0):
                        s_deltaH[ind_l_crit]=0.0
                        s_duringH[ind_free[ind_l_crit]]=0.0
                    out_of_ubound=find_((ub[ind_free] - s_duringH[ind_free]) < 0.0)
                    out_of_lbound=find_((s_duringH[ind_free] - lb[ind_free]) < 0.0)
                    if (stratLam == 1 or verbose > 0):
                        if ((length_(out_of_ubound) != 0) or (length_(out_of_lbound) != 0)):
                            outside=1
                            if (length_(out_of_ubound) != 0):
                                diff_b_u,ind_b_u=max_(abs(ub[ind_free[out_of_ubound]] - s_duringH[ind_free[out_of_ubound]]),nargout=2)
                                norms_b=abs(s_deltaH[out_of_ubound[ind_b_u]])
                                delta_b=abs(ub[ind_free[out_of_ubound[ind_b_u]]] - s[ind_free[out_of_ubound[ind_b_u]]])
                                ind_b=out_of_ubound[ind_b_u]
                                sign_b=sign_(ub[ind_free[out_of_ubound[ind_b_u]]] - s[ind_free[out_of_ubound[ind_b_u]]])
#                                out_of_ubound_init=matlabarray([[out_of_ubound_init],[out_of_ubound]])
                                out_of_ubound_init=concatenate_([out_of_ubound_init,out_of_ubound],axis=0)
                            if (length_(out_of_lbound) != 0):
                                diff_b_l,ind_b_l=max_(abs(s_duringH[ind_free[out_of_lbound]] - lb[ind_free[out_of_lbound]]),nargout=2)
                                norms_b=abs(s_deltaH[out_of_lbound[ind_b_l]])
                                delta_b=abs(lb[ind_free[out_of_lbound[ind_b_l]]] - s[ind_free[out_of_lbound[ind_b_l]]])
                                ind_b=out_of_lbound[ind_b_l]
                                sign_b=sign_(lb[ind_free[out_of_lbound[ind_b_l]]] - s[ind_free[out_of_lbound[ind_b_l]]])
#                                out_of_lbound_init=matlabarray([[out_of_lbound_init],[out_of_lbound]])
                                out_of_lbound_init=concatenate_([out_of_lbound_init,out_of_lbound],axis=0)
                            if ((length_(out_of_ubound) != 0) and (length_(out_of_lbound) != 0)):
                                if (diff_b_u > diff_b_l):
                                    norms_b=abs(s_deltaH[out_of_ubound[ind_b_u]])
                                    delta_b=abs(ub[ind_free[out_of_ubound[ind_b_u]]] - s[ind_free[out_of_ubound[ind_b_u]]])
                                    ind_b=out_of_ubound[ind_b_u]
                                    sign_b=sign_(ub[ind_free[out_of_ubound[ind_b_u]]] - s[ind_free[out_of_ubound[ind_b_u]]])
                                else:
                                    norms_b=abs(s_deltaH[out_of_lbound[ind_b_l]])
                                    delta_b=abs(lb[ind_free[out_of_lbound[ind_b_l]]] - s[ind_free[out_of_lbound[ind_b_l]]])
                                    ind_b=out_of_lbound[ind_b_l]
                                    sign_b=sign_(lb[ind_free[out_of_lbound[ind_b_l]]] - s[ind_free[out_of_lbound[ind_b_l]]])
                        else:
                            outside=0
                            if (length_(out_of_ubound_init) != 0):
                                diff_b_u,ind_b_u=min_(abs(ub[ind_free[out_of_ubound_init]] - s_duringH[ind_free[out_of_ubound_init]]),nargout=2)
                                norms_b=abs(s_deltaH[out_of_ubound_init[ind_b_u]])
                                delta_b=abs(ub[ind_free[out_of_ubound_init[ind_b_u]]] - s[ind_free[out_of_ubound_init[ind_b_u]]])
                                ind_b=out_of_ubound_init[ind_b_u]
                                sign_b=sign_(ub[ind_free[out_of_ubound_init[ind_b_u]]] - s[ind_free[out_of_ubound_init[ind_b_u]]])
                            if (length_(out_of_lbound_init) != 0):
                                diff_b_l,ind_b_l=min_(abs(s_duringH[ind_free[out_of_lbound_init]] - lb[ind_free[out_of_lbound_init]]),nargout=2)
                                norms_b=abs(s_deltaH[out_of_lbound_init[ind_b_l]])
                                delta_b=abs(lb[ind_free[out_of_lbound_init[ind_b_l]]] - s[ind_free[out_of_lbound_init[ind_b_l]]])
                                ind_b=out_of_lbound_init[ind_b_l]
                                sign_b=sign_(lb[ind_free[out_of_lbound_init[ind_b_l]]] - s[ind_free[out_of_lbound_init[ind_b_l]]])
                            if ((length_(out_of_ubound_init) != 0) and (length_(out_of_lbound_init) != 0)):
                                if (diff_b_u < diff_b_l):
                                    norms_b=abs(s_deltaH[out_of_ubound_init[ind_b_u]])
                                    delta_b=abs(ub[ind_free[out_of_ubound_init[ind_b_u]]] - s[ind_free[out_of_ubound_init[ind_b_u]]])
                                    ind_b=out_of_ubound_init[ind_b_u]
                                    sign_b=sign_(ub[ind_free[out_of_ubound_init[ind_b_u]]] - s[ind_free[out_of_ubound_init[ind_b_u]]])
                                else:
                                    norms_b=abs(s_deltaH[out_of_lbound_init[ind_b_l]])
                                    delta_b=abs(lb[ind_free[out_of_lbound_init[ind_b_l]]] - s[ind_free[out_of_lbound_init[ind_b_l]]])
                                    ind_b=out_of_lbound_init[ind_b_l]
                                    sign_b=sign_(lb[ind_free[out_of_lbound_init[ind_b_l]]] - s[ind_free[out_of_lbound_init[ind_b_l]]])
                    if (verbose):
                        lambda_save.append(_lambda)
                        norms_b_save.append(norms_b)
                        if (outside == 0):
                            fprintf_('%s%d%s %12.8e %s %12.8e %s\n' % (' bcdfo_solve_TR_MS_bc (',i,'): |s_i| = ',norms_b,'  |bound_i| = ',delta_b,'   s < bounds'))
                        else:
                            fprintf_('%s%d%s %12.8e %s %12.8e\n' % (' bcdfo_solve_TR_MS_bc (',i,'): |s_i| = ',norms_b,'  |bound_i| = ',delta_b))
                    out_of_uEpsbound=find_((ub[ind_free] - s_duringH[ind_free]) < - eps_bound)
                    out_of_lEpsbound=find_((s_duringH[ind_free] - lb[ind_free]) < - eps_bound)
                    if (isempty_(out_of_uEpsbound) and isempty_(out_of_lEpsbound)):
                        if (verbose >= 2):
                            disp_('all components inside the bounds + eps_bound')
                        back_inside=1
                        ind_u_active=find_(abs(ub[ind_free] - s_duringH[ind_free]) <= eps_bound)
                        ind_l_active=find_(abs(s_duringH[ind_free] - lb[ind_free]) <= eps_bound)
                        if ((length_(ind_u_active) + length_(ind_l_active)) != 0):
                            if (verbose >= 2):
                                disp_('all components inside the bounds + eps_bound, ',num2str_(length_(ind_u_active) + length_(ind_l_active)),' component/s close to one of its bounds')
                            s_afterH=s + I[:,ind_free] .dot( s_deltaH )
                            if (length_(ind_u_active) > 0):
                                s_afterH[ind_free[ind_u_active]]=ub[ind_free[ind_u_active]]
                            if (length_(ind_l_active) > 0):
                                s_afterH[ind_free[ind_l_active]]=lb[ind_free[ind_l_active]]
                            msg='boundary solution'
                            break
                    if (stratLam == 0):
                        if (back_inside == 0):
                            _lambda=2 * _lambda
                            if (upper < _lambda):
                                upper=2 * _lambda
                        else:
                            if (isempty_(out_of_ubound) and isempty_(out_of_lbound)):
                                upper=copy(_lambda)
                            else:
                                lower=copy(_lambda)
                            new_lambda=(_lambda + old_lambda) / 2
                            theta_range=theta * (upper - lower)
                            if (new_lambda > lower + theta_range and new_lambda < upper - theta_range):
                                _lambda=copy(new_lambda)
                            else:
                                _lambda=max_(sqrt_(lower * upper),lower + theta_range)
                    else:
                        if (isempty_(out_of_ubound) and isempty_(out_of_lbound)):
                            upper=copy(_lambda)
                        else:
                            lower=copy(_lambda)
                        unitvec=zeros_(nfree,1)
                        unitvec[ind_b]=1
                        es=unitvec.T .dot( s_deltaH )
                        if (sign_(es) != sign_b):
                            new_lambda=(lower + upper) / 2
                        else:
                            w1=numpy.linalg.solve(R.T,unitvec)
                            w2=numpy.linalg.solve(R.T,s_deltaH)
                            new_lambda=_lambda + ((norms_b - delta_b) / delta_b) * (norms_b ** 2 / (es .dot( (w1.T .dot( w2)))))
                            if (back_inside == 0 and upper <= new_lambda):
                                upper=2 * new_lambda
                        theta_range=theta * (upper - lower)
                        if (new_lambda > lower + theta_range and new_lambda <= upper - theta_range):
                            _lambda=copy(new_lambda)
                        else:
                            _lambda=real_(max_(sqrt_(lower * upper),lower + theta_range))
                else:
                    if (verbose):
                        disp_('unsuccessful factorization')
                    hardcase=0
                    lower=copy(_lambda)
                    t=0.5
                    _lambda=(1 - t) * lower + t * upper
                if (i >= 100):
                    s[0:n]=0.0
                    norms=0
                    msg='limit in bc-MS exceeded'
                    if (verbose):
                        disp_('Error in bcdfo_solve_TR_MS_bc: iteration limit in bc-MS exceeded!')
                        disp_('bcdfo_solve_TR_MS_bc: exit!')
                    return s,_lambda,norms,value,gplus,nfact,neigd,msg

        else:
            if (verbose >= 2):
                disp_('step inside bounds!')
            msg='(partly) interior solution'
            s=copy(s_after_reduced_ms)
            norms=norm_(s)
            value=0.5 * s.T.dot( H .dot( s )) + s.T .dot( g0 )
            if (verbose):
                disp_('bcdfo_solve_TR_MS_bc: exit!')
            return s,_lambda,norms,value,gplus,nfact,neigd,msg
        s=copy(s_afterH)
        norms=norm_(s)
        value=0.5 * s.T .dot( H .dot( s )) + s.T .dot( g0 )
        Delta=Delta0 - norms
        if (Delta < - eps_bound):
            disp_('Error in bcdfo_solve_TR_MS_bc: delta smaller than zero !!!!!!')
            msg='error7'
            if (verbose):
                disp_('bcdfo_solve_TR_MS_bc: exit!')
            return s,_lambda,norms,value,gplus,nfact,neigd,msg
        else:
            if (Delta < 0):
                if (verbose):
                    disp_('bcdfo_solve_TR_MS_bc: exit!')
                return s,_lambda,norms,value,gplus,nfact,neigd,msg
        g=g0 + H .dot( s )
#        ind_active=find_((ub - s) <= eps_bound or (s - lb) <= eps_bound)
        ind_active=find_(logical_or((ub - s) <= eps_bound, (s - lb) <= eps_bound))
        ind_active=ind_active.T
        ind_free=setdiff_(range(0,n),ind_active)
        nfree=length_(ind_free)
        if (nfree > 0):
            ng_reduced=norm_(g[ind_free],inf)
            if (ng_reduced <= 1e-05):
                if (verbose >= 2):
                    disp_('point first order critical - return')
                    ng_reduced
                if (verbose):
                    disp_('bcdfo_solve_TR_MS_bc: exit!')
                return s,_lambda,norms,value,gplus,nfact,neigd,msg
#            ind_g_crit=find_((abs(lb[ind_free]) <= 1e-10 and g[ind_free] > 0) or (ub[ind_free] <= 1e-10 and g[ind_free] < 0))
            ind_g_crit=find_(logical_or(logical_and(abs(lb[ind_free]) <= 1e-10, g[ind_free] > 0), logical_and(ub[ind_free] <= 1e-10, g[ind_free] < 0)))
            if (length_(ind_g_crit) != 0):
#                ind_active=matlabarray([ind_active,ind_free[ind_g_crit]])
                ind_active=concatenate_([ind_active,ind_free[ind_g_crit].T],axis=1)
                ind_free=setdiff_(ind_free,ind_active)
                nfree=length_(ind_free)

    return s,_lambda,norms,value,gplus,nfact,neigd,msg
