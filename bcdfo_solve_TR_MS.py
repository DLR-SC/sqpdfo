# -*- coding: utf-8 -*-
#from __future__ import division
from runtime import *
from numpy import *
from copy import copy
#from numpy import inf
def bcdfo_solve_TR_MS_(g=None,H=None,Delta=None,eps_D=None,*args,**kwargs):
    """
#
#%  A simple implementation of exact trust-region minimization based on the
#%  MorÃ©-Sorensen algorithm.
#
#%  INPUT: 
#
#%  g        : the model's gradient
#%  H        : the model's Hessian
#%  Delta    : the trust-region's radius
#%  eps_D    : the accuracy required on the equation ||s|| = Delta for a
#%             boundary solution
#
#%  OUTPUT:
#
#%  s        : the trust-region step
#%  lambda   : the Lagrange multiplier corresponding to the trust-region constraint
#%  norms    : the norm of s
#%  value    : the value of the model at the optimal solution
#%  gplus    : the value of the model's gradient at the optimal solution
#%  nfact    : the number of Cholesky factorization required
#%  neigd    : the number of eigenvalue decompositions required
#%  msg      : an information message
#
#%  DEPENDENCIES: -
#
#%  PROGRAMMING: Ph. Toint and S. Gratton, April 2009. (This version 14 I 2010)
#
#%  TEST:
#
#%  bcdfo_solve_TR_MS( [ 2 ; 3 ], [ 4 6; 6 5 ], 1.0, 0.001 )
#%  should give
#%    0.5153
#%   -0.8575
#%  bcdfo_solve_TR_MS( [ 2 ; 0 ], [ 4 0; 0 -15 ], 1.0, 0.001 )
#%  should give
#%   -0.1053
#%    0.9944
#
#%  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
    """

#    varargin = cellarray(args)
#    nargin = 4-[g,H,Delta,eps_D].count(None)+len(args)

    verbose=0
    theta=1e-13
    epsilon=1e-12
    nitmax=300
    n=length_(g)
    s=zeros_(n,1)
    norms=0
    _lambda=0
    value=0
    gplus=copy(g)
    nfact=0
    neigd=0
    hardcase=0
    if (verbose):
        disp_(' bcdfo_solve_TR_MS : ============ enter')
    if (length_(find_(isnan(H))) != 0):
        disp_(' bcdfo_solve_TR_MS : H contains NaNs!')
        msg='error1'
        if (verbose):
            disp_(' bcdfo_solve_TR_MS : ============ error exit')
        return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
    if (length_(find_( ~isreal(H))) != 0):
        disp_(' bcdfo_solve_TR_MS : H contains imaginary parts!')
        msg='error2'
        if (verbose):
            disp_(' bcdfo_solve_TR_MS : ============ error exit')
        return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
    if (length_(find_(isinf(H))) != 0):
        disp_(' bcdfo_solve_TR_MS : H contains infinite elements!')
        msg='error3'
        if (verbose):
            disp_(' bcdfo_solve_TR_MS : ============ error exit')
        return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
    gnorm=norm_(g)
    goverD=gnorm / Delta
    Hnorminf=norm_(H,inf)
    if (Hnorminf > 0):
        HnormF=norm_(H,'fro')
    else:
        HnormF=0
    lower=max_(0,goverD - min_(Hnorminf,HnormF))
    upper=max_(0,goverD + min_(Hnorminf,HnormF))
    Dlower=(1 - eps_D) * Delta
    Dupper=(1 + eps_D) * Delta
    if Delta == 0:
        msg='bcdfo_solve_TR_MS : trust region is zero - exit !'
        if verbose:
            disp_(msg)
        sfound=1
        norms=norm_(s)
        _lambda=0
        return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
    if (gnorm ** 2 < epsilon):
        msg='zero gradient'
        if (verbose):
            disp_(' bcdfo_solve_TR_MS : ============ zero gradient:')
        sfound=1
        norms=norm_(s)
        _lambda=0
        return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
    else:
        if (verbose):
            disp_(' bcdfo_solve_TR_MS : ============ nonzero gradient:')
        if (lower == 0):
            _lambda=0
        else:
            _lambda=max_(sqrt_(lower * upper),lower + theta * (upper - lower))
        for i in range(0,nitmax):
            new_lambda=- 1
            sfound=0
            if (verbose):
                disp_(' bcdfo_solve_TR_MS (',int2str_(i),'): lower = ',num2str_(lower),' lambda = ',num2str_(_lambda),' upper = ',num2str_(upper))
            R,p=chol_(H + _lambda * eye_(n),nargout=2)
            if (length_(find_(isnan(R))) != 0):
                disp_(' bcdfo_solve_TR_MS : NaNs in Cholesky factorization')
                msg='error4'
                if (verbose):
                    disp_(' bcdfo_solve_TR_MS : ============ error exit')
                return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
            nfact=nfact + 1
            if (p == 0):
                s=numpy.linalg.solve(- R,(numpy.linalg.solve(R.T,g)))
                sfound=1
                norms=norm_(s)
                if (verbose):
                    disp_(' bcdfo_solve_TR_MS (',int2str_(i),'): ||s|| = ',num2str_(norms),' Delta  = ',num2str_(Delta))
                if ((_lambda <= epsilon and norms <= Dupper) or (norms >= Dlower and norms <= Dupper)):
                    w=H.dot(s)
                    value=g.T.dot(s) + 0.5 *s.T.dot(w)
                    gplus=g + w
                    norms=norm_(s)
                    if (norms < (1 - eps_D) * Delta):
                        msg='interior solution'
                    else:
                        msg='boundary solution'
                    if (verbose):
                        disp_(' bcdfo_solve_TR_MS : ============ successful exit')
                    return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
                w=numpy.linalg.solve(R.T,s)
                normw2=w.T.dot(w)
                new_lambda=_lambda + ((norms - Delta) / Delta) * (norms ** 2 / normw2)
                if (norms > Dupper):
                    lower=copy(_lambda)
                else:
                    upper=copy(_lambda)
                theta_range=theta * (upper - lower)
                if (new_lambda > lower + theta_range and new_lambda < upper - theta_range):
                    _lambda=copy(new_lambda)
                else:
                    _lambda=max_(sqrt_(lower * upper),lower + theta_range)
            else:
                lower=copy(_lambda)
                t=0.5
                _lambda=(1 - t) * lower + t * upper
            if (upper - lower < theta * max_(1,upper)):
                break
    V,D=eig_(H,nargout=2)
    neigd=neigd + 1
    mu,imu=min_(diag(D),nargout=2)
    if (verbose):
        gamma=abs_(V[:,imu].reshape(1,-1).dot( g))
        disp_(' bcdfo_solve_TR_MS : ============ pseudo hard case: gamma = ',num2str_(gamma),' ||g|| = ',num2str_(norm_(g)))
    D=D - mu * eye_(n)
    maxdiag=max_(diag(D))
    ii=find_(abs(diag(D)) < 1e-10 * maxdiag)
    if (length_(ii) < n and  not(isempty_(ii))):
        D[ii,ii.T]=0.5 * maxdiag * eye_(length_(ii))
        Dinv=inv_(D)
        Dinv[ii,ii.T]=0
        scri=- V .dot( Dinv.dot( V.T .dot(g)))
        nscri=norm_(scri)
    else:
        scri=zeros_(n,1)
        nscri=0
    if (nscri <= Delta):
        p2=poly1d([norm_(V[:,imu]) ** 2,2 * V[:,imu].reshape(1,-1).dot( scri),nscri ** 2 - Delta ** 2])
        root=max(p2.r)
        s=scri + root * V[:,imu].reshape(-1,1)
    else:
        s=Delta * scri / nscri
    _lambda=- mu
    if (verbose):
        disp_(' bcdfo_solve_TR_MS : ============ ||scri|| = ',num2str_(norm_(scri)),' lambda = ',num2str_(_lambda))
    hardcase=1
    w=H.dot(s)
    value=g.T .dot( s ) + 0.5 *s.T.dot(w)
    gplus=g + w
    norms=norm_(s)
    if abs(value) <= 1e-15:
        s=zeros_(size_(s))
    if (norms < (1 - eps_D) * Delta):
        msg = 'interior solution ( '+ str( nfact )+' factorizations,  lambda = '+ str( _lambda )+ ')' 
    else:
        msg = 'boundary solution ( '+ str( nfact )+ ' factorizations, '+str( neigd )+ ' eigen decomposition, lambda = '+ str( _lambda )+ ' )'
        if (verbose):
            disp_(' bcdfo_solve_TR_MS : ============ hard case exit')
    return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
