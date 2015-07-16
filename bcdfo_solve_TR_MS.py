# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 12:41:45 2015

@author: lien_ol
"""

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#from __future__ import division
from runtime import *
#from numpy import inf
def bcdfo_solve_TR_MS_(g=None,H=None,Delta=None,eps_D=None,*args,**kwargs):
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
    gplus=copy_(g)
    nfact=0
    neigd=0
    hardcase=0
    if (verbose):
        disp_(char(' bcdfo_solve_TR_MS : ============ enter'))
    if (length_(find_(isnan_(H))) != 0):
        disp_(char(' bcdfo_solve_TR_MS : H contains NaNs!'))
        msg=char('error1')
        if (verbose):
            disp_(char(' bcdfo_solve_TR_MS : ============ error exit'))
        return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
    if (length_(find_(~ isreal_(H))) != 0):
        disp_(char(' bcdfo_solve_TR_MS : H contains imaginary parts!'))
        msg=char('error2')
        if (verbose):
            disp_(char(' bcdfo_solve_TR_MS : ============ error exit'))
        return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
    if (length_(find_(isinf_(H))) != 0):
        disp_(char(' bcdfo_solve_TR_MS : H contains infinite elements!'))
        msg=char('error3')
        if (verbose):
            disp_(char(' bcdfo_solve_TR_MS : ============ error exit'))
        return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
    gnorm=norm_(g)
    goverD=gnorm / Delta
    Hnorminf=norm_(H,inf)
    if (Hnorminf > 0):
        HnormF=norm_(H,char('fro'))
    else:
        HnormF=0
    lower=max_(0,goverD - min_(Hnorminf,HnormF))
    upper=max_(0,goverD + min_(Hnorminf,HnormF))
    Dlower=(1 - eps_D) * Delta
    Dupper=(1 + eps_D) * Delta
    if Delta == 0:
        msg=char('bcdfo_solve_TR_MS : trust region is zero - exit !')
        if verbose:
            disp_(msg)
        sfound=1
        norms=norm_(s)
        _lambda=0
        return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
    if (gnorm ** 2 < epsilon):
        msg=char('zero gradient')
        if (verbose):
            disp_(char(' bcdfo_solve_TR_MS : ============ zero gradient:'))
        sfound=1
        norms=norm_(s)
        _lambda=0
        return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
    else:
        if (verbose):
            disp_(char(' bcdfo_solve_TR_MS : ============ nonzero gradient:'))
        if (lower == 0):
            _lambda=0
        else:
            _lambda=max_(sqrt_(lower * upper),lower + theta * (upper - lower))
        for i in arange_(1,nitmax).reshape(-1):
            new_lambda=- 1
            sfound=0
            if (verbose):
                disp_([char(' bcdfo_solve_TR_MS ('),int2str_(i),char('): lower = '),num2str_(lower),char(' lambda = '),num2str_(_lambda),char(' upper = '),num2str_(upper)])
            try:
                R = chol_( H + _lambda * eye_( n ) ).T
                p = 0
            except:
                R = matlabarray([[]])
                p = 1
            if (length_(find_(isnan_(R))) != 0):
                H
                _lambda
                norm_(g)
                R
                p
                disp_(char(' bcdfo_solve_TR_MS : NaNs in Cholesky factorization'))
                msg=char('error4')
                if (verbose):
                    disp_(char(' bcdfo_solve_TR_MS : ============ error exit'))
                return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
            nfact=nfact + 1
            if (p == 0):
                s=solve_(- R,(solve_(R.T,g)))
                sfound=1
                norms=norm_(s)
                if (verbose):
                    disp_([char(' bcdfo_solve_TR_MS ('),int2str_(i),char('): ||s|| = '),num2str_(norms),char(' Delta  = '),num2str_(Delta)])
                if ((_lambda <= epsilon and norms <= Dupper) or (norms >= Dlower and norms <= Dupper)):
                    w=H * s
                    value=g.T * s + 0.5 * s.T * w
                    gplus=g + w
                    norms=norm_(s)
                    if (norms < (1 - eps_D) * Delta):
                        msg=char('interior solution')
                    else:
                        msg=char('boundary solution')
                    if (verbose):
                        disp_(char(' bcdfo_solve_TR_MS : ============ successful exit'))
                    return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase
                w=solve_(R.T,s)
                normw2=w.T * w
                new_lambda=_lambda + ((norms - Delta) / Delta) * (norms ** 2 / normw2)
                if (norms > Dupper):
                    lower=copy_(_lambda)
                else:
                    upper=copy_(_lambda)
                theta_range=theta * (upper - lower)
                if (new_lambda > lower + theta_range and new_lambda < upper - theta_range):
                    _lambda=copy_(new_lambda)
                else:
                    _lambda=max_(sqrt_(lower * upper),lower + theta_range)
            else:
                lower=copy_(_lambda)
                t=0.5
                _lambda=(1 - t) * lower + t * upper
            if (upper - lower < theta * max_(1,upper)):
                break
    D,V=eig_(H,nargout=2)
    D=diag_(D)
    neigd=neigd + 1
    mu,imu=min_(diag_(D),nargout=2)
    if (verbose):
        gamma=abs_(V[:,imu].T * g)
        disp_([char(' bcdfo_solve_TR_MS : ============ pseudo hard case: gamma = '),num2str_(gamma),char(' ||g|| = '),num2str_(norm_(g))])
    D=D - mu * eye_(n)
    maxdiag=max_(diag_(D))
    ii=find_(abs_(diag_(D)) < 1e-10 * maxdiag)
    if (length_(ii) < n and isempty_(ii) != 1):
        D[ii,ii]=0.5 * maxdiag * eye_(length_(ii))
        Dinv=inv_(D)
        Dinv[ii,ii]=0
        scri=- V * Dinv * V.T * g
        nscri=norm_(scri)
    else:
        scri=zeros_(n,1)
        nscri=0
    if (nscri <= Delta):
        p2=poly1d_([norm_(V[:,imu]) ** 2,2 * V[:,imu]*scri,nscri ** 2 - Delta ** 2],r=0)
        root = max_(p2.r)
        s=scri + root* V[:,imu].T
    else:
        s=Delta * scri / nscri
    _lambda=- mu
    if (verbose):
        disp_([char(' bcdfo_solve_TR_MS : ============ ||scri|| = '),num2str_(norm_(scri)),char(' lambda = '),num2str_(_lambda)])
    hardcase=1
    w=H * s
    value=g.T * s + 0.5 * s.T * w
    gplus=g + w
    norms=norm_(s)
    if abs_(value) <= 1e-15:
        s=zeros_(size_(s))
    if ( norms < ( 1 - eps_D ) * Delta ):
      msg = [ 'interior solution ( '+ str( nfact )+' factorizations,  lambda = '+ str( _lambda )+ ')' ]
    else:
      msg = [ 'boundary solution ( '+ str( nfact )+ ' factorizations, '+str( neigd )+ ' eigen decomposition, lambda = '+ str( _lambda )+ ' )']
    if (verbose):
        disp_(char(' bcdfo_solve_TR_MS : ============ hard case exit'))
    return s,_lambda,norms,value,gplus,nfact,neigd,msg,hardcase