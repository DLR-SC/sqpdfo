# -*- coding: utf-8 -*-
from __future__ import division
#try:
from runtime import *
from ecdfo_global_variables import *
from numpy import array
#except ImportError:
#from smop.runtime import *

class tcgInfo():
	def __init__(self):
		self.flag = None# 0
		self.iter = None# 2
		self.prec = None# 0
		self.curv = None# 1
  
  
def sqplab_tcg_(A=None,b=None,delta=None,max_iter=None,tol=None,plevel=None,fout=None,*args,**kwargs):
    """
    % [x,info] = sqplab_tcg (A,b,delta,max_iter,tol,plevel,fout);
%
% Solves A x = b for x, by Steihaug's conjugate gradient (CG) method,
% where A is a symmetric (possibly indefinite) matrix. CG iterations
% starts from x = 0, up to one of the situations described below by info
% occurs.
%
% On entry
%   A: symmetric matrix of the linear system
%   b: RHS of the linear system
%   delta: trust region radius
%   max_iter: maximum iterations allowed
%   tol: x is considered a solution if the Euclidean norm of the
%     gradient A*x-b is less than tol.
%   plevel: printing level
%     =0: don't print anything
%     >0; print on fout
%   fout: printing channel
%
% On return
%   info is a structure providing information on the run
%     info.curv is the last encountered curvature d'*Ad/(d'*d) (not a
%       field if no iteration has been performed)
%     info.flag is the return code
%       = 0: convergence is obtained, up to the given tolerance tol: the
%            final x is the solution of the LS,
%       =-1: stop on max_iter,
%       = 1: the boundary of the ball of radius delta is encountered:
%            the final x is on the boundary of the ball
%       = 2: a negative curvature direction has been encountered: the
%            final x is on the boundary of the ball.
%     info.iter is the number of CG iterations
%     info.prec is the final precision, the l2-norm of the residual A*x-b
%   x is the computed approximate solution

%-----------------------------------------------------------------------
%
% Author: Jean Charles Gilbert, INRIA.
%
% Copyright 2008, 2009, INRIA.
%
% SQPlab is distributed under the terms of the Q Public License version
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
%-----------------------------------------------------------------------
    """
#    varargin = cellarray(args)
#    nargin = 7-[A,b,delta,max_iter,tol,plevel,fout].count(None)+len(args)

    info=tcgInfo()

    x=zeros_(size_(b))
    cost=0
    g=- b
    g2=g.T.dot(g)
    tol2=tol*tol
    delta2=delta * delta
    dAd=array([])
    if plevel:
        fprintf_(fout,'    TCG solver; required tolerance %8.2e\n'%(tol))
        fprintf_(fout,'    iter       cost        |res|   curvature  stepsize   |step|\n')
    set_iter(0)
    while 1:

        set_iter(get_iter() + 1)
        if plevel:
            fprintf_(fout,'    %4i  %14.7e  %7.1e'%(get_iter(),cost,sqrt_(g2)))
        if g2 <= tol2:
            info.flag=0
            break
        if get_iter() > max_iter:
            set_iter(copy(max_iter))
            info.flag=- 1
            break
        if get_iter() == 1:
            d=- g
        else:
            d=- g + (g2 / g2_) * d
        Ad=A.dot(d)
        dAd=d.T.dot(Ad)
        if plevel:
            fprintf_(fout,'  %9.2e'%(dAd / (d.T.dot(d))))
        if dAd <= 0:
            x,alpha=dogleg_(x,x + d,delta,nargout=2)
            info.flag=2
            if plevel:
                fprintf_(fout,'  %8.2e  %8.2e\n'%(alpha,norm_(x)))
                cost=0.5 * (x.T.dot(A.dot(x))) - b.T.dot(x)
                fprintf_(fout,'    %4i  %14.7e\n'%(get_iter() + 1,cost))
            break
        alpha=- (g.T.dot(d)) / dAd
        xx=x + alpha * d
        if plevel:
            fprintf_(fout,'  %8.2e'%(alpha))
        if xx.T.dot(xx) > delta2:
            x,alpha=dogleg_(x,xx,delta,nargout=2)
            info.flag=1
            if plevel:
                fprintf_(fout,'  %8.2e  %8.2e\n'%(alpha,norm_(x)))
                cost=0.5 * (x.T.dot(A.dot(x))) - b.T.dot(x)
                fprintf_(fout,'    %4i  %14.7ee\n'%(get_iter()+ 1,cost))
            break
        else:
            x=copy(xx)
        if plevel:
            fprintf_(fout,'  %8.2e\n'%(norm_(x)))
        g=g + alpha * Ad
        g2_=copy(g2)
        g2=g.T.dot(g)
        if plevel:
            cost=0.5 * (x.T.dot((g - b)))

    info.iter=get_iter()
    info.prec=sqrt_(g2)
    if not isempty_(dAd):
        info.curv=dAd / (d.T.dot(d))
    return x,info

def dogleg_(dc=None,dn=None,delta=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 3-[dc,dn,delta].count(None)+len(args)

    dd=dn - dc
    aa=dd.T.dot(dd)
    if aa == 0:
        dd=copy(dc)
        t=0
        return dd,t
    bb=dc.T.dot(dd)
    cc=dc.T.dot(dc) - delta ** 2
    if cc >= 0:
        dd=copy(dc)
        t=0
        return dd,t
    t=(sqrt_(bb ** 2 - aa * cc) - bb) / aa
    dd=dc + t*(dd)
    return dd,t
