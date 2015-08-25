# -*- coding: utf-8 -*-
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
    
from bcdfo_find_new_yj import bcdfo_find_new_yj_
    
def bcdfo_poisedness_Y_(QZ=None,RZ=None,Y=None,eps_L=None,xbase=None,lSolver=None,whichmodel=None,hardcons=None,xl=None,xu=None,indfree=None,stratLam=None,scale=None,shift_Y=None,*args,**kwargs):
    
    """
   Computes the poisedness of the interpolation set Y in a ball of radius
   Delta centered at Y(:,1), assuming that a QR factorization of the
   interpolation system is known (and given by QZ, RZ).
   Poisedness is defined here as the maximum aboslute value of
   the Lagrange polynomials taken over the ball and for all polynomials.

   INPUT:

   QZ          : the Q matrix of the QR decomposition of Z(Y)
   RZ          : the R matrix of the QR decomposition of Z(Y)
   Y           : a matrix whose columns contain the current interpolation points
   eps_L       : the relative accuracy on the trust-region constraint for
                maximization of the Lagrange polynomials
   xbase       : the current base point
   lSolver     : linear solver used for the minimization of the model
   whichmodel  : kind of model/Lagrange polynomial to compute
   scale       : the current interpolation set scaling
   shift_Y     : 0 if no shift in interpolation points, 1 otherwise

   OUTPUT:

   lambda      : the poisedness of the interpolation set Y
   Y_radius    : the poisedness radius for Y, that is the largest distance
                from Y(:,i) to Y(:,1) (the base point)

   PROGRAMMING: Ph. Toint and S. Gratton, April 2009. (This version 22 VI 2009)

   DEPENDENCIES: bcdfo_find_new_yj, bcdfo_find_new_yj_bc

   CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
    """
#    varargin = cellarray(args)
#    nargin = 14-[QZ,RZ,Y,eps_L,xbase,lSolver,whichmodel,hardcons,xl,xu,indfree,stratLam,scale,shift_Y].count(None)+len(args)

    _lambda=0
    n,p1=size_(Y,nargout=2)
    Y_radius=0
    for j in range(1,p1): #it is indexed range(2,p1) in matlab
        Y_radius=max_(Y_radius,norm_(Y[:,j] - Y[:,0]))
    for j in range(1,p1): #it is indexed range(2,p1) in matlab
        if (hardcons == 1):
            y,improvement, msgTR=bcdfo_find_new_yj_bc_(QZ,RZ,Y,j,Y_radius,eps_L,xbase,lSolver,whichmodel,xl,xu,indfree,stratLam,scale,shift_Y,nargout=2)
        else:
            y,improvement, msgTR=bcdfo_find_new_yj_(QZ,RZ,Y,j,Y_radius,eps_L,xbase,lSolver,whichmodel,scale,shift_Y,nargout=2)
        _lambda=max_(improvement,_lambda)
    return _lambda,Y_radius