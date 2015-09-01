# -*- coding: utf-8 -*-
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from bcdfo_computeLj import bcdfo_computeLj_
from bcdfo_gradP import bcdfo_gradP_
from bcdfo_hessP import bcdfo_hessP_
from bcdfo_solve_TR_MS import bcdfo_solve_TR_MS_
from numpy import *



def bcdfo_find_new_yj_(QZ=None,RZ=None,Y=None,j=None,Delta=None,eps_L=None,xbase=None,lSolver=None,whichmodel=None,scale=None,shift_Y=None,*args,**kwargs):
    """
#
#  Computes a point which is best to replace yj, the j-th (j>1) column of Y (the
#  base for the j-th polynomial) in a ball of radius Delta centered at the
#  first column of Y.  This is achieved by maximizing the absolute value of
#  the j-th Lagrange polynomial in that ball.
#
#  For conventions on how polynomals are represented, see the documentation of
#  evalZ.
#
#  INPUT:
#
#  QZ          : the Q matrix of the QR decomposition of Z(Y)
#  RZ          : the R matrix of the QR decomposition of Z(Y)
#  Y           : a matrix whose columns contain the current interpolation points
#  j           : the index of the interpolation point one wishes to replace (j > 1)
#  Delta       : the radius of the ball centered at Y(:,1) in which the
#                replacement point must be found
#  eps_L       : the relative accuracy on the trust-region constraint for
#                maximization of the Lagrange polynomials
#  xbase       : the current base point
#  lSolver     : linear solver used for the minimization of the model
#  whichmodel  : kind of model/Lagrange polynomial to compute
#  scale       : the current interpolation set scaling
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#
#  OUTPUT:
#
#  ynew        : the best replacement for Y(:,j)
#  improvement : the improvement in poisedness obtained by the update, which
#                is equal to |L_j(new y)|. If this value is smaller than the
#                threshold input parameter, L and X are unchanged by the
#                procedure.
#
#  PROGRAMMING: Ph. Toint, February 2009. (This version 22 VI 2009)
#
#  USES: bcdfo_gradP, bcdfo_hessP, bcdfo_solve_TR_MS
#
#  TEST:
#  Y = [ 3 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; whichmodel = 0;
#  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y , whichmodel, 0 );
#  [ ynew, improvement ] = bcdfo_find_new_yj( QZ, RZ, Y, 5, 1.0, 0.001, xbase, 1, ...
#       whichmodel, scale, 0 )
#  should give
#  ynew =
#
#    3.2808
#   -0.9598
#
#  improvement =
#
#  314.8825
#
#  the same must be obtained by the shifted and scaled version:
#  Y = [ 3 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; whichmodel = 0;
#  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y , whichmodel, 1 );
#  [ ynew, improvement ] = bcdfo_find_new_yj( QZ, RZ, Y, 5, 1.0, 0.001, xbase, 1, ...
#       whichmodel, scale, 1 )
#
    """
#    varargin = cellarray(args)
#    nargin = 11-[QZ,RZ,Y,j,Delta,eps_L,xbase,lSolver,whichmodel,scale,shift_Y].count(None)+len(args)

    verbose=0 # 1 for debug
    n=size_(Y,1)
    ynew=zeros_(1,n)
    improvement=0
    msgTR=''
    if (verbose):
        disp_('--------- enter find_new_yj ')
    if (j < 1): # never attempt to replace the current iterate.
        return ynew,improvement,msgTR
    #  Get the j-th Lagrange polynomial 
    Lj=bcdfo_computeLj_(QZ,RZ,j,Y,whichmodel,scale,shift_Y)
    if (length_(find_(isnan(Lj))) != 0 or length_(find_(~ isreal(Lj))) != 0 or length_(find_(isinf(Lj))) != 0):
        msgTR='Error0: Lagrange polynomial contains NaN or Inf or nonreal components!!'
        if (verbose):
            disp_(msgTR)
        return ynew,improvement,msgTR
#  Maximize Lj in a larger 2-norm TR if using infty-norm in the local solver (CG)
    if (lSolver == 2):
        Delta=sqrt_(n) * Delta
#  Get the polynomial's gradient and Hessian at the current iterate.
    if (shift_Y):


#     When shifted, the origin in the scaled variables corresponds 
#     to Y(:,0) in the unscaled space
        g=bcdfo_gradP_(Lj,zeros_(n,1),xbase,scale,0)
        H=bcdfo_hessP_(Lj,zeros_(n,1),xbase,scale,0)
#     Minimize this polynomial and its opposite.
        pstep,_lambda,norms,pvalue,gplus,nfact,neigd,msgTR,hardcase=bcdfo_solve_TR_MS_(g,H,Delta * scale[1],eps_L,nargout=9)
        pstep=pstep / scale[1]
        mstep,_lambda,norms,mvalue,gplus,nfact,neigd,msgTR,hardcase=bcdfo_solve_TR_MS_(- g,- H,Delta * scale[1],eps_L,nargout=9)
        mstep=mstep / scale[1]
    else:
#     When no shift occurs, the current iterate is Y(:,1)
        g=bcdfo_gradP_(Lj,Y[:,[0]],xbase,scale,0)
        H=bcdfo_hessP_(Lj,Y[:,[0]],xbase,scale,0)
#     Minimize this polynomial and its opposite.
        pstep,_lambda,norms,pvalue,gplus,nfact,neigd,msgTR,hardcase=bcdfo_solve_TR_MS_(g,H,Delta,eps_L,nargout=9)
        mstep,_lambda,norms,mvalue,gplus,nfact,neigd,msgTR,hardcase=bcdfo_solve_TR_MS_(- g,- H,Delta,eps_L,nargout=9)
    if (verbose):
        disp_(' === find_new_yj: j = ',str(j),' positive value = ',str(pvalue),' step:')
        pstep.T
        disp_(' === find_new_yj: j = ',str(j),' negative value = ',str(mvalue),' step:')
        mstep.T
#  Select the maximum in absolute value.
    if (mvalue < pvalue):
        improvement=abs(mvalue)
        ynew=Y[:,0].reshape(-1,1) + mstep
    else:
        improvement=abs(pvalue)
        ynew=Y[:,0].reshape(-1,1) + pstep
    if (verbose):
        disp_('--------- exit find_new_yj ')
    return ynew,improvement,msgTR