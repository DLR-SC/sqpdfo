# -*- coding: utf-8 -*-
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from bcdfo_evalZ import *
from numpy import array
from copy import copy

def bcdfo_computeP_(QZ=None,RZ=None,Y_=None,fY=None,whichmodel=None,P_old=None,ind_Y=None,i_xold=None,i_xplus=None,g=None,scale=None,shift_Y=None,Delta0=None,*args,**kwargs):
    """
#%
#%  Computes the polynomial P, where P is then represented by a row vector
#%  containing its coefficients for the successive monomials.
#%  More specifically, these values are:
#%  P(1)        : constant coefficient,
#%  P(2:n+1)    : coefficients for the linear terms in x(1)... x(n),
#%  P(n+2:2n+2) : coefficients for the squared terms in x(1)^2 ... x(n)^2
#%  P(2n+3,3n+2): coefficients for the quadratic terms of the first subdiagonal:
#%                in x(1)*x(2) ... x(n-1)*x(n)
#%  P(3n+3,4n+1): coefficients for the quadratic terms of the second subdiagonal:
#%                in x(1)*x(3) ... x(n-2)*x(n)
#%  etc.
#%
#%  INPUTS:
#%
#%  QZ, RZ      : the QR factors of the (possibly shifted) matrix Z containing
#%                the polynomial expansion of the interpolation points
#%  Y           : current interpolation set
#%  fY          : function values of the interpolation points
#%  whichmodel  : the kind of model to build
#%  P_old       : polynomial coefficients of the former polynomial
#%  scale       : scaling of the interpolation matrix
#%  shift_Y     : shift of the matrix of interpolation points
#%
#%  OUTPUT:
#%
#%  P           : a row vector containing the coefficients of the polynomial
#%
#%  DEPENDENCIES: -
#%
#%  PROGRAMMING: A. Troeltzsch, September 2010.
#%
#%  TEST:
#%  Y = [ 0 1 0 2 0 ; 0 0 1 0 2 ]; fY = [ 1 2 3 4 5 ];
#%  whichmodel = 0;
#%  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y(  Y, whichmodel, 1 );
#%  P = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, [0, 0, 0, 0, 0], 0, ...
#%      1, 1, [0, 0], scale, 1, 1 )
#%  should give:
#%  P =
#%   1.0000   1.0000   4.0000   4.0000       0
#%
#%  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#%
    """
#    varargin = cellarray(args)
#    nargin = 13-[QZ,RZ,Y,fY,whichmodel,P_old,ind_Y,i_xold,i_xplus,g,scale,shift_Y,Delta0].count(None)+len(args)

    #Copy only if necessary
    if (shift_Y):
        Y=copy(Y_)
    else:
        Y=Y_

    n,p1=size_(Y,nargout=2)
    badcond=0
    q=((n + 1) * (n + 2)) / 2
    if (whichmodel == 3 and p1 < q):
        whichmodel=2
    if (whichmodel == 0):
        P=(QZ.dot(numpy.linalg.solve(RZ.T,fY.T))).T
    else:
        if (whichmodel == 1):
            n_rhs=size_(fY,1)
            if (p1 == n + 1 or p1 == q):
                P[0:n_rhs,0:p1]=(numpy.linalg.solve(RZ,(QZ.T.dot(fY.T)))).T
                if (p1 == n + 1):
                    P[0:n_rhs,n + 1:q]=0.0
            else:
                if (shift_Y):
                    xbase=copy(Y[:,0])
                    scaleY=0
                    for i in range(0,p1):
                        Y[:,i]=Y[:,i] - xbase
                        scaleY=max_(scaleY,norm_(Y[:,i]))
                    Y=Y / scaleY
                P=array([])
                for i in range(0,n_rhs):
                    rhs=concatenate_([fY[i,:],zeros_(1,n + 1)], axis=1)
#                    rhs=array([fY[i,:],zeros_(1,n + 1)])
                    mualpha=(numpy.linalg.solve(RZ,(QZ.T.dot(rhs.T)))).T
                    P_i[0:n + 1]=mualpha[p1:p1 + n + 1].T
                    M=bcdfo_evalZ_(Y,q).T
                    P_i[n + 1:q]=M[:,n + 1:q].T.dot(mualpha[0:p1].T)
                    P=concatenate_([P, P_i])
#                    P=array([[P],[P_i]])
        else:
            if (whichmodel == 2):
                P=(QZ.dot((numpy.linalg.solve(RZ.T,fY.T)))).T
            else:
                if (whichmodel == 3):
                    P=(pinv_(RZ).dot(QZ.T.dot(fY.T))).T
    return P