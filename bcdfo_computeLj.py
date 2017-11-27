# -*- coding: utf-8 -*-

from runtime import *

from bcdfo_evalZ import *
from copy import copy
from numpy import array
import numpy
def bcdfo_computeLj_(QZ=None,RZ=None,j=None,Y_=None,whichmodel=None,scale=None,shift_Y=None,*args,**kwargs):
    """
#
#  Computes the coefficients of the j-th Lagrange polynomial.
#
#  INPUTS:
#
#  QZ, RZ      : the QR factors of the (possibly shifted) matrix Z containing
#                the polynomial expansion of the interpolation points
#  j           : the index of the interpolation point one wishes to compute
#                the associated Lagrange polynomial
#  Y           : current interpolation set Y
#  whichmodel  : the method to compute the polynomial
#  scale       : scaling factor of the interpolation matrix
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#
#  OUTPUT:
#
#  Lj          : a row vector containing the coefficients of the polynomial
#
#  DEPENDENCIES: -
#
#  PROGRAMMING: A. Troeltzsch, September 2010.
#
#  TEST:
#  Y = [ 0 1 0 2 0 ; 0 0 1 0 2 ];  whichmodel = 0;
#  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y(  Y, whichmodel, 1 );
#  Lj = bcdfo_computeLj( QZ, RZ, 1, Y, whichmodel, scale, 1 )
#  should give:
#  Lj =
#   1.0000  -3.0000  -3.0000   4.0000   4.0000
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
    """
#    varargin = cellarray(args)
#    nargin = 7-[QZ,RZ,j,Y,whichmodel,scale,shift_Y].count(None)+len(args)

    #Copy only if necessary
    if (shift_Y):
        Y=copy(Y_)
    else:
        Y=Y_

    j=j+1
    n,p1=size_(Y,nargout=2)
    q=((n + 1) * (n + 2)) / 2
    if (whichmodel == 3 and p1 < q):
#          for underdetermined regression model use min l2-norm model
        whichmodel=2
    if (whichmodel == 0):
#    build (sub-basis) Lagrange polynomials (p1 = q)
#    (QZ and RZ are the factors of Z = M')
        Lj=(QZ.dot((numpy.linalg.solve(RZ.T,concatenate_([zeros_(j - 1,1),array([[1]]),zeros_(p1 - j,1)],axis=0))))).T
    elif (whichmodel == 1):
#    build mixed Lagrange polynomials:
#    Minimum Frobenius norm Lagrange polynomials (p1 <= q) or
#    L2-norm Lagrange polynomials (p1 == n+1 or p1 == q )
        if (p1 == n + 1 or p1 == q):
#       build linear/quadratic Lagrange polynomials 
           Lj=(numpy.linalg.solve(RZ,QZ.T).dot(concatenate_([zeros_(j - 1,1),array([[1]]),zeros_(p1 - j,1)],axis=0))).T
           if (p1 == n + 1):
                Lj[n + 1:q]=0
        else:
#       build Minimum Frobenius-norm Lagrange polynomials (p1 <= q)
#
#        If shifting is active, the matrix of the interpoaltion points is first
#        shifted and scaled (to compute M correctly).
            if (shift_Y):
                xbase=copy(Y[:,0])
                scaleY=0
                for i in range(0,p1):
                    Y[:,i]=Y[:,i] - xbase
                    scaleY=max_(scaleY,norm_(Y[:,i]))
                Y=Y / scaleY
#       compute the right-hand side of the system
            rhs=concatenate_([zeros_(j - 1,1),array([[1]]),zeros_(p1 + n +1 - j,1)],axis=0)
            mualpha=(numpy.linalg.solve(RZ,(QZ.T.dot(rhs)))).T
#       constant and linear part of P
            Lj[0:n + 1]=mualpha[p1:p1 + n + 1].T
#       quadratic part of P
            M=bcdfo_evalZ_(Y,q).T
            Lj[n + 1:q]=M[:,n + 1:q].T.dot( mualpha[0:p1]).T
    elif (whichmodel == 2):
#    build Minimum L2 norm Lagrange polynomials (p1 <= q)
#    (QZ and RZ are the factors of Z = M')
#    Take pseudo-inverse for underdetermined system because the result
#    is different from using backslash-operator
        if (p1 < q):
            Lj=(QZ.dot((pinv_(RZ.T).dot(concatenate_([zeros_(j - 1,1),array([[1]]),zeros_(p1 - j,1)],axis=0))))).T
        else:
            Lj=(QZ.dot((numpy.linalg.solve(RZ.T,concatenate_([zeros_(j - 1,1),array([[1]]),zeros_(p1 - j,1)],axis=0))))).T
    elif (whichmodel == 3):
#    build Regression Lagrange polynomials (p1 >= q)
#    (QZ and RZ are the factors of Z = M)
        Lj=(pinv_(RZ).dot( QZ.T.dot( concatenate_([zeros_(j - 1,1),array([[1]]),zeros_(p1 - j,1)],axis=0)))).T
    return Lj