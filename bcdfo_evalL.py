# -*- coding: utf-8 -*-
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
    
from bcdfo_evalP import bcdfo_evalP_
import numpy
from copy import copy
    
def bcdfo_evalL_(QZ=None,RZ=None,Y_=None,choice_set=None,x=None,xbase=None,whichmodel=None,scale=None,shift_Y=None,*args,**kwargs):
    """
#
#  Computes the values of the Lagrange polynomials at x.
#
#  INPUTS:
#
#  QZ, RZ      : the QR factors of the (possibly shifted) matrix Z containing
#                the polynomial expansion of the interpolation points
#  Y           : current interpolation set Y
#  choice_set  : the indices (of Y's columns) for which to calculate the 
#                Lagrange polynomial values
#  x           : the point at which the model must be evaluated
#  xbase       : the current base point
#  whichmodel  : kind of model/Lagrange polynomial to compute
#  scale       : the current model scaling
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#
#  OUTPUT:
#
#  values      : the values of the Lagrange polynomials at x.
#
#  PROGRAMMING: A. Troeltzsch, September 2010. 
#
#  DEPENDENCIES: bcdfo_evalP
#
#  TEST:
#  Y = [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; 
#  [ QZ, RZ, xbase, scale] = bcdfo_build_QR_of_Y( Y, 0, 1 );
#  values = bcdfo_evalL( QZ, RZ, Y, [2:6], [-1;1], xbase, 0, scale, 1 )
#  should give 
#     values =
#         0   97.0000    2.9900    1.0000 -100.0000   -0.4950 
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
    """
#    varargin = cellarray(args)
#    nargin = 9-[QZ,RZ,Y,choice_set,x,xbase,whichmodel,scale,shift_Y].count(None)+len(args)

    #Copy only if necessary
    if (shift_Y):
        Y=copy(Y_)
    else:
        Y=Y_
    
    choice_set=choice_set.reshape(-1)
    n,p1=size_(Y,nargout=2)
    I=eye_(p1)
    lc=length_(choice_set)
    q=((n + 1) * (n + 2)) / 2
    values=zeros_(p1,1)
    if (whichmodel == 3 and p1 < q):
#     for underdetermined regression model use min l2-norm model
        whichmodel=2
    if (whichmodel == 0):
#    evaluate (sub-basis) Lagrange polynomial (p1 = q)
        values[choice_set]=bcdfo_evalP_(I[choice_set,:].dot( (numpy.linalg.solve(RZ,QZ.T))),x,xbase,scale,shift_Y)
    elif (whichmodel == 1):
        if (p1 == n + 1 or p1 == q):
#        evaluate Minimum L2 norm Lagrange polynomial (p1 == n+1 or p1 == q)

            values[choice_set]=bcdfo_evalP_(I[choice_set,:] .dot( (QZ / RZ.T)),x,xbase,scale,shift_Y)
        else:
#       evaluate Minimum Frobenius norm Lagrange polynomial (p1 <= q)
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
            M=bcdfo_evalZ_(Y,q).T
            MQ=M[:,n + 1:q]
            if (shift_Y):
                phi=bcdfo_evalZ_((x - xbase) * scale[1],q)
            else:
                phi=bcdfo_evalZ_(x,q)
            values[choice_set]=concatenate_([  I[choice_set,:],  zeros_(lc,n + 1).dot(QZ .dot(numpy.linalg.solve(RZ.T,concatenate_([MQ.dot(phi[n + 1:q]),phi[0:n + 1]]))))],axis=1)
    elif (whichmodel == 2):
#    evaluate Minimum L2 norm Lagrange polynomial (p1 <= q)
        if (p1 < q):
            values[choice_set]=bcdfo_evalP_(I[choice_set,:].dot( (pinv_(RZ).dot(QZ.T))),x,xbase,scale,shift_Y)
        else:
            values[choice_set]=bcdfo_evalP_(I[choice_set,:].dot( (numpy.linalg.solve(RZ,QZ.T))),x,xbase,scale,shift_Y)
    elif (whichmodel == 3):
#    evaluate Regression Lagrange polynomial (p1 >= q)
        values[choice_set]=bcdfo_evalP_(I[choice_set,:].dot((QZ .dot(pinv_(RZ.T)))),x,xbase,scale,shift_Y)
    return values