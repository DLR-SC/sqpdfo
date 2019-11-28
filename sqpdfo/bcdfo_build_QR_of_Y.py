# -*- coding: utf-8 -*-

from runtime import *
from bcdfo_evalZ import *
from copy import copy
from numpy import *
    
def bcdfo_build_QR_of_Y_(Y_=None,whichmodel=None,shift_Y=None,Delta=None,normgx=None,kappa_ill=None,*args,**kwargs):
    
    """
#
#  Computes the QR factorization of the (possibly shifted) matrix containing
#  the polynomial expansion of the interpolation points. If shifting is
#  required, (i.e. if shift_Y is true) the matrix being factorized has columns
#  containing (Y(:,j)-Y(:,1))/ scale, where scale is the max( norm(Y(:,j)-Y(:,1)).
#
#  INPUT:
#
#  Y           : a matrix whose columns contain the current interpolation points
#  whichmodel  : kind of model to build
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#  Delta       : trust-region radius
#  normgx      : infinity norm of the projected gradient
#  kappa_ill   : threshold to declare a system matrix as ill-conditioned
#
#  OUTPUT:
#
#  QZ, RZ      : the QR factors of the (possibly shifted) matrix containing
#                the polynomial expansion of the interpolation points,
#  xbase       : the base point,
#  scale       : the model diagonal scaling.
#
#  PROGRAMMING: A. Troeltzsch, Ph. Toint, S. Gratton, 2009-2011.
#               (This version 14 I 2011)
#
#  DEPENDENCIES: bcdfo_evalZ
#
#  TEST:
#  Y = [ 1 2 1 3 3 1 ; 1 2 2 1 1.01 3 ];
#  [QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y, 0, 1, 1, 1, 1e15 );
#  model = ( QZ * ( RZ' \ [1 2 3 4 5 6 ]' ) )';
#  model * bcdfo_evalZ( ([1;3]-xbase)*scale(2),6)
#  should give 6.0
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
    """
#    varargin = cellarray(args)
#    nargin = 6-[Y,whichmodel,shift_Y,Delta,normgx,kappa_ill].count(None)+len(args)

    #Copy only if necessary
    if (shift_Y):
        Y=copy(Y_)
    else:
        Y=Y_
        
    n,p1=size_(Y,nargout=2)
    badcond=0

#  Compute and check the size of the tolerance.
    
    
    if (normgx == 0.0):
        _del=Delta ** 2 / 1e-12
    else:
        _del=Delta ** 2 / normgx
    if (_del > 0.1):
        _del=0.1
    if (_del < 1e-10):
        _del=1e-10

#  Define the size of the polynomial        
    
    if (whichmodel == 0):
        q=copy(p1)
    else:
        q=int(((n + 1) * (n + 2)) / 2)
    if (whichmodel == 3 and p1 < q):
#         for underdetermined regression model use min l2-norm model

        whichmodel=2
#  If shifting is active, the matrix of the interpoaltion points is first
#  shifted and scaled, and the base point and scaling factors defined.
    if (shift_Y and (p1 > 1)):
        xbase=copy(Y[:,0])
        scaleY=0
        for i in range(0,p1):
            Y[:,i]=Y[:,i] - xbase
            scaleY=max_(scaleY,norm_(Y[:,i]))
        scale=concatenate_([array([[1]]),scaleY ** - 1 * ones_(1,min_(n,q - 1)),scaleY ** - 2 * ones_(1,q - n - 1)], axis=1).T
        Y=Y / scaleY
    #  Otherwise, the base point and scaling factors are trivial.
    else:
        scale=ones_(q,1)
        xbase=zeros_(size_(Y,1),1)
        

#  Perform factorization of the Z matrix corresponding to the (possibly
#  shifted and scaled) interpolation matrix.
    if (whichmodel == 0):
#         QR of (sub-basis) interpolation matrix (p1 = q) or
#         (factorize matrix Z = M')
        
        Z=bcdfo_evalZ_(Y,q)
        
#         Check condition of Z and cure if ill-conditioned

        if (length_(find_(isnan(Z))) == 0 and length_(find_(isinf(Z))) == 0):
            condZ=cond_(Z)
            if (condZ > kappa_ill):
                badcond=1
        else:
            badcond=1
        if (badcond):
            U,Sdiag,V=linalg.svd(Z,full_matrices=0)
            indices=find_(Sdiag < _del)
            Sdiag[indices]=_del
            S=diag(Sdiag)
            M = np.dot(U, np.dot(S, V))
            QZ,RZ=qr_(M,nargout=2)
        else:
            QZ,RZ=qr_(Z,nargout=2)

    elif whichmodel == 1:
        # Mixed model: minimum Frobenius-norm model (when underdetermined) and
        # minimum l2-norm model (at linear and quadratic degree)
        if (p1 == n + 1 or p1 == q):
            F=bcdfo_evalZ_(Y,p1).T
        else:
            # QR of Minimum Frobenius norm interpolation matrix (p1 <= q)
            # (factorize matrix F = [MQMQ' ML; ML' 0])
            M=bcdfo_evalZ_(Y,q).T
            ML=M[:,0:n + 1]
            MQ=M[:,n + 1:q]
            F1=concatenate([MQ.dot(MQ.T),ML], axis=1)
            F2=concatenate([ML.T,zeros_(n + 1,n + 1)], axis=1)
            F=concatenate([F1,F2],axis=0)
        # Check condition of Z and cure if ill-conditioned
        if (length_(find_(isnan(F))) == 0 and length_(find_(isinf(F))) == 0):
            condZ=cond_(F)
            if (condZ > kappa_ill):
                badcond=1
        else:
            badcond=1
        if (badcond):
            U,Sdiag,V=linalg.svd(F,full_matrices=0)
            indices=find_(Sdiag < _del)
            Sdiag[indices]=_del
            S=diag(Sdiag)
            M = np.dot(U, np.dot(S, V))
            QZ,RZ=qr_(M,nargout=2)
        else:
            QZ,RZ=qr_(F,nargout=2)
    elif (whichmodel == 2):
#    QR of Minimum L2 norm interpolation matrix (p1 < q)
#    (factorize matrix Z = M')
        Z=bcdfo_evalZ_(Y,q)
#     Check condition of Z and cure if ill-conditioned
        if (length_(find_(isnan(Z))) == 0 and length_(find_(isinf(Z))) == 0):
            condZ=cond_(Z)
            if (condZ > kappa_ill):
                badcond=1
        else:
            badcond=1
        if (badcond):
            U,Sdiag,V=linalg.svd(Z,full_matrices=0)
            indices=find_(Sdiag < _del)
            Sdiag[indices]=_del
            S=diag(Sdiag)
            M = np.dot(U, np.dot(S, V))
            QZ,RZ=qr_(M,nargout=2)
        else:
            QZ,RZ=qr_(Z,nargout=2)
    elif (whichmodel == 3):
#    QR of Regression interpolation matrix (p1 >= q)
#    (factorize matrix Z = M)
        Z=bcdfo_evalZ_(Y,q).T
#     Check condition of Z and cure if ill-conditioned
        if (length_(find_(isnan(Z))) == 0 and length_(find_(isinf(Z))) == 0):
            condZ=cond_(Z)
            if (condZ > kappa_ill):
                badcond=1
        else:
            badcond=1
        if (badcond):
            U,Sdiag,V=linalg.svd(Z,full_matrices=0)
            indices=find_(Sdiag < _del)
            Sdiag[indices]=_del
            S=diag(Sdiag)
            M = np.dot(U, np.dot(S, V))
            QZ,RZ=qr_(M,nargout=2)
        else:
            QZ,RZ=qr_(Z,nargout=2)
    return QZ,RZ,xbase.reshape(-1,1),scale