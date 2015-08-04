# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 13:09:54 2015

@author: lien_ol
"""
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from bcdfo_evalZ import *
from numpy import copy
    
def bcdfo_build_QR_of_Y_(Y_=None,whichmodel=None,shift_Y=None,Delta=None,normgx=None,kappa_ill=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 6-[Y,whichmodel,shift_Y,Delta,normgx,kappa_ill].count(None)+len(args)

    #Copy only if necessary
    if (shift_Y):
        Y=copy_(Y_)
    else:
        Y=Y_
        
    n,p1=size_(Y,nargout=2)
    badcond=0
    if (normgx == 0.0):
        _del=Delta ** 2 / 1e-12
    else:
        _del=Delta ** 2 / normgx
    if (_del > 0.1):
        _del=0.1
    if (_del < 1e-10):
        _del=1e-10
    if (whichmodel == 0):
        q=copy_(p1)
    else:
        q=((n + 1) * (n + 2)) / 2
    if (whichmodel == 3 and p1 < q):
        whichmodel=2
            
    if (shift_Y and (p1 > 1)):
        xbase=Y[:,1]
        scaleY=0
        for i in arange_(1,p1).reshape(-1):
            Y[:,i]=Y[:,i] - xbase
            scaleY=max_(scaleY,norm_(Y[:,i]))
        scale=concatenate_([matlabarray([[1]]),scaleY ** - 1 * ones_(1,min_(n,q - 1)),scaleY ** - 2 * ones_(1,q - n - 1)], axis=1).T
        Y=Y / scaleY
    else:
        scale=ones_(q,1)
        xbase=zeros_(size_(Y,1),1)
    if (whichmodel == 0):
        Z=bcdfo_evalZ_(Y,q)
        if (length_(find_(isnan_(Z))) == 0 and length_(find_(isinf_(Z))) == 0):
            condZ=cond_(Z)
            if (condZ > kappa_ill):
                badcond=1
        else:
            badcond=1
        if (badcond):
            U,S,V=svd_(Z,0,nargout=3)
            Sdiag=diag_(S)
            indices=find_(Sdiag < _del)
            Sdiag[indices]=_del
            S=diag_(Sdiag)
            M=(V * S * U.T).T
            QZ,RZ=qr_(M,nargout=2)
        else:
            QZ,RZ=qr_(Z,nargout=2)
    else:
        if (whichmodel == 1):
            if (p1 == n + 1 or p1 == q):
                F=bcdfo_evalZ_(Y,p1).T
            else:
                M=bcdfo_evalZ_(Y,q).T
                ML=M[:,1:n + 1]
                MQ=M[:,n + 2:q]
                F=matlabarray([[MQ * MQ.T,ML],[ML.T,zeros_(n + 1,n + 1)]])
            if (length_(find_(isnan_(F))) == 0 and length_(find_(isinf_(F))) == 0):
                condZ=cond_(F)
                if (condZ > kappa_ill):
                    badcond=1
            else:
                badcond=1
            if (badcond):
                U,S,V=svd_(F,0,nargout=3)
                Sdiag=diag_(S)
                indices=find_(Sdiag < _del)
                Sdiag[indices]=_del
                S=diag_(Sdiag)
                M=(V * S.T * U.T).T
                QZ,RZ=qr_(M,nargout=2)
            else:
                QZ,RZ=qr_(F,nargout=2)
        else:
            if (whichmodel == 2):
                Z=bcdfo_evalZ_(Y,q)
                if (length_(find_(isnan_(Z))) == 0 and length_(find_(isinf_(Z))) == 0):
                    condZ=cond_(Z)
                    if (condZ > kappa_ill):
                        badcond=1
                else:
                    badcond=1
                if (badcond):
                    U,S,V=svd_(Z,'econ',nargout=3)
                    Sdiag=diag_(S)
                    indices=find_(Sdiag < _del)
                    Sdiag[indices]=_del
                    S=diag_(Sdiag)
                    M=(V * S.T * U.T).T
                    QZ,RZ=qr_(M,nargout=2)
                else:
                    QZ,RZ=qr_(Z,nargout=2)
            else:
                if (whichmodel == 3):
                    Z=bcdfo_evalZ_(Y,q).T
                    if (length_(find_(isnan_(Z))) == 0 and length_(find_(isinf_(Z))) == 0):
                        condZ=cond_(Z)
                        if (condZ > kappa_ill):
                            badcond=1
                    else:
                        badcond=1
                    if (badcond):
                        U,S,V=svd_(Z,0,nargout=3)
                        Sdiag=diag_(S)
                        indices=find_(Sdiag < _del)
                        Sdiag[indices]=_del
                        S=diag_(Sdiag)
                        M=(V * S * U.T).T
                        QZ,RZ=qr_(M,nargout=2)
                    else:
                        QZ,RZ=qr_(Z,nargout=2)
    return QZ,RZ,xbase,scale