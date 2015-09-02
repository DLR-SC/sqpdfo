# -*- coding: utf-8 -*-
from __future__ import division
#try:
from runtime import *
from numpy import inf,diag, isnan, isinf
from copy import copy
import numpy
#except ImportError:
#    from smop.runtime import *

def ecdfo_check_cond_(A_=None,cthreshold=None,options=None,*args,**kwargs):
    """
    
#
# checks the condition of matrix A and if a threshold for the condition
# number is exceeded, matrix is perturbed by exchanging very small singular
# values
#
# programming: A. Troeltzsch, 2014
#
    """
#    varargin = cellarray(args)
#    nargin = 3-[A,cthreshold,options].count(None)+len(args)

    A=copy(A_)

    badcond=0
    eps=1e-14
    if (isempty_(find_(isnan(A))) and isempty_(find_(isinf(A)))):
        condA=cond_(A)
        if (condA > cthreshold):
            badcond=1
    else:
        badcond=1
    if (badcond):
        U,Sdiag,V=numpy.linalg.svd(A)
        V=V.T
        a=Sdiag < 1e-07
        Sdiag[a]=1e-07
        S=diag(Sdiag)
        A=(V.dot(S.dot(U.T))).T
      # Make sure it is symmetric
        
        if norm_(A - A.T,inf) > eps:
            if options.verbose >= 3:
                disp_('### ecdfo_check_cond: ',"matrix is non symmetric. Resetting A=(A+A')/2.")
            A=0.5*(A + A.conj().T)
    return A,badcond
