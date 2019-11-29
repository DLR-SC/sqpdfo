# -*- coding: utf-8 -*-

from sqpdfo.runtime import *
from numpy import inf,diag, isnan, isinf, nan_to_num, linalg


def sqpdfo_check_cond_(A=None,cthreshold=None,options=None,*args,**kwargs):
    """
    
#
# checks the condition of matrix A and if a threshold for the condition
# number is exceeded, matrix is perturbed by exchanging very small singular
# values
#
# programming: A. Troeltzsch, 2014 (nan_to_num: 06/2019)
#
    """

    badcond=0
    eps=1e-14
    
    # check whether A has Nan or Inf entries
    if (isnan(A).any() or isinf(A).any()):
        if options.verbose >= 3:
            print('### sqpdfo_check_cond: Matrix A has Nan or Inf entries !!!')
        
        # repair matrix
        A = nan_to_num(A)
    
    # compute condition number
    condA=cond_(A)
    if (condA > cthreshold):
        badcond=1
    
    if (badcond):
        U,Sdiag,V=linalg.svd(A)
        V=V.T
        a=Sdiag < 1e-07
        Sdiag[a]=1e-07
        S=diag(Sdiag)
        A=(V.dot(S.dot(U.T))).T
      # Make sure it is symmetric
        
        if norm_(A - A.T,inf) > eps:
            if options.verbose >= 3:
                disp_('### sqpdfo_check_cond: ',"matrix is non symmetric. Resetting A=(A+A')/2.")
            A=0.5*(A + A.conj().T)
    return A,badcond
