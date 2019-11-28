# -*- coding: utf-8 -*-

from runtime import *
from copy import copy
from numpy import diag, array,isreal,real, isnan, isinf, linalg, nan_to_num

def sqpdfo_check_convex_(A_=None,options=None,*args,**kwargs):
    """
#
# Checks whether matrix A is not positive definite (has negative eigenvalues).
# If A has negative eigenvalues, it is made positive definite by perturbing small
# eigenvalues.
#
# The bigger the value EPS, the more different the matrix gets.
#
# programming: A. Troeltzsch, 2014 (nan_to_num: 06/2019)
#
    """
    A = copy(A_)

    # check whether A has Nan or Inf entries
    if (isnan(A).any() or isinf(A).any()):
        if options.verbose >= 2:
            print('### sqpdfo_check_convex: Matrix A has Nan or Inf entries !!!')
               
    # compute eigenvalues of symmetrc matrix A
    ev=linalg.eigvalsh(A)

    # check whether A has negative Eigenvalues
    evneg=ev[ev < 0]
    
    # remove neg. eigvals by eigendecomposition
    if not isempty_(evneg):
        #print('eigendecomposition')
        ZERO=array([1e-10])
        EPS=array([1e-09])
        d,v=linalg.eigh(A)   # CALCULATE EIGENVECTOR AND EIGENVALUES 
        d[d < ZERO]=EPS  # FIND ALL EIGENVALUES<=ZERO AND CHANGE THEM FOR EPS 
        d=diag(d)    # CONVERT VECTOR d INTO A MATRIX 
        A=v.dot( d .dot( v.T ))  # RECOMPOSE MATRIX USING EIGENDECOMPOSITION
    
    # return A
    return real(A)
