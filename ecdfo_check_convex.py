# -*- coding: utf-8 -*-

from runtime import *
from numpy import diag, array,isreal,real, isnan, isinf, linalg, nan_to_num

def ecdfo_check_convex_(A=None,options=None,*args,**kwargs):
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

    # check whether A has Nan or Inf entries
    if (isnan(A).any() or isinf(A).any()):
        if options.verbose >= 3:
            print('### ecdfo_check_convex: Matrix A has Nan or Inf entries !!!')
               
    # compute eigenvalues
    ev=linalg.eigvals(A)
    
    # check whether A has negative Eigenvalues
    evneg=ev[ev < 0]
    
    # remove neg. eigvals by eigendecomposition
    if not isempty_(evneg):
        ZERO=array([1e-10])
        EPS=array([1e-09])
        d,v=linalg.eig(A)   # CALCULATE EIGENVECTOR AND EIGENVALUES 
        d[d < ZERO]=EPS  # FIND ALL EIGENVALUES<=ZERO AND CHANGE THEM FOR EPS 
        d=diag(d)    # CONVERT VECTOR d INTO A MATRIX 
        A=v.dot( d .dot( v.T ))  # RECOMPOSE MATRIX USING EIGENDECOMPOSITION

    # check for complex entries (non-symmetric matrix)
    if not isempty_(not isreal(A).any()):
        if options.verbose >= 3:
            print('### ecdfo_check_convex: Matrix A is non symmetric. Resetting A !')
            
        # repair matrix
        A = (A+A.T)*0.5
    
    # remove still possible complex entries       
    return real(A)
