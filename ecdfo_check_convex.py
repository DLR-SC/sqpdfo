# -*- coding: utf-8 -*-

from runtime import *
from numpy import diag, array,isreal,real
import numpy

def ecdfo_check_convex_(A=None,options=None,*args,**kwargs):
    """
#
# checks whether matrix A is convex.
# if A is not convex, matrix is convexified by perturbing small
# eigenvalues.
#
# the bigger the value EPS, the more different the matrix gets.
#
# programming: A. Troeltzsch, 2014
#
    """

    ev=numpy.linalg.eigvals(A)
    evneg=ev[ev < 0]
    if not isempty_(evneg):
        ZERO=array([1e-10])
        EPS=array([1e-09])
        d,v=numpy.linalg.eig(A)   # CALCULATE EIGENVECTOR AND EIGENVALUES 
        d[d < ZERO]=EPS  # FIND ALL EIGENVALUES<=ZERO AND CHANGE THEM FOR EPS 
        d=diag(d)    # CONVERT VECTOR d INTO A MATRIX 
        A=v.dot( d .dot( v.T ))  # RECOMPOSE MATRIX USING EIGENDECOMPOSITION

    # remove possible complex entries
       
    return real(A)
