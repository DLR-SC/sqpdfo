# -*- coding: utf-8 -*-
from __future__ import division
#try:
from runtime import *
#print "type eig", type(eig_)
#except ImportError:
#    from smop.runtime import *
from numpy import diag, array,isreal,real
from copy import copy
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
#    varargin = cellarray(args)
#    nargin = 2-[A,options].count(None)+len(args)

    ev=numpy.linalg.eigvals(A)
    evneg=ev[ev < 0]
    if not isempty_(evneg):
        ZERO=array([1e-10])
        EPS=array([1e-09])
        d,v=numpy.linalg.eig(A)   # CALCULATE EIGENVECTOR AND EIGENVALUES 
        d[d < ZERO]=EPS  # FIND ALL EIGENVALUES<=ZERO AND CHANGE THEM FOR EPS 
        d=diag(d)    # CONVERT VECTOR d INTO A MATRIX 
        A=v.dot( d .dot( v.T ))  # RECOMPOSE MATRIX USING EIGENDECOMPOSITION 
    # check for complex entries
       
        if not isempty_(find_(~ isreal(A),1)):
            if options.verbose >= 3:
                #This is a little bit weird since we did not test if it was symmetric. Maybe a error of copy/paste with ecdfo_check_cond ?
                disp_('### ecdfo_check_convex: matrix is non symmetric. Resetting A.')
            A=0.5*(A+A.conj().T)

    #On account of the test asking if there is a non real in A, I suppose that we only want the real parts (also having imaginary parts raise errors on problem 5 in ecdfo_func),
    #hence the real conversion
    return real(A)
