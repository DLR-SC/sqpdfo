# -*- coding: utf-8 -*-
from runtime import isempty_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from copy import copy

def sqpdfo_swap_in_Y_(i=None,j=None,QZ=None,RZ=None,Y_=None,ind_Y_=None,fY_=None,ciY_=None,ceY_=None,xbase=None,whichmodel=None,scale=None,shift_Y=None,Delta=None,normgx=None,kappa_ill=None,*args,**kwargs):
    """
#
#  Swaps the position of interpolation points i and j in Y, and updates the
#  factorization of Z(Y) accordingly.
#
#  INPUT:
#
#  i           : the position in which Y(:,j) should be inserted
#  j           : the position in which Y(:,i) should be inserted
#  QZ, RZ      : the QR decomposition of Z(Y)
#  Y           : a matrix whose columns contain the current interpolation points
#  ind_Y		   : set of indices which points are in Y out of set X
#  fY				: corresponding function values
#  xbase       : the current base point
#  whichmodel  : kind of model/Lagrange polynomial to compute
#  scale       : the current interpolation set scaling
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#  Delta       : trust-region radius
#  normgx      : infinity norm of the projected gradient
#  kappa_ill   : threshold to declare a system matrix as ill-conditioned
#
#  OUTPUT:
#
#  QZ          : the updated Q matrix of the QR decomposition of Z(Y)
#  RZ          : the updated R matrix of the QR decomposition of Z(Y)
#  Y           : the updated set of interpolation points (column-wise)
#  xbase       : the base point after the exchange
#  scale       : the interpolation set scaling after the exchange
#
#  PROGRAMMING: Ph. Toint and A. Troeltzsch, January 2009. 
#               A. Troeltzsch, Feb 2013.
#
#  DEPENDENCIES: bcdfo_build_QR_of_Y
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
    """
#    varargin = cellarray(args)
#    nargin = 16-[i,j,QZ,RZ,Y,ind_Y,fY,ciY,ceY,xbase,whichmodel,scale,shift_Y,Delta,normgx,kappa_ill].count(None)+len(args)

    Y=copy(Y_)
    ind_Y=copy(ind_Y_)
    fY=copy(fY_)
    ciY=copy(ciY_)
    ceY=copy(ceY_)

    
#  Ensure that ii is smaller than jj.


    if (i > j):
        ii=copy(j)
        jj=copy(i)
    else:
        if (i < j):
            ii=copy(i)
            jj=copy(j)
        else:
            return QZ,RZ,Y,ind_Y,fY,ciY,ceY,xbase,scale
#  Permute the columns of Y, the indices in ind_Y and the values in fY.

    y=copy(Y[:,ii])
    Y[:,ii]=copy(Y[:,jj])
    Y[:,jj]=y
    ind=copy(ind_Y[ii])
    ind_Y[ii]=copy(ind_Y[jj])
    ind_Y[jj]=ind
    f=fY[ii]
    fY[ii]=copy(fY[jj])
    fY[jj]=f
    if not isempty_(ciY):
        ci=copy(ciY[:,ii])
        ciY[:,ii]=copy(ciY[:,jj])
        ciY[:,jj]=ci
    if not isempty_(ceY):
        ce=copy(ceY[:,ii])
        ceY[:,ii]=copy(ceY[:,jj])
        ceY[:,jj]=ce
#  a new factorization must be computed

    QZ,RZ,xbase,scale=bcdfo_build_QR_of_Y_(Y,whichmodel,shift_Y,Delta,normgx,kappa_ill,nargout=4)
    return QZ,RZ,Y,ind_Y,fY,ciY,ceY,xbase,scale
