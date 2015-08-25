# -*- coding: utf-8 -*-
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from bcdfo_evalL import bcdfo_evalL_
from bcdfo_replace_in_Y import bcdfo_replace_in_Y_
from copy import copy

def bcdfo_include_in_Y_(x=None,QZ_=None,RZ_=None,Y_=None,choice_set=None,poisedness_threshold=None,criterion=None,xbase_=None,whichmodel=None,succ=None,scale_=None,shift_Y=None,Delta=None,normgx=None,kappa_ill=None,*args,**kwargs):

    """
#  Attempts to include x in the interpolation set by replacing an existing
#  interpolation point. This point is chosen in the set defined by the
#  intersection of choice_set and the set of indices such that the absolute
#  value of the associated polynomial at x is larger than
#  poisedness_threshold. If more than one point remains, the choice is then
#  made my maximizing a criterion.  If criterion is 'weighted', then one
#  chooses to maximize ||Y(:,j)-x||^2 |L_j(x)|.  If criterion is 'standard',
#  then one chooses to maximize |L_j(x)|, while the value of ||Y(:,j)-Y(:,1)||
#  alone is used when criterion is 'distance'. (Ties are broken arbitrarily). 
#  If no point satisfies the desired properties, nothing is done.
#
#  INPUT:
#
#  x           : the new vector one attempts to include in Y
#  QZ          : the Q matrix of the QR decomposition of Z(Y)
#  RZ          : the R matrix of the QR decomposition of Z(Y)
#  Y           : a matrix whose columns contain the current interpolation points
#  choice_set  : the indices (of Y's columns) amongst which the replaced point 
#                may be chosen
#  poisedness_threshold : the minimal acceptable improvement in volume simplex
#                (poisedness) for the new point to be acceptable. 
#  criterion   : either 'weighted' or 'standard' (see above)
#  xbase       : the current base point
#  whichmodel  : kind of model/Lagrange polynomial to compute
#  succ        : indicates a successful iteration
#  scale       : the current interpolation set scaling
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#  Delta       : trust-region radius
#  normgx      : infinity norm of the projected gradient
#  kappa_ill   : threshold to declare a system matrix as ill-conditioned
#
#  OUTPUT:
#
#  QZ          : the Q matrix of the QR decomposition of Z(Y) for the updated Y
#  RZ          : the R matrix of the QR decomposition of Z(Y) for the updated Y
#  Y           : the updated set of interpolation points (column-wise)
#  pos         : the position of the new point, if included, 0 otherwise
#  xbase       : the updated base point
#  scale       : the updated interpolation set scaling
#
#  PROGRAMMING: Ph. Toint, S. Gratton, April 2009. (This version 22 VI 2009)
#
#  DEPENDENCIES: bcdfo_evalP, bcdfo_replace_in_Y.
#
#  TEST:
#  Y = [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; whichmodel = 0;
#  [ QZ, RZ, xbase, scale] = bcdfo_build_QR_of_Y( Y, whichmodel, 1, 1, 1, 1e15 );
#  [QZplus, RZplus, Yplus, pos, xbase, scale] = ...
#     bcdfo_include_in_Y([-1;1], QZ, RZ, Y, [2:6], 0.01, ...
#     'weighted', xbase, whichmodel, 0, scale, 1, 1, 1, 1e15 )
#  QZplus =
#
#   -1.0000         0         0         0         0         0
#         0    0.9701    0.0000   -0.2425    0.0000   -0.0000
#         0    0.0000   -0.9701   -0.0000    0.0000   -0.2425
#         0    0.2425   -0.0000    0.9701   -0.0000    0.0000
#         0    0.0000   -0.2425    0.0000   -0.0000    0.9701
#         0   -0.0000    0.0000    0.0000    1.0000    0.0000
#
#  RZplus =
#
#   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000
#         0    0.5154         0    1.0914   -0.4548    0.0000
#         0         0   -0.5154   -0.0000   -0.5154   -1.0914
#         0         0         0    0.2425    0.2425   -0.0000
#         0         0         0         0   -0.2500   -0.0000
#         0         0         0         0         0    0.2425
#
#  Yplus =
#
#     0     1     0     2    -1     0
#     0     0     1     0     1     2
#
#  pos =
#
#     5
#
#  xbase =
#
#     0
#     0
#
#  scale =
#
#    1.0000
#    0.5000
#    0.5000
#    0.2500
#    0.2500
#    0.2500
#
#  and
#  RZplus\(QZplus'*bcdfo_evalZ((Yplus-Y(:,1)*ones(1,6))*scale(2),6))
#  should give the identity matrix.
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
    """  
    
#    varargin = cellarray(args)
#    nargin = 15-[x,QZ,RZ,Y,choice_set,poisedness_threshold,criterion,xbase,whichmodel,succ,scale,shift_Y,Delta,normgx,kappa_ill].count(None)+len(args)
    Y=copy(Y_)
    QZ=copy(QZ_)
    RZ=copy(RZ_)
    xbase=copy(xbase_)
    scale=copy(scale_)

    if (length_(choice_set) == 0):
        pos=-1
        return QZ,RZ,Y,pos,xbase,scale
    Lvals=bcdfo_evalL_(QZ,RZ,Y,choice_set,x,xbase,whichmodel,scale,shift_Y)
    choice=find_(abs(Lvals) > poisedness_threshold)
    lc=length_(choice)
    if (lc == 0):
        pos=-1
        return QZ,RZ,Y,pos,xbase,scale
    crit_val=0
    pos=-1
    for i in range(0,lc):
        j=choice[i] # since choice is a vector column, choice[i] will return something like array([a]), and dimensions below are correct
        if (criterion == 'weighted'):
            if (succ == 1):
                cv=norm_(Y[:,j] - x) ** 2 * abs(Lvals[j])
            else:
                cv=norm_(Y[:,j] - Y[:,[0]]) ** 2 * abs(Lvals[j])
        else:
            if (criterion == 'standard'):
                cv=abs(Lvals[j])
            else:
                if (criterion == 'distance'):
                    if (succ == 1):
                        cv=norm_(Y[:,j] - x)
                    else:
                        cv=norm_(Y[:,j] - Y[:,[0]])
        if (cv > crit_val):
            pos=copy(j)
            crit_val=copy(cv)
    if (int(pos) == -1):
        return QZ,RZ,Y,int(pos),xbase,scale
    QZ,RZ,Y,xbase,scale=bcdfo_replace_in_Y_(QZ,RZ,x,Y,int(pos),xbase,whichmodel,scale,shift_Y,Delta,normgx,kappa_ill,nargout=5)
    return QZ,RZ,Y,int(pos),xbase,scale