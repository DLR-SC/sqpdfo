# -*- coding: utf-8 -*-
from __future__ import division
from runtime import *

from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from copy import copy

def bcdfo_replace_in_Y_(QZ=None,RZ=None,ynew=None,Y_=None,j=None,xbase=None,whichmodel=None,scale=None,shift_Y=None,Delta=None,normgx=None,kappa_ill=None,*args,**kwargs):
    """
#
#  Updates the interpolation set for a transformation of Y in Yplus where the
#  vector ynew replaces Y(:,j).  Also update the factorization of Z(Y)
#  accordingly.
#
#  INPUT:
#
#  QZ          : the Q matrix of the QR decomposition of Z(Y)
#  RZ          : the R matrix of the QR decomposition of Z(Y)
#  Y           : the current interpolation set (by columns)
#  ynew        : the vector hich is to replace Y(:,j) in the interpolation set
#  j           : the index of the interpolation point to be replaced
#  xbase       : the current base point
#  whichmodel  : kind of model to build
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
#  Y           : the updated Y
#  xbase       : the base point after the replacement
#  scale       : the interpolation set scaling after the replacement
#
#  PROGRAMMING: Ph. Toint, S. Gratton, A. Troeltzsch April 2009.
#               (This version 15 IX 2010)
#
#  USES: bcdfo_build_QR_of_Y
#
#  TEST:
#  Y = [ 1 1 0 0 0 -1; 1 0 -1 1 0 0 ];
#  [QZ,RZ,xbase,scale] = bcdfo_build_QR_of_Y( Y , 0, 0, 1, 1, 1e15 );
#  ynew = 0.25*Y(:,3);
#  [ QZplus, RZplus, Yplus ] = bcdfo_replace_in_Y( QZ, RZ, ynew, Y, 3, ...
#            xbase, whichmodel, scale, 0, 1, 1, 1e15)
#  QZplus =
#
#   -0.4714   -0.4714    0.7205    0.1527    0.1149    0.0000
#   -0.4714   -0.4714   -0.5764   -0.1221   -0.0919    0.4472
#   -0.4714    0.4714   -0.1891    0.6333    0.3446    0.0000
#   -0.2357   -0.2357   -0.2882   -0.0611   -0.0459   -0.8944
#   -0.2357    0.2357    0.1081    0.1813   -0.9189   -0.0000
#   -0.4714    0.4714    0.1351   -0.7240    0.1149   -0.0000
#
#  RZplus =
#
#   -2.1213   -1.0607   -0.3609   -1.0607   -0.4714   -0.1179
#         0   -1.0607   -0.5819    0.1179   -0.4714   -0.1179
#         0         0    0.7711    0.5854    0.7205    1.1527
#         0         0         0    0.8766    0.1527    0.2442
#         0         0         0         0    0.1149    0.1838
#         0         0         0         0         0   -0.8944
#
#  Yplus =
#
#    1.0000    1.0000         0         0         0   -1.0000
#    1.0000         0   -0.2500    1.0000         0         0
#
#  The scaled version:
#  Y = [ 1 1 0 0 0 -1; 1 0 -1 1 0 0 ];
#  [QZ,RZ,xbase,scale] = bcdfo_build_QR_of_Y( Y , 0, 1, 1, 1, 1e15 );
#  ynew = 0.25*Y(:,3);
#  [ QZplus, RZplus, Yplus ] = bcdfo_replace_in_Y( QZ, RZ, ynew, Y, 3, ...
#            xbase, whichmodel, scale, 1, 1, 1, 1e15 )
#  QZplus =
#
#   -1.0000         0         0         0         0         0
#         0    0.0000   -0.8552    0.4700    0.0000    0.2182
#         0   -0.9759    0.0127    0.0232   -0.2166    0.0000
#         0   -0.0000    0.1912   -0.1051    0.0000    0.9759
#         0    0.2182    0.0569    0.1036   -0.9687    0.0000
#         0   -0.0000    0.4781    0.8699    0.1211    0.0000
#
#  RZplus =
#
#   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000
#         0    0.4583    0.5796   -0.0000    0.4583    0.4583
#         0         0    0.5229    0.4016    0.4972    1.0327
#         0         0         0   -0.2207   -0.0467   -0.1145
#         0         0         0         0    0.0242    0.0484
#         0         0         0         0         0    0.1952
#
#  Yplus =
#
#    1.0000    1.0000         0         0         0   -1.0000
#    1.0000         0   -0.2500    1.0000         0         0
#
    """  
    
#    varargin = cellarray(args)
#    nargin = 12-[QZ,RZ,ynew,Y,j,xbase,whichmodel,scale,shift_Y,Delta,normgx,kappa_ill].count(None)+len(args)


    Y=copy(Y_)

# replace new point in Y
    Y[:,j]=ynew.reshape(-1)
#  a new factorization must be computed. 

    QZ,RZ,xbase,scale=bcdfo_build_QR_of_Y_(Y,whichmodel,shift_Y,Delta,normgx,kappa_ill,nargout=4)
    return QZ,RZ,Y,xbase,scale