# -*- coding: utf-8 -*-
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from bcdfo_find_new_yj import bcdfo_find_new_yj_
from bcdfo_replace_in_Y import bcdfo_replace_in_Y_
from copy import copy
    
def bcdfo_repair_Y_(QZ_=None,RZ_=None,Y_=None,Delta=None,farfact=None,farthr=None,closethr=None,eps_L=None,xbase_=None,lSolver=None,whichmodel=None,hardcons=None,xl=None,xu=None,indfree=None,stratLam=None,scale_=None,shift_Y=None,normgx=None,kappa_ill=None,*args,**kwargs):
    """
#
#  Repairs the interpolation set Y by first replacing the interpolation points
#  at a distance larger that farfact*Delta from Y(:,0), and then replaces 
#  interpolation points with the best possible point obtained by taking the best
#  optimal replacement, until this improvement no longer causes a relative simplex
#  volume increase by at least threshold.  For conventions on how polynomials
#  are represented, see the documentation of evalZ.
#
#  INPUT:
#
#  QZ          : the Q matrix of the QR decomposition of Z(Y)
#  RZ          : the R matrix of the QR decomposition of Z(Y)
#  Y           : a matrix whose columns contain the current interpolation points
#  Delta       : the radius of the ball centered at Y(:,1) in which the
#                replacement point must be found
#  farfact     : multiple of TR radius defining far points
#  farthr      : the minimal acceptable improvement in volume simplex
#                (poisedness) for the new point to be acceptable in
#                replacement of a far away point.
#  closethr    : the minimal acceptable improvement in volume simplex
#                (poisedness) for the new point to be acceptable in
#                replacement of a close point.
#  eps_L       : the relative accuracy on the trust-region constraint for
#                maximization of the Lagrange polynomials
#  xbase       : the current base point
#  lSolver     : linear solver used for the minimization of the model
#  whichmodel  : kind of model/Lagrange polynomial to compute
#  scale       : the current interpolation set scaling
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#  normgx      : infinity norm of the projected gradient
#  kappa_ill   : threshold to declare a system matrix as ill-conditioned
#
#  OUTPUT:
#
#  QZ          : the Q matrix of the QR decomposition of Z(Y) for updated Y
#  RZ          : the R matrix of the QR decomposition of Z(Y) for updated Y
#  Y           : the updated set of interpolation points (column-wise)
#  replaced    : a vector containing the indices of the points of Y which have
#                been replaced
#  maximprove  : the best improvement still possible by replacing a single
#                interpolation point (gives an estimate of poisedness)
#  Y_radius    : the maximum distance from any interpolation point to the 
#                base point
#  xbase       : the updated base point
#  scale       : the updated interpolation set scaling
#
#  PROGRAMMING: Ph. Toint, S. Gratton, April 2009. (This version 22 VI 2009)
#
#  DEPENDENCIES: bcdfo_find_new_yj, bcdfo_gradP, bcdfo_hessP, bcdfo_replace_in_Y
#                bcdfo_poisedness_Y
#
#  TEST:
#  Y = [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ];
#  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y(  Y, 0, 0, 1, 1, 1e15 );
#  [ QZplus, RZplus, Yplus, replaced, maximprove, Y_radius, xbase, scale ] = ...
#            bcdfo_repair_Y( QZ, RZ, Y, 0.7, 10, 1.0e-10, 1.1, 0.001, xbase,...
#            1, 0, 0, [-10,-10], [10,10], [1,2], 1, scale, 0, 1, 1e15 )
#  QZplus =
#
#   -1.0000         0         0         0         0         0
#         0    0.6596    0.5831    0.2590   -0.3837   -0.1032
#         0    0.6401   -0.7030    0.1160    0.1723   -0.2299
#         0   -0.1657    0.0581    0.9007    0.3837    0.1032
#         0   -0.1560    0.2178   -0.0539    0.2743   -0.9220
#         0   -0.3216   -0.3391    0.3244   -0.7750   -0.2752
#
#  RZplus =
#
#   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000
#         0   -0.7616   -0.1207    0.9878   -0.1412    0.9682
#         0         0    0.7473    1.2823   -0.5045   -0.9705
#         0         0         0    2.3194   -0.0241    0.1242
#         0         0         0         0    0.5459    0.8931
#         0         0         0         0         0   -2.3039
#
#  Yplus =
#
#         0   -0.5023    0.3561    2.0000   -0.6030         0
#         0   -0.4875   -0.6026         0    0.3555    2.0000
#
#  replaced =
#
#     2     3     5
#
#  maximprove =
#
#    1.0002
# 
#  Y_radius =
#
#     2
#
#  xbase =
#
#     0
#     0
#
#  scale =
#
#     1
#     1
#     1
#     1
#     1
#     1
#
#  The scaled version: 
#  Y = [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ];
#  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y(  Y, 0, 1, 1, 1, 1e15 );
#  [ QZplus, RZplus, Yplus, replaced, maximprove, Y_radius, xbase, scale ] = ...
#           bcdfo_repair_Y( QZ, RZ, Y, 0.7, 10, 1.0e-10, 1.1, 0.001, xbase,...
#           1, 0, 0, [-10,-10], [10,10], [1,2], 1, scale, 1, 1, 1e15 )
#  QZplus =
#
#   -1.0000         0         0         0         0         0
#         0    0.7017    0.6633    0.1499   -0.2060   -0.0529
#         0    0.6810   -0.7144    0.0589    0.0925   -0.1178
#         0   -0.0881    0.0262    0.9003    0.4120    0.1058
#         0   -0.0830    0.1119   -0.0352    0.2945   -0.9448
#         0   -0.1710   -0.1910    0.4028   -0.8322   -0.2820
#
#
#  RZplus =
#
#   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000
#         0   -0.3580   -0.0762    0.6576   -0.0867    0.6395
#         0         0    0.3491    0.6764   -0.3138   -0.6584
#         0         0         0    0.6001   -0.0159    0.0413
#         0         0         0         0    0.1465    0.2397
#         0         0         0         0         0   -0.5902
#
#  Yplus =
#
#         0   -0.5023    0.3561    2.0000   -0.6030         0
#         0   -0.4875   -0.6026         0    0.3555    2.0000
#
#  replaced =
#
#     2     3     5
#
#  maximprove =
#
#    1.0002
#
#  Y_radius =
#
#     2
#
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
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
    """    
    
#    varargin = cellarray(args)
#    nargin = 20-[QZ,RZ,Y,Delta,farfact,farthr,closethr,eps_L,xbase,lSolver,whichmodel,hardcons,xl,xu,indfree,stratLam,scale,shift_Y,normgx,kappa_ill].count(None)+len(args)

    Y=copy(Y_)
    QZ=copy(QZ_)
    RZ=copy(RZ_)
    xbase=copy(xbase_)
    scale=copy(scale_)

    replaced=np.array([])
    verbose=0
    max_improve_loops=20
    if (verbose):
        disp_('--- enter bcdfo_repair_Y  for Delta = ',num2str_(Delta),' ---')
    n,p1=size_(Y,nargout=2)
    d=zeros_(p1,1)
    for j in range(1,p1): #it was initially in matlab range(2,p1)
        d[j]=norm_(Y[:,j] - Y[:,0])
#    dsorted,jsorted=sort_(d,'descend',nargout=2)
    dsorted=np.sort(d.reshape(-1))[::-1]
    jsorted=np.argsort(d.reshape(-1))[::-1]
    if (verbose):
      print 'd='
      print d
      print 'dsorted='
      print dsorted
      print 'jsorted='
      print jsorted
    for j in range(0,p1):
        if (dsorted[j] > farfact * (1 + eps_L) * Delta):
            jmax=jsorted[j]
            if (hardcons == 1):
                y,improvement,msgTR=bcdfo_find_new_yj_bc_(QZ,RZ,Y,jmax,Delta,eps_L,xbase,lSolver,whichmodel,xl,xu,indfree,stratLam,scale,shift_Y,nargout=2)
            else:
                y,improvement,msgTR=bcdfo_find_new_yj_(QZ,RZ,Y,jmax,Delta,eps_L,xbase,lSolver,whichmodel,scale,shift_Y,nargout=2)
            if (verbose):
                disp_('lambda = ',num2str_(improvement),' at j=',num2str_(jmax))
            QZ,RZ,Y,xbase,scale=bcdfo_replace_in_Y_(QZ,RZ,y,Y,jmax,xbase,whichmodel,scale,shift_Y,Delta,normgx,kappa_ill,nargout=5)
#            replaced=matlabarray([replaced,jmax])
            replaced  = np.append( replaced, jmax )
            d[jmax]=norm_(y - Y[:,[0]])
        else:
            Y_radius=dsorted[j]
            break
    if (verbose):
        print 'replaced='
        print replaced
        poised,Y_radius=bcdfo_poisedness_Y_(QZ,RZ,Y,eps_L,xbase,lSolver,whichmodel,hardcons,xl,xu,indfree,stratLam,scale,shift_Y,nargout=2)
        disp_('after distant replace: poisedness(Y) = ',num2str_(poised),' Y_radius  = ',num2str_(Y_radius))
    for k in range(1,max_improve_loops+1):
        maximprove=0
        for j in range(1,p1):
            if (hardcons == 1):
                y,improvement, msgTR=bcdfo_find_new_yj_bc_(QZ,RZ,Y,j,Delta,eps_L,xbase,lSolver,whichmodel,xl,xu,indfree,stratLam,scale,shift_Y,nargout=2)
            else:
                y,improvement, msgTR=bcdfo_find_new_yj_(QZ,RZ,Y,j,Delta,eps_L,xbase,lSolver,whichmodel,scale,shift_Y,nargout=2)
            if (verbose > 1):
                disp_(' ==> j = ',int2str_(j),' improve = ',num2str_(improvement))
                print 'y='
                print y
            if (improvement > maximprove):
                maximprove=copy(improvement)
                jmax=copy(j)
                ymax=copy(y)
        if (maximprove < closethr or jmax == 0):
            Y_radius=max_(d)
            if (verbose):
                print 'replaced='
                print replaced
                disp_('maximprove(small)= ',num2str_(maximprove),', jmax= ',int2str_(jmax))
                poised,Y_radius=bcdfo_poisedness_Y_(QZ,RZ,Y,eps_L,xbase,lSolver,whichmodel,hardcons,xl,xu,indfree,stratLam,scale,shift_Y,nargout=2)
                disp_('after everything: poisedness(Y) = ',num2str_(poised),' Y_radius  = ',num2str_(Y_radius))
                disp_('--- exit 1 bcdfo_repair_Y in round k=',num2str_(k),' ---')
            if (isempty_(replaced)):
                maximprove=0
            return QZ,RZ,Y,replaced,maximprove,Y_radius,xbase,scale
        if (verbose):
            disp_('maximprove= ',num2str_(maximprove),', jmax= ',int2str_(jmax))
        QZ,RZ,Y,xbase,scale=bcdfo_replace_in_Y_(QZ,RZ,ymax,Y,jmax,xbase,whichmodel,scale,shift_Y,Delta,normgx,kappa_ill,nargout=5)
        d[jmax]=norm_(ymax - Y[:,[0]])
        if (length_(find_(replaced == jmax)) == 0):
            replaced  = np.append( replaced, jmax )
    Y_radius=max_(d)
    if (verbose):
        print 'replaced='
        print replaced
        poised,Y_radius=bcdfo_poisedness_Y_(QZ,RZ,Y,eps_L,xbase,lSolver,whichmodel,hardcons,xl,xu,indfree,stratLam,scale,shift_Y,nargout=2)
        disp_('after everything: poisedness(Y) = ',num2str_(poised),' Y_radius  = ',num2str_(Y_radius))
        disp_('--- exit 2 bcdfo_repair_Y after round k=',num2str_(k),'  ---')
    return QZ,RZ,Y,replaced,maximprove,Y_radius,xbase,scale    