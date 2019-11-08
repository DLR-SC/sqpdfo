# -*- coding: utf-8 -*-

from runtime import *
from copy import copy
from sqpdfo_swap_in_Y import sqpdfo_swap_in_Y_
from numpy import array, zeros, concatenate
import sqpdfo_global_variables as glob

def sqpdfo_find_smallf_(c_=None,QZ_=None,RZ_=None,Y_=None,fY_=None,ciY_=None,ceY_=None,\
    ind_Y_=None,i_xbest_=None,cur_degree_=None,indfree_=None,x_=None,\
    xl_=None,xu_=None,fx_=None,dstatus_=None,whichmodel_=None,scale_=None,\
    shift_Y_=None,Delta_=None,normgx_=None,kappa_ill_=None,sigma_=None,\
    info_=None,*args,**kwargs):
    
#############################################################################
#  Subroutine finds the smallest value in fY, which are the associated function
#  values of the set Y. There are only points inside the bounds considered
#  and those which are and no dummy points!
#  Exchanges the best point in the current interpolation set Y if a smaller
#  function value was found.
#
#  INPUT:
#
#  c           : contains a bunch of constants
#  QZ, RZ      : the QR factors of the (possibly shifted) matrix containing
#                the polynomial expansion of the interpolation points
#  Y           : interpolation set
#  fY          : function values associated to the current interpolation points
#  ind_Y       : indices of the points in Y out of X, the set of all points
#  i_xbest     : index of the best point
#  cur_degree  : number of interpolation points
#  indfree     : number of free variables
#  x           : best point
#  xl, xu      : lower/upper bounds on the variables
#  fx          : best function value
#  dstatus     : status vector of dummy points in X
#  whichmodel  : kind of model to build
#  scale       : model diagonal scaling
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#  Delta       : trust-region radius
#  normgx      : infinity norm of the projected gradient
#  kappa_ill   : threshold to declare a system matrix as ill-conditioned
#
#  OUTPUT:
#
#  (possibly) updated INPUT values
#
#  PROGRAMMING: A. Troeltzsch, March 2013.
#
#  DEPENDENCIES: sqpdfo_swap_in_Y
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#############################################################################
    
    c=copy(c_)
    QZ=copy(QZ_)
    RZ=copy(RZ_)
    Y=copy(Y_)
    fY=copy(fY_)
    ciY=copy(ciY_)
    ceY=copy(ceY_)
    ind_Y=copy(ind_Y_)
    i_xbest=copy(i_xbest_)
    cur_degree=copy(cur_degree_)
    indfree=copy(indfree_)
    x=copy(x_)
    xl=copy(xl_)
    xu=copy(xu_)
    fx=copy(fx_)
    dstatus=copy(dstatus_)
    whichmodel=copy(whichmodel_)
    scale=copy(scale_)
    shift_Y=copy(shift_Y_)
    Delta=copy(Delta_)
    normgx=copy(normgx_)
    kappa_ill=copy(kappa_ill_)
    sigma=copy(sigma_)
    info=copy(info_)
    
    nbr_slacks = glob.get_nbr_slacks()
    sl = glob.get_slacks()    

    norm_ceY = zeros_(1,cur_degree).reshape(-1)
    
    #  Select indices of all points inside the bounds and which are not a dummy point

    dummy_set=find_(dstatus == c.dummy)
    ind_insideBounds=array([])
    
    for i in range(0,cur_degree):
        if((isempty_(find_(logical_or_(Y[:,i] < xl[indfree] , Y[:,i] > xu[indfree]),1)))\
            and (isempty_(find_(dummy_set == ind_Y[i],1)))):
            ind_insideBounds=concatenate_([ind_insideBounds,[i]],axis=1)
        else:
            ind_insideBounds=concatenate_([ind_insideBounds,[0]],axis=1)
            
    #  Compute merit function values for all points

    if ceY.any():
        if nbr_slacks > 0:
            me = size_(ceY,1)
            for i in range(0,cur_degree):
                norm_ceY[i]=norm_(ceY[:,i] - \
                           concatenate((zeros((me-nbr_slacks,1)),sl**2)).T[0])   
        else:
            for i in range(0,cur_degree):
                norm_ceY[i]=norm_(ceY[:,i])

    meritY=fY + int(sigma) * norm_ceY
    
    #  Find the smallest function value among them

    fmin,imin=min_(meritY[ind_insideBounds],nargout=2)
    
    #  Exchange best point only if the new smallest fmin is strictly smaller than fx
    #  (to avoid a change of the center of the interpolation at equal function values)

    if (imin != 1 and fmin < meritY[0]):
    
        QZ,RZ,Y,ind_Y,fY,ciY,ceY,x,scale=\
        sqpdfo_swap_in_Y_(0,imin,QZ,RZ,Y,ind_Y,fY,ciY,ceY,x,whichmodel,\
            scale,shift_Y,Delta,normgx,kappa_ill,nargout=9)
            
        fx=copy(fY[0])
        i_xbest=copy(ind_Y[0])
        
        if (not shift_Y):
            x=copy(Y[:,[0]])
            
    info.f=copy(fY[0])
    
    if length_(ceY) > 0:
        info.ce=copy(ceY[:,[0]])
        
    if length_(ciY) > 0:
        info.ci=copy(ciY[:,[0]])
        
    return x,fx,QZ,RZ,Y,fY,ciY,ceY,ind_Y,i_xbest,scale,info
