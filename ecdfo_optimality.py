# -*- coding: utf-8 -*-

from runtime import *
import ecdfo_global_variables as glob
from copy import copy
from numpy import arange, logical_and, zeros, diag, ones, concatenate, array


def ecdfo_optimality_(x=None,lm=None,lb=None,ub=None,info_=None,options=None,*args,**kwargs):

    ##############################################################
    #  check for optimality
    ##############################################################

    info=copy(info_)
    nbr_slacks = glob.get_nbr_slacks()
    sl = glob.get_slacks()

    # Dimensions

    n=length_(info.g)
    mi=length_(info.ci)
    me=length_(info.ce)

    # Gradient of the Lagrangian

    info.glag=copy(info.g)

    # consider only active and violated bounds for gradient computation
    bounds=(abs(lb - x) < 1e-05) + (abs(ub - x) < 1e-05)

    # contribution of the variable bounds

    if options.verbose >= 4:
        fprintf_('\n     lb             x            ub             glag            lm\n')
        for i in range(0,n):
            fprintf_('%12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n'%(lb[i],x[i],ub[i],\
            info.glag[i],lm[i]))
            
    I=find_(bounds[0:n])
    info.glag[I]=info.glag[I] + lm[I]

    # contribution of the equality constraints
    
    if me > 0:
        if nbr_slacks > 0:
            hh1 = concatenate((zeros((me-nbr_slacks,nbr_slacks)),-2*diag(sl.T[0])))
            gconstraints = concatenate((info.ae,hh1),axis=1)
            glag = concatenate((info.glag,zeros((nbr_slacks,1))))
            glag = glag + gconstraints.T.dot(lm[n:n+me])
            info.glag=array([glag.T[0][0:n]]).T
        else:
            info.glag=info.glag + info.ae.T.dot(lm[n + mi:n + mi + me])
           
    # Feasibility (a vector)

    if nbr_slacks > 0:
        feas = info.ce - concatenate((zeros((me-nbr_slacks,1)),sl**2))
    else:
        feas = info.ce

    # Complementarity (a vector)

    #v=concatenate_([x,info.ci])
    compl=zeros_(n + mi,1)

    # lower bounds on v
    
    #v(I) is inactive at its lower bounds
    #I=find_(logical_and(lb > - options.inf,abs(lb - v) > options.dxmin)) 
    #if not isempty_(I):  
    #    compl[I]=max_(compl[I],max_(0,- lm[I])) # lm(I) must be >= 0

    # upper bounds on v

    # v(I) is inactive at its upper bounds
    #I=find_(logical_and(ub < options.inf,abs(ub - v) > options.dxmin)) 
    #if not isempty_(I):  
    #    compl[I]=max_(compl[I],max_(0,lm[I])) # lm(I) must be <= 0

    return feas,compl,info
