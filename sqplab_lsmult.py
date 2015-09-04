# -*- coding: utf-8 -*-
from __future__ import division
#try:
from runtime import *
from ecdfo_check_convex import *
from ecdfo_check_cond import *
from blls import *
from copy import copy
from ecdfo_global_variables import get_check_condition
from numpy import array
#except ImportError:
#    from smop.runtime import *

def sqplab_lsmult_(x=None,lb=None,ub=None,info_=None,options=None,values=None,*args,**kwargs):
    """
# [lm,info] = sqplab_lsmult (x,lb,ub,info,options,values);
#
# This procedure computes an exact 
# least-squares multiplier 'lm'. It solves in lm the
# quadratic optimization problem:
#
#   min || g+A'*lm ||^2
#   subject to possible bounds on lm,
#
# where g = info.g, lm = lm(1:n+mi+me) and A' = [ones(n) info.ai'
# info.ae']. 
#
# A multiplier associated with an inequality constraint having 
# - infinite lower and upper bounds vanishes,
# - infinite lower bound and finite upper bound is nonnegative,
# - finite lower bound and infinite upper bound is nonpositive,
# - finite lower bound and finite upper bound can have any sign.
# A lower (resp. upper) bound is considered as infinite if its value is
# empty or <= -options.inf (resp. empty or >= options.inf).
#
# On entry:
#   info.ci = ci(x)
#   lb: (optional or (n+mi) x 1) lower bounds on the n variables and mi
#       inequality constraints, is considered as empty if not present
#   ub: (optional or (n+mi) x 1) upper bounds on the n variables and mi
#       inequality constraints, is considered as empty if not present
#
# On return:
#   lm: computed least-squares multiplier, more precisely
#       lm(1:n): multiplier associated with the bounds on the variables
#       lm(n+1:n+mi): multiplier associated with the mi inequality
#           constraints
#       lm(n_mi+1:n+mi+me): multiplier associated with the me equality
#           constraints

#-----------------------------------------------------------------------
#
# Author: Jean Charles Gilbert, INRIA, and Anke Troeltzsch, DLR.
#
# Copyright 2008, 2009, INRIA.
#
# SQPlab is distributed under the terms of the Q Public License version
# 1.0.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the Q Public
# License version 1.0 for more details.
#
# You should have received a copy of the Q Public License version 1.0
# along with this program.  If not, see
# <http://doc.trolltech.com/3.0/license.html>.
#
#-----------------------------------------------------------------------
    """
#    varargin = cellarray(args)
    info=copy(info_)

    nargin = 6-[x,lb,ub,info,options,values].count(None)+len(args)
    
    #If somewhere we have done set_check_condition (in the tests for instance we have set_check_condition(0)), then
    #we get this value, otherwise we take '1' by default.
    try:
        check_condition=get_check_condition()
    except:
        check_condition=1    
# Initialization

    lm=array([])
    info.flag=values.success
    badcond=0
# Dimensions

    n=length_(info.g)
    me=0
    if (nargin >= 3):
        me=size_(info.ae,1)
# Check input arguments

    if (nargin < 4) or isempty_(lb):
        lb=- options.inf * ones_(n,1)
    else:
        lb=lb[:]
        if any_(size_(lb) != [n,1]):
            fprintf_('\n### sqplab_lsmult: incorrect size of lb\n\n')
            info.flag=values.fail_strange
            return lm,info
    if (nargin < 5) or isempty_(ub):
        ub=options.inf * ones_(n,1)
    else:
        ub=ub[:]
        if any_(size_(ub) != [n,1]):
            fprintf_('\n### sqplab_lsmult: incorrect size of ub\n\n')
            info.flag=values.fail_strange
            return lm,info
# Form the matrix A

    A=concatenate_([eye_(n),info.ae])
# Compute the lower (lo) and upper (up) bounds on lm

    lo=- inf * ones_(n + me,1)
    up=inf * ones_(n + me,1)
    for i in arange(0,n):
        if (lb[i] <= - options.inf):
            lo[i]=0
        if (ub[i] >= options.inf):
            up[i]=0
    # at this point, if there are lower and upper bounds on [x,ci], the multiplier can take any sign;
    # take additional constraint in case [x,ci] is active at a bound

        if (lb[i] > - options.inf) and (abs(x[i] - lb[i]) < options.dxmin):
            up[i]=0
            #disp_(['up ',str(i),' = 0']);
            # the multiplier must be <= 0
        if (ub[i] < options.inf) and (abs(x[i] - ub[i]) < options.dxmin):
            lo[i]=0  
            #disp_(['lo ',str(i),' = 0']);
            # the multiplier must be >= 0

    AA=A.dot(A.T)
# check condition of matrix AA

    if check_condition:
        cthreshold=1e+17
        AA,badcond=ecdfo_check_cond_(AA,cthreshold,options,nargout=2)
    check_convex=1
#  check that matrix AA is convex (otherwise convexify)

    if check_convex:
        AA=ecdfo_check_convex_(AA,options)
    Ag=A.dot(info.g)
    AAn=copy(AA)
    Agn=copy(Ag)
    lon=copy(lo)
    upn=copy(up)
    ifree=ones_(size_(lo))
    k=0
    for i in arange(0,length_(lo)):
        if lo[i] == up[i]:
#            AAn[k,:]=[]
#            AAn[:,k]=[]
#            Agn[k]=[]
#            lon[k]=[]
#            upn[k]=[]
            AAn=np.delete(AAn, k, 0)
            AAn=np.delete(AAn, k, 1)
            Agn=np.delete(Agn, k, 0)
            lon=np.delete(lon, k, 0)
            upn=np.delete(upn, k, 0)
            ifree[i]=0
        else:
            k=k + 1
    if not isempty_(ifree[ifree > 0]):
        sn,rn,op,exitc=blls_(AAn,- Agn,lon,upn,nargout=4)
        I=eye_(length_(lo))
        lm=np.delete(I, find_(ifree<=0), 1).dot(sn)
    else:
        lm=zeros_(size_(lo))
    return lm,info