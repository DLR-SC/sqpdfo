# -*- coding: utf-8 -*-
from __future__ import division
#try:
from runtime import *
#except ImportError:
#    from smop.runtime import *
from ecdfo_global_variables import get_prob
from numpy import array, zeros

def ecdfo_func_(x=None,*args,**kwargs):
    """
    #-----------------------------------------------------------------------
    # Computation of f, ci, ce
    #-----------------------------------------------------------------------
    """
#    varargin = cellarray(args)
#    nargin = 1-[x].count(None)+len(args)
    # Initialization
    prob=get_prob()
    ci=array([])
    if prob == 1:
        f=- (5 - (x[0] - 2) ** 2 - 2 * (x[1] - 1) ** 2)
        ce=zeros(1)
        ce=x[0] + 4 * x[1] - 3
        ce=ce.reshape(-1,1)
    elif prob == 2:
        f=2 * x[0] ** 2 + x[1] ** 2
        ce=zeros(1)
        ce=x[0] + x[1] - 1
        ce=ce.reshape(-1,1)
    elif prob == 3:
        f=x[0] ** 2 + x[1] ** 2 + x[2] ** 2
        ce=zeros(2)
        ce[0]=x[0] + x[1] + x[2]
        ce[1]=x[0] + 2 * x[1] + 3 * x[2] - 1
        ce=ce.reshape(-1,1)
    elif prob == 4:
        f=x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3]
        ce=zeros(3)
        ce[0]=x[0] + x[1] + x[2]
        ce[1]=x[0] + 2 * x[1] + 3 * x[2] - 1
        ce[2]=x[3] ** 3 - 1
        ce=ce.reshape(-1,1)
    elif prob == 5:
        #  Powells function from solnp - manual
        #  x* = (-1.717, 1.5957, 1.8272, -0.7636, -0.7636)
   
        f=exp_(x[0] * x[1] * x[2] * x[3] * x[4])
        ce=zeros(3)
        ce[0]=x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2 + x[4] ** 2 - 10
        ce[1]=x[1] * x[2] - 5 * x[3] * x[4]
        ce[2]=x[0] ** 3 + x[1] ** 3 + 1
        ce=ce.reshape(-1,1)
    elif prob ==6:
        ce=array([])
        f=-(0.592*((exp_(1)-1)*x[0])/((-0.408*x[0]+1)*(exp_(x[0])-1)) -1)
    elif prob==7:  #alkyl problem found here :http://www.gamsworld.org/global/globallib/alkyl.htm
        f=x[0]
        ce=zeros(8)
        ce[0]=6.3*x[4]*x[7]+x[0]-5.04*x[1]-0.35*x[2]-x[3]-3.36*x[5]
        ce[1]=-0.819672131147541*x[1]+x[4]-0.819672131147541*x[5]
        ce[2]=0.98*x[3]-x[6]*(0.01*x[4]*x[9]+x[3])
        ce[3]=x[1]*x[8]+10*x[2]+x[5]
        ce[4]=x[4]*x[11]-x[1]*(1.12+0.13167*x[8]-0.067*x[8]*x[8]) 
        ce[5]=x[7]*x[12]-0.01*(1.098*x[8]- 0.038*x[8]*x[8]) -0.325*x[6]  - 0.57425
        ce[6]=x[9]*x[13]+22.2*x[10]-35.82
        ce[7]=x[10]*x[14]-3*x[7]+1.33
        ce=ce.reshape(-1,1)
    msg=0
    return msg,f,ci,ce
   