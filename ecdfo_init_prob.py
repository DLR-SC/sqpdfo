# -*- coding: utf-8 -*-
from __future__ import division
#try:
from runtime import *
from numpy import inf, array
from ecdfo_global_variables import *
#except ImportError:
    #from smop.runtime import *
def ecdfo_init_prob_(prob=None,*args,**kwargs):
    """
    # This function returns the dimensions of the problem:
# . n  = number of variables,
# . nb = number of variables with bounds,
# . mi = number of inequality constraints,
# . me = number of equality constraints.
    """
#    varargin = cellarray(args)
#    nargin = 1-[prob].count(None)+len(args)

    eps = 2.220446049250313e-16
# Set output variables
   
    x0=array([])
    lx=array([])
    ux=array([])
    li=array([])
    ui=array([])
# fid of the output file

    set_fileoutput(1)
    set_simul_not_initialized(1)
    #dxmin = sqrt(eps);

    dxmin=1e-06
    dcimin=sqrt_(eps)
    infb=1e+20
    if prob == 1:
        n=2
        nb=2
        mi=0
        me=1
        x0=array([[4.6],[0.0]]).T
        lx=array([1.95,- 1e+20]).reshape(-1,1)
        ux=array([1e+20,0.3]).reshape(-1,1)
    elif prob == 2:
        n=2
        nb=0
        mi=0
        me=1
        x0=array([[- 1],[2.54378]]).T
        lx=- inf * ones_(n,1)
        ux=inf * ones_(n,1)
    elif prob == 3:
        nb=2
        mi=0
        me=2
        x0=array([[0.0],[0.0],[0.5]]).T
        n=length_(x0)
        lx=array([- 0.5,0.0,- inf]).reshape(-1,1)
        ux=array([inf,inf,inf]).reshape(-1,1)
    elif prob == 4:
        nb=0
        mi=0
        me=3
        x0=array([[1.0],[1.0],[1.0],[0.0]]).T
        n=length_(x0)
        lx=- inf * ones_(n,1)
        ux=inf * ones_(n,1)
    elif prob == 5:
        nb=0
        mi=0
        me=3
        x0=array([[- 2.0],[2.0],[2.0],[1.0],[1.0]]).T
        n=5
        lx=- inf * ones_(n,1)
        ux=inf * ones_(n,1)
    elif prob==6:
        n=1
        nb=1
        mi=0
        me=0
        x0=array([[0.6]])
        lx=array([0.5]).reshape(-1,1)
        ux=array([0.8]).reshape(-1,1)
    elif prob==7: #alkyl problem found here :http://www.gamsworld.org/global/globallib/alkyl.htm
        n=15
        nb=14
        me=8
        mi=0
        x0=array([[-0.9,1.745,1.2,1.1,3.048,1.974,0.893,0.928,8,3.6,1.50,1,1,1,1]]).T
        lx=array([-inf,0,0,0,0,0,0.85,0.9,3,1.2,1.45,0.99,0.99,0.9,0.99]).reshape(-1,1)
        ux=array([inf,2,1.6,1.2,5,2,0.93,0.95,12,4,1.62,1.01010101010101,1.01010101010101,1.11111111111111,1.01010101010101]).reshape(-1,1)

    else:
        #Warning : here the CUTEr interface from this website has to be installed in order to use CUTEr problems :
        #http://fides.fe.uni-lj.si/~arpadb/software-pycuter.html Thanks to Prof. Dr. Árpád Bűrmen
        set_prob_cuter(prob)
        
        cproblem=get_prob_cuter()
        info=cproblem.getinfo()
        n=info['n']
        mi=0
        me=info['m']
        x0=info['x'].reshape(-1,1)
        lx=info['bl'].reshape(-1,1)
        ux=info['bu'].reshape(-1,1)
        nb=sum_(min_((lx[0:n]>-inf)+(inf > ux[0:n]),1))
    info=0
    return x0,lx,ux,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info
