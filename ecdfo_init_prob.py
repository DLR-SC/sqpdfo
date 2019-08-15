# -*- coding: utf-8 -*-

from runtime import *
from numpy import inf, array
from ecdfo_global_variables import *

def ecdfo_init_prob_(prob=None,*args,**kwargs):
    """
    # This function returns the dimensions of the problem:
# . n  = number of variables,
# . nb = number of variables with bounds,
# . mi = number of inequality constraints,
# . me = number of equality constraints.
    """

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
    elif prob == 10:  # problem 19 from Hock-Schittkowskis collection
        n=2
        nb=4
        me=0
        mi=2
        x0 = array([[20.1,5.84]])
        lx = array([[13.0,0.0]]).reshape(-1,1)
        ux = array([[100.0,100.0]]).reshape(-1,1)
    elif prob == 11:  # problem 21 from Hock-Schittkowskis collection
        n=2
        nb=4
        me=0
        mi=1
        x0 = array([[-1.0,-1.0]])
        lx = array([[2.0,-50.0]]).reshape(-1,1)
        ux = array([[50.0,50.0]]).reshape(-1,1)
    elif prob == 12:  # problem 35 (Beale) from HS collection
        n=3
        nb=3
        me=0
        mi=1
        x0 = array([[0.5,0.5,0.5]])
        lx = array([[0.0,0.0,0.0]]).reshape(-1,1)
        ux = array([[1e20,1e20,1e20]]).reshape(-1,1)
    elif prob == 13:  # problem 76 from Hock-Schittkowskis collection
        n=4
        nb=4
        me=0
        mi=3
        x0 = array([[0.5,0.5,0.5,0.5]])
        lx = array([[0.0,0.0,0.0,0.0]]).reshape(-1,1)
        ux = array([[1e20,1e20,1e20,1e20]]).reshape(-1,1)
    elif prob == 14:  # problem 44 from Hock-Schittkowskis collection
        n=4
        nb=4
        me=0
        mi=6
        x0 = array([[0.0,0.0,0.0,0.0]])
        lx = array([[0.0,0.0,0.0,0.0]]).reshape(-1,1)
        ux = array([[1e20,1e20,1e20,1e20]]).reshape(-1,1)
    elif prob == 15:
        n=2
        nb=4
        me=2
        mi=1
        x0 = array([[-1.2,1.0]])
        lx = array([[-5.0,-5.0]]).reshape(-1,1)
        ux = array([[10.0,10.0]]).reshape(-1,1)
    elif prob == 16:
        n=2
        nb=4
        me=1
        mi=2
        x0 = array([[0.3,0.3]])
        lx = array([[-5.0,-5.0]]).reshape(-1,1)
        ux = array([[10.0,10.0]]).reshape(-1,1)
    elif prob==1000:
        #Warning : here the CUTEst interface from this website has to be 
        #installed in order to use CUTEst problems :
        #https://jfowkes.github.io/pycutest/_build/html/index.html 
        #Thanks to Jaroslav Fowkes and Lindon Roberts
        
        cproblem=get_prob_cuter()
        n=cproblem.n
        m=cproblem.m
        me= sum(cproblem.is_eq_cons)
        mi=m-me
        #x0=info['x'].reshape(-1,1)
        #lx=info['bl'].reshape(-1,1)
        #ux=info['bu'].reshape(-1,1)
        x0=cproblem.x0.reshape(-1,1)
        lx=cproblem.bl.reshape(-1,1)
        ux=cproblem.bu.reshape(-1,1)
        nb=sum_(min_((lx[0:n]>-inf)+(inf > ux[0:n]),1))
    else:
        sys.exit('Problem number is not set!\n'\
            'Please import ecdfo_global_variables and use set_prob(nbr)\n'\
            'where nbr is 1,...,5 for test examples or 1000 for a CUTEst problem.')
    info=0
    return x0,lx,ux,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info
