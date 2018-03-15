# -*- coding: utf-8 -*-

from runtime import *
import ecdfo_global_variables as glob
from numpy import array, zeros
import time

def ecdfo_func_(x=None,*args,**kwargs):
    """
    #-----------------------------------------------------------------------
    # Computation of f, ci, ce
    #-----------------------------------------------------------------------
    """
#    varargin = cellarray(args)
#    nargin = 1-[x].count(None)+len(args)
    # Initialization
    prob=glob.get_prob()
    ci=array([])
    ce=array([])
    msgf = 0
    msgce = 0
    msg = 0
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
        f=-(0.592*((exp_(1)-1)*x[0])/((-0.408*x[0]+1)*(exp_(x[0])-1)) -1)
    elif prob==7:
        # alkyl problem found here :http://www.gamsworld.org/global/globallib/alkyl.htm
        # best known solution found by Baron f*=-1.76499964633
        # x1=-1.76499964633; x2=1.70370291772; x3=1.58470999786; x4=0.543084200389; x5=3.03582206371; x6=2.0; x6=-0.882499871995;
        # x7=0.901319381076; x8=0.95; x8=-17.3201583536; x9=10.4754765591; x10=1.56163636364; x11=1.53535353535;
        # x12=0.99; x12=2.06361745003; x13=0.99; x13=21.9374937958; x14= 1.11111111111; x14=-0.488775780565;
        # x15=0.99; x15=10.7759591879;
        # e1=1.0; e2=-5.13314985138; e3=12.7096102404; e4=0.035; e5=-0.679755732296; e6=-23.0920987324; e7=0.312989497393; e8=-7.01855236578
        f=x[0]
        ce=zeros(8)
        ce[0]=6.3*x[4]*x[7]+x[0]-5.04*x[1]-0.35*x[2]-x[3]-3.36*x[5]
        ce[1]=-0.819672131147541*x[1]+x[4]-0.819672131147541*x[5]
        ce[2]=0.98*x[3]-x[6]*(0.01*x[4]*x[9]+x[3])
        ce[3]=-x[1]*x[8]+10*x[2]+x[5]
        ce[4]=x[4]*x[11]-x[1]*(1.12+0.13167*x[8]-0.0067*x[8]*x[8]) 
        ce[5]=x[7]*x[12]-0.01*(1.098*x[8]- 0.038*x[8]*x[8]) -0.325*x[6]  - 0.57425
        ce[6]=x[9]*x[13]+22.2*x[10]-35.82
        ce[7]=x[10]*x[14]-3*x[7]+1.33
        ce=ce.reshape(-1,1)
    elif prob==100:
        ffunc = glob.get_filename_f()
        cfunc = glob.get_filename_ce()
        # GTlab problem
        #scalefacX = array([[1000., 0.1, 10., 1., 1.]]).T
        #x = copy(x) / scalefacX
        xc = copy(x.T)
        xc = list(xc[0])

        try:
            f = ffunc(xc)
        except:
            f = np.inf
            msgf = 1
        if cfunc != '':
            try:
                ce = cfunc(xc)
                ce = array(ce).reshape(-1, 1)
            except:
                msgce = 1
        if msgf > 0 or msgce > 0:
            msg = 'Error: error during calculating f (and/or ce) in user-defined problem !'
    elif prob==1000:
        # CUTEr problems
        cproblem=glob.get_prob_cuter()
        (f, c)=cproblem.objcons(x.reshape(-1))
        ce=c.reshape(-1,1)
    else:
        msg = 'Error: Problem number prob='+str(prob)+' is not defined in ecdfo_func.py ! \nUserdefined problems have problem number prob=100.'
        f = np.inf

    return msg,f,ci,ce
   
