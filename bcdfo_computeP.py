# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 09:14:31 2015

@author: lien_ol
"""

from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from bcdfo_evalZ import *

def bcdfo_computeP_(QZ=None,RZ=None,Y=None,fY=None,whichmodel=None,P_old=None,ind_Y=None,i_xold=None,i_xplus=None,g=None,scale=None,shift_Y=None,Delta0=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 13-[QZ,RZ,Y,fY,whichmodel,P_old,ind_Y,i_xold,i_xplus,g,scale,shift_Y,Delta0].count(None)+len(args)

    n,p1=size_(Y,nargout=2)
    badcond=0
    q=((n + 1) * (n + 2)) / 2
    if (whichmodel == 3 and p1 < q):
        whichmodel=2
    if (whichmodel == 0):
        P=(QZ * (numpy.linalg.solve(RZ.T,fY.T))).T
    else:
        if (whichmodel == 1):
            n_rhs=size_(fY,1)
            if (p1 == n + 1 or p1 == q):
                P[1:n_rhs,1:p1]=(numpy.linalg.solve(RZ,(QZ.T * fY.T))).T
                if (p1 == n + 1):
                    P[1:n_rhs,n + 2:q]=0.0
            else:
                if (shift_Y):
                    xbase=Y[:,1]
                    scaleY=0
                    for i in arange_(1,p1).reshape(-1):
                        Y[:,i]=Y[:,i] - xbase
                        scaleY=max_(scaleY,norm_(Y[:,i]))
                    Y=Y / scaleY
                P=matlabarray([])
                for i in arange_(1,n_rhs).reshape(-1):
                    rhs=concatenate_([fY[i,:],zeros_(1,n + 1)], axis=1)
#                    rhs=matlabarray([fY[i,:],zeros_(1,n + 1)])
                    mualpha=(numpy.linalg.solve(RZ,(QZ.T * rhs.T))).T
                    P_i[1:n + 1]=mualpha[p1 + 1:p1 + n + 1].T
                    M=bcdfo_evalZ_(Y,q).T
                    P_i[n + 2:q]=M[:,n + 2:q].T * mualpha[1:p1].T
                    P=concatenante_([P, P_i])
#                    P=matlabarray([[P],[P_i]])
        else:
            if (whichmodel == 2):
                P=(QZ * (numpy.linalg.solve(RZ.T,fY.T))).T
            else:
                if (whichmodel == 3):
                    P=(pinv_(RZ) * QZ.T * fY.T).T
    return P