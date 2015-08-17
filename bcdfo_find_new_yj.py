# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 13:23:34 2015

@author: lien_ol
"""

from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from bcdfo_computeLj import bcdfo_computeLj_
from bcdfo_gradP import bcdfo_gradP_
from bcdfo_hessP import bcdfo_hessP_
from bcdfo_solve_TR_MS import bcdfo_solve_TR_MS_
from numpy import *



def bcdfo_find_new_yj_(QZ=None,RZ=None,Y=None,j=None,Delta=None,eps_L=None,xbase=None,lSolver=None,whichmodel=None,scale=None,shift_Y=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 11-[QZ,RZ,Y,j,Delta,eps_L,xbase,lSolver,whichmodel,scale,shift_Y].count(None)+len(args)

    verbose=0
    n=size_(Y,1)
    ynew=zeros_(1,n)
    improvement=0
    if (verbose):
        disp_('--------- enter find_new_yj ')
    if (j < 2):
        return ynew,improvement,msgTR
    Lj=bcdfo_computeLj_(QZ,RZ,j,Y,whichmodel,scale,shift_Y)
    if (length_(find_(isnan(Lj))) != 0 or length_(find_(~ isreal(Lj))) != 0 or length_(find_(isinf(Lj))) != 0):
        msgTR='Error0: Lagrange polynomial contains NaN or Inf or nonreal components!!'
        if (verbose):
            disp_(msgTR)
        return ynew,improvement,msgTR
    if (lSolver == 2):
        Delta=sqrt_(n) * Delta
    if (shift_Y):
        g=bcdfo_gradP_(Lj,zeros_(n,1),xbase,scale,0)
        H=bcdfo_hessP_(Lj,zeros_(n,1),xbase,scale,0)
        pstep,_lambda,norms,pvalue,gplus,nfact,neigd,msgTR,hardcase=bcdfo_solve_TR_MS_(g,H,Delta * scale[1],eps_L,nargout=9)
        pstep=pstep / scale[1]
        mstep,_lambda,norms,mvalue,gplus,nfact,neigd,msgTR,hardcase=bcdfo_solve_TR_MS_(- g,- H,Delta * scale[1],eps_L,nargout=9)
        mstep=mstep / scale[1]
    else:
        g=bcdfo_gradP_(Lj,Y[:,0],xbase,scale,0)
        H=bcdfo_hessP_(Lj,Y[:,0],xbase,scale,0)
        pstep,_lambda,norms,pvalue,gplus,nfact,neigd,msgTR,hardcase=bcdfo_solve_TR_MS_(g,H,Delta,eps_L,nargout=9)
        mstep,_lambda,norms,mvalue,gplus,nfact,neigd,msgTR,hardcase=bcdfo_solve_TR_MS_(- g,- H,Delta,eps_L,nargout=9)
    if (verbose):
        disp_(' === find_new_yj: j = ',int2str_(j),' positive value = ',num2str_(pvalue),' step:')
        pstep.T
        disp_(' === find_new_yj: j = ',int2str_(j),' negative value = ',num2str_(mvalue),' step:')
        mstep.T
    if (mvalue < pvalue):
        improvement=abs(mvalue)
        ynew=Y[:,0].reshape(-1,1) + mstep
    else:
        improvement=abs(pvalue)
        ynew=Y[:,0].reshape(-1,1) + pstep
    if (verbose):
        disp_('--------- exit find_new_yj ')
    return ynew,improvement,msgTR