# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 13:04:32 2015

@author: lien_ol
"""
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
    
from bcdfo_evalP import bcdfo_evalP_
import numpy
from copy import copy
    
def bcdfo_evalL_(QZ=None,RZ=None,Y_=None,choice_set=None,x=None,xbase_=None,whichmodel=None,scale=None,shift_Y=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 9-[QZ,RZ,Y,choice_set,x,xbase,whichmodel,scale,shift_Y].count(None)+len(args)

    #Copy only if necessary
    if (shift_Y):
        Y=copy(Y_)
        xbase=copy(xbase_)
    else:
        Y=Y_
        xbase=xbase_        
    
    choice_set=choice_set.reshape(-1)
    n,p1=size_(Y,nargout=2)
    I=eye_(p1)
    lc=length_(choice_set)
    q=((n + 1) * (n + 2)) / 2
    values=zeros_(p1,1)
    if (whichmodel == 3 and p1 < q):
        whichmodel=2
    if (whichmodel == 0):
        values[choice_set]=bcdfo_evalP_(I[choice_set,:].dot( (numpy.linalg.solve(RZ,QZ.T))),x,xbase,scale,shift_Y)
    else:
        if (whichmodel == 1):
            if (p1 == n + 1 or p1 == q):
                values[choice_set]=bcdfo_evalP_(I[choice_set,:] .dot( (QZ / RZ.T)),x,xbase,scale,shift_Y)
            else:
                if (shift_Y):
                    xbase=copy(Y[:,0])
                    scaleY=0
                    for i in range(0,p1):
                        Y[:,i]=Y[:,i] - xbase
                        scaleY=max_(scaleY,norm_(Y[:,i]))
                    Y=Y / scaleY
                M=bcdfo_evalZ_(Y,q).T
                MQ=M[:,n + 1:q]
                if (shift_Y):
                    phi=bcdfo_evalZ_((x - xbase) * scale[1],q)
                else:
                    phi=bcdfo_evalZ_(x,q)
                values[choice_set]=concatenate_([  I[choice_set,:],  zeros_(lc,n + 1).dot(QZ .dot(numpy.linalg.solve(RZ.T,concatenate_([MQ.dot(phi[n + 1:q]),phi[0:n + 1]]))))],axis=1)
        else:
            if (whichmodel == 2):
                if (p1 < q):
                    values[choice_set]=bcdfo_evalP_(I[choice_set,:].dot( (pinv_(RZ).dot(QZ.T))),x,xbase,scale,shift_Y)
                else:
                    values[choice_set]=bcdfo_evalP_(I[choice_set,:].dot( (numpy.linalg.solve(RZ,QZ.T))),x,xbase,scale,shift_Y)
            else:
                if (whichmodel == 3):
                    values[choice_set]=bcdfo_evalP_(I[choice_set,:].dot((QZ .dot(pinv_(RZ.T)))),x,xbase,scale,shift_Y)
    return values