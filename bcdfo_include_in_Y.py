# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 15:32:33 2015

@author: lien_ol
"""
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from bcdfo_evalL import bcdfo_evalL_
from bcdfo_replace_in_Y import bcdfo_replace_in_Y_

def bcdfo_include_in_Y_(x=None,QZ_=None,RZ_=None,Y_=None,choice_set=None,poisedness_threshold=None,criterion=None,xbase_=None,whichmodel=None,succ=None,scale_=None,shift_Y=None,Delta=None,normgx=None,kappa_ill=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 15-[x,QZ,RZ,Y,choice_set,poisedness_threshold,criterion,xbase,whichmodel,succ,scale,shift_Y,Delta,normgx,kappa_ill].count(None)+len(args)
    
    Y=copy_(Y_)
    QZ=copy_(QZ_)
    RZ=copy_(RZ_)
    xbase=copy_(xbase_)
    scale=copy_(scale_)

    if (length_(choice_set) == 0):
        pos=0
        return QZ,RZ,Y,pos,xbase,scale
    Lvals=bcdfo_evalL_(QZ,RZ,Y,choice_set,x,xbase,whichmodel,scale,shift_Y)
    choice=find_(abs_(Lvals) > poisedness_threshold)
    lc=length_(choice)
    if (lc == 0):
        pos=0
        return QZ,RZ,Y,pos,xbase,scale
    crit_val=0
    pos=0
    for i in arange_(1,lc).reshape(-1):
        j=choice[i]
        if (criterion == char('weighted')):
            if (succ == 1):
                cv=norm_(Y[:,j] - x) ** 2 * abs_(Lvals[j])
            else:
                cv=norm_(Y[:,j] - Y[:,1]) ** 2 * abs_(Lvals[j])
        else:
            if (criterion == char('standard')):
                cv=abs_(Lvals[j])
            else:
                if (criterion == char('distance')):
                    if (succ == 1):
                        cv=norm_(Y[:,j] - x)
                    else:
                        cv=norm_(Y[:,j] - Y[:,1])
        if (cv > crit_val):
            pos=copy_(j)
            crit_val=copy_(cv)
    if (pos == 0):
        return QZ,RZ,Y,pos,xbase,scale
    QZ,RZ,Y,xbase,scale=bcdfo_replace_in_Y_(QZ,RZ,x,Y,pos,xbase,whichmodel,scale,shift_Y,Delta,normgx,kappa_ill,nargout=5)
    return QZ,RZ,Y,pos,xbase,scale