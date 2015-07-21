# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 09:14:30 2015

@author: lien_ol
"""

from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from bcdfo_evalZ import *

def bcdfo_computeLj_(QZ=None,RZ=None,j=None,Y=None,whichmodel=None,scale=None,shift_Y=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 7-[QZ,RZ,j,Y,whichmodel,scale,shift_Y].count(None)+len(args)

    n,p1=size_(Y,nargout=2)
    q=((n + 1) * (n + 2)) / 2
    if (whichmodel == 3 and p1 < q):
        whichmodel=2
    if (whichmodel == 0):
        Lj=(QZ * (numpy.linalg.solve(RZ.T,[[zeros_(j - 1,1)],[1],[zeros_(p1 - j,1)]]))).T
    else:
        if (whichmodel == 1):
            if (p1 == n + 1 or p1 == q):
                Lj=(numpy.linalg.solve(RZ,QZ.T) * [[zeros_(j - 1,1)],[1],[zeros_(p1 - j,1)]]).T
                if (p1 == n + 1):
                    Lj[n + 2:q]=0
            else:
                if (shift_Y):
                    xbase=Y[:,1]
                    scaleY=0
                    for i in arange_(1,p1).reshape(-1):
                        Y[:,i]=Y[:,i] - xbase
                        scaleY=max_(scaleY,norm_(Y[:,i]))
                    Y=Y / scaleY
                rhs=matlabarray([[zeros_(j - 1,1)],[1],[zeros_(p1 + n + 1 - j,1)]])
                mualpha=(numpy.linalg.solve(RZ,(QZ.T * rhs))).T
                Lj[1:n + 1]=mualpha[p1 + 1:p1 + n + 1].T
                M=bcdfo_evalZ_(Y,q).T
                Lj[n + 2:q]=M[:,n + 2:q].T * mualpha[1:p1].T
        else:
            if (whichmodel == 2):
                if (p1 < q):
                    Lj=(QZ * (pinv_(RZ.T) * [[zeros_(j - 1,1)],[1],[zeros_(p1 - j,1)]])).T
                else:
                    Lj=(QZ * (numpy.linalg.solve(RZ.T,[[zeros_(j - 1,1)],[1],[zeros_(p1 - j,1)]]))).T
            else:
                if (whichmodel == 3):
                    Lj=(pinv_(RZ) * QZ.T * [[zeros_(j - 1,1)],[1],[zeros_(p1 - j,1)]]).T
    return Lj