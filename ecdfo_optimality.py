# -*- coding: utf-8 -*-
from __future__ import division
#try:
from runtime import *
from copy import copy
from numpy import arange
#except ImportError:
    #from smop.runtime import *
def ecdfo_optimality_(x=None,lm=None,lb=None,ub=None,info_=None,options=None,*args,**kwargs):
    """    
%-----------------------------------------------------------------------
%
% Author: Jean Charles Gilbert, INRIA
%     and Anke Troeltzsch, DLR.
%
% Copyright 2008, 2009, INRIA. 2013, DLR.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the Q Public
% License version 1.0 for more details.
%
% You should have received a copy of the Q Public License version 1.0
% along with this program.  If not, see
% <http://doc.trolltech.com/3.0/license.html>.
%
%-----------------------------------------------------------------------
    """
#    varargin = cellarray(args)
#    nargin = 6-[x,lm,lb,ub,info,options].count(None)+len(args)

    info=copy(info_)
    n=length_(info.g)
    mi=length_(info.ci)
    me=length_(info.ce)
    info.glag=copy(info.g)
    gradlag=info.g
    bounds=(abs(lb - x) < 1e-05) + (abs(ub - x) < 1e-05)
    if options.verbose >= 4:
        fprintf_('\n     lb             x            ub             glag            lm\n')
        for i in range(0,n):
            fprintf_('%12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n'%(lb[i],x[i],ub[i],info.glag[i],lm[i]))
    I=find_(bounds[0:n])
    info.glag[I]=info.glag[I] + lm[I]
    boundsmult=lm[I]
    gradlag=info.glag
    if mi > 0:
        I=find_(bounds[n:n + mi])
        info.glag=info.glag + info.ai(I,arange()).T.dot(lm[n + I])
    if me > 0:
        derivequcons=info.ae
        equmult=lm[n + mi:n + mi + me]
        info.glag=info.glag + info.ae.T.dot(lm[n + mi:n + mi + me])
        gradlag=info.glag
#    feas=array([[max_(0,max_([[x],[info.ci]] - ub,lb - [[x],[info.ci]]))],[info.ce]])
#    v=array([[x],[info.ci]])
    feas=concatenate_([max_(0,max_(concatenate_([x, info.ci]) - ub, lb - concatenate_([x,info.ci]))) ,
										info.ce])
    v=concatenate_([x,info.ci])
    compl=zeros_(n + mi,1)
    I=find_(logical_and(lb > - options.inf,abs(lb - v) > options.dxmin))
    if not isempty_(I):
        compl[I]=max_(compl[I],max_(0,- lm[I]))
    I=find_(logical_and(ub < options.inf,abs(ub - v) > options.dxmin))
    if not isempty_(I):
        compl[I]=max_(compl[I],max_(0,lm[I]))
    return feas,compl,info
