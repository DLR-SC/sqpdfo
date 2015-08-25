# -*- coding: utf-8 -*-
from __future__ import division
#try:
from runtime import *
from bcdfo_computeP import bcdfo_computeP_
from bcdfo_hessP import bcdfo_hessP_
from copy import copy
from numpy import zeros
#except ImportError:
    #from smop.runtime import *

def ecdfo_computeHessian_(simul=None,x=None,null_step=None,constrained_pbl=None,lm=None,M_=None,n=None,me=None,mi=None,s=None,gx=None,gci=None,gce=None,info_=None,options=None,values_=None,fcmodel=None,Y=None,fY=None,ciY=None,ceY=None,sigma=None,scale=None,shift_Y=None,QZ=None,RZ=None,whichmodel=None,ind_Y=None,i_xbest=None,m=None,*args,**kwargs):
    """
    %-----------------------------------------------------------------------
    % approximation of the Hessian of the Lagrangian
    %-----------------------------------------------------------------------
    """
#    varargin = cellarray(args)
#    nargin = 30-[simul,x,null_step,constrained_pbl,lm,M,n,me,mi,s,gx,gci,gce,info,options,values,fcmodel,Y,fY,ciY,ceY,sigma,scale,shift_Y,QZ,RZ,whichmodel,ind_Y,i_xbest,m].count(None)+len(args)


    M=copy(M_)
    info=copy(info_)
    values=copy(values_)

    norm_ceY=zeros(size_(Y,2))

    

    pc=1.0
    if options.algo_method == values.newton:
        info.nsimul[2]=info.nsimul[2] + 1
        outdic,tmp,tmp,tmp,tmp,info.g,info.ai,info.ae=simul[3,x]
        if outdic:
            info=sqplab_badsimul_(outdic,info,options,values)
            return M,pc,info
        info.nsimul[4]=info.nsimul[4] + 1
        outdic,M=simul[5,x,lm]
        if outdic:
            info=sqplab_badsimul_(outdic,info,options,values)
            return M,pc,info
    else:
        if not null_step:
            if options.hess_approx == values.bfgs:
                if options.verbose >= 4:
                    if constrained_pbl:
                        fprintf_(options.fout,'\nBFGS update:\n')
                    else:
                        fprintf_(options.fout,'\nBFGS inverse update:\n')
                y=- info.g
                if me:
                    y=y - info.ae.T.dot(lm[n + mi:n + mi + me])
                info.g=gx
                info.ai=gci
                info.ae=gce
                y=y + info.g
                if me:
                    y=y + info.ae.T.dot(lm[n + mi:n + mi + me])
                if options.bfgs_restart > 0 and mod_(info.nsimul[1],options.bfgs_restart) == 0:
                    M=eye_(size_(M))
                    pc=2.0
                else:
                    first=0
                    if info.niter == 1:
                        first=1
                    M,pc,info,values=sqplab_bfgs_(M,y,s,first,info,options,values,nargout=4)
                    if info.flag == values.fail_strange:
                        M=eye_(size_(M))
                        M,pc,info,values=sqplab_bfgs_(M,y,s,first,info,options,values,nargout=4)
                        if info.flag == values.fail_strange:
                            return M,pc,info
            else:
                if options.hess_approx == values.model:
                    info.g=gx
                    info.ai=gci
                    info.ae=gce
                    cur_degree=size_(Y,2)
                    if me + mi > 0:
                        if length_(ceY) > 0:
                            for i in range(0,cur_degree):
                                norm_ceY[i]=norm_(ceY[:,i])
                        else:
                            norm_ceY=zeros(cur_degree)
                        meritfY=fY + sigma *norm_ceY
                        model=bcdfo_computeP_(QZ,RZ,Y,meritfY.reshape(1,-1),whichmodel,fcmodel[[0],:],ind_Y,i_xbest,m,gx,scale,shift_Y)
                        M=bcdfo_hessP_(model,x,x,scale,shift_Y)
                    else:
                        M=bcdfo_hessP_(fcmodel[[0],:],x,x,scale,shift_Y)
                else:
                    if options.verbose:
                        fprintf_(options.fout,'\n### ecdfo: options.hess_approx not recognized\n\n')
                    info.flag=values.fail_on_argument
                    return M,pc,info
    return M,pc,info
