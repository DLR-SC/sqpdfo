# -*- coding: utf-8 -*-
from __future__ import division
#try:
from runtime import *
from numpy import inf, arange
from copy import copy

from sqplab_lsmult import sqplab_lsmult_
from ecdfo_optimality import ecdfo_optimality_
from ecdfo_iter_printout import ecdfo_iter_printout_
from ecdfo_solve_TR_bc import ecdfo_solve_TR_bc_
from ecdfo_augmX_evalf import ecdfo_augmX_evalf_
from bcdfo_augment_Y import bcdfo_augment_Y_
from ecdfo_swap_in_Y import ecdfo_swap_in_Y_
from bcdfo_computeP import bcdfo_computeP_
from bcdfo_gradP import bcdfo_gradP_
from bcdfo_projgrad import bcdfo_projgrad_
from ecdfo_computeHessian import ecdfo_computeHessian_
from bcdfo_poisedness_Y import bcdfo_poisedness_Y_
from bcdfo_repair_Y import bcdfo_repair_Y_
from ecdfo_find_smallf import ecdfo_find_smallf_
from bcdfo_include_in_Y import bcdfo_include_in_Y_
from numpy import array, zeros
#except ImportError:
    #from smop.runtime import *

def ecdfo_main_(func_=None,n_=None,nb_=None,mi_=None,me_=None,lm_=None,nitold_=None,nit_=None,i_xbest_=None,lb_=None,ub_=None,m_=None,X_=None,fX_=None,ciX_=None,ceX_=None,ind_Y_=None,QZ_=None,RZ_=None,delta_=None,cur_degree_=None,neval_=None,maxeval_=None,maxit_=None,fcmodel_=None,gx_=None,normgx_=None,show_errg_=None,pquad_=None,pdiag_=None,plin_=None,stallfact_=None,eps_rho_=None,Deltamax_=None,rep_degree_=None,epsilon_=None,verbose_=None,eta1_=None,eta2_=None,gamma1_=None,gamma2_=None,gamma3_=None,interpol_TR_=None,factor_CV_=None,Lambda_XN_=None,Lambda_CP_=None,factor_FPU_=None,factor_FPR_=None,Lambda_FP_=None,criterion_S_=None,criterion_FP_=None,criterion_CP_=None,mu_=None,theta_=None,eps_TR_=None,eps_L_=None,lSolver_=None,stratLam_=None,eps_current_=None,vstatus_=None,xstatus_=None,sstatus_=None,dstatus_=None,ndummyY_=None,sspace_save_=None,xspace_save_=None,xfix_=None,fxmax_=None,poised_model_=None,M_=None,kappa_ill_=None,kappa_th_=None,eps_bnd_=None,poised_=None,Y_radius_=None,c_=None,level_=None,whichmodel_=None,hardcons_=None,noisy_=None,scaleX_=None,scalefacX_=None,CNTsin_=None,shrink_Delta_=None,scale_=None,shift_Y_=None,info_=None,options_=None,values_=None,*args,**kwargs):
    """
% Realizes the optimization loop for the following algorithms:
% - quasi-Newton method
% - with trust regions.
%-----------------------------------------------------------------------

% Authors: Jean Charles Gilbert, INRIA.
%      and Anke Troeltzsch, DLR.
%
% Copyright 2008, 2009, INRIA. 2013, DLR.
%
% SQPLAB is distributed under the terms of the Q Public License version
% 1.0.
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
#    nargin = 89-[func,n,nb,mi,me,lm,nitold,nit,i_xbest,lb,ub,m,X,fX,ciX,ceX,ind_Y,QZ,RZ,delta,cur_degree,neval,maxeval,maxit,fcmodel,gx,normgx,show_errg,pquad,pdiag,plin,stallfact,eps_rho,Deltamax,rep_degree,epsilon,verbose,eta1,eta2,gamma1,gamma2,gamma3,interpol_TR,factor_CV,Lambda_XN,Lambda_CP,factor_FPU,factor_FPR,Lambda_FP,criterion_S,criterion_FP,criterion_CP,mu,theta,eps_TR,eps_L,lSolver,stratLam,eps_current,vstatus,xstatus,sstatus,dstatus,ndummyY,sspace_save,xspace_save,xfix,fxmax,poised_model,M,kappa_ill,kappa_th,eps_bnd,poised,Y_radius,c,level,whichmodel,hardcons,noisy,scaleX,scalefacX,CNTsin,shrink_Delta,scale,shift_Y,info,options,values].count(None)+len(args)

    func=copy(func_)
    n=copy(n_)
    nb=copy(nb_)
    mi=copy(mi_)
    me=copy(me_)
    lm=copy(lm_)
    nitold=copy(nitold_)
    nit=copy(nit_)
    i_xbest=copy(i_xbest_)
    lb=copy(lb_)
    ub=copy(ub_)
    m=copy(m_)
    X=copy(X_)
    fX=copy(fX_)
    ciX=copy(ciX_)
    ceX=copy(ceX_)
    ind_Y=copy(ind_Y_)
    QZ=copy(QZ_)
    RZ=copy(RZ_)
    delta=copy(delta_)
    cur_degree=copy(cur_degree_)
    neval=copy(neval_)
    maxeval=copy(maxeval_)
    maxit=copy(maxit_)
    fcmodel=copy(fcmodel_)
    gx=copy(gx_)
    normgx=copy(normgx_)
    show_errg=copy(show_errg_)
    pquad=copy(pquad_)
    pdiag=copy(pdiag_)
    plin=copy(plin_)
    stallfact=copy(stallfact_)
    eps_rho=copy(eps_rho_)
    Deltamax=copy(Deltamax_)
    rep_degree=copy(rep_degree_)
    epsilon=copy(epsilon_)
    verbose=copy(verbose_)
    eta1=copy(eta1_)
    eta2=copy(eta2_)
    gamma1=copy(gamma1_)
    gamma2=copy(gamma2_)
    gamma3=copy(gamma3_)
    interpol_TR=copy(interpol_TR_)
    factor_CV=copy(factor_CV_)
    Lambda_XN=copy(Lambda_XN_)
    Lambda_CP=copy(Lambda_CP_)
    factor_FPU=copy(factor_FPU_)
    factor_FPR=copy(factor_FPR_)
    Lambda_FP=copy(Lambda_FP_)
    criterion_S=copy(criterion_S_)
    criterion_FP=copy(criterion_FP_)
    criterion_CP=copy(criterion_CP_)
    mu=copy(mu_)
    theta=copy(theta_)
    eps_TR=copy(eps_TR_)
    eps_L=copy(eps_L_)
    lSolver=copy(lSolver_)
    stratLam=copy(stratLam_)
    eps_current=copy(eps_current_)
    vstatus=copy(vstatus_)
    xstatus=copy(xstatus_)
    sstatus=copy(sstatus_)
    dstatus=copy(dstatus_)
    ndummyY=copy(ndummyY_)
    sspace_save=copy(sspace_save_)
    xspace_save=copy(xspace_save_)
    xfix=copy(xfix_)
    fxmax=copy(fxmax_)
    poised_model=copy(poised_model_)
    M=copy(M_)
    kappa_ill=copy(kappa_ill_)
    kappa_th=copy(kappa_th_)
    eps_bnd=copy(eps_bnd_)
    poised=copy(poised_)
    Y_radius=copy(Y_radius_)
    c=copy(c_)
    level=copy(level_)
    whichmodel=copy(whichmodel_)
    hardcons=copy(hardcons_)
    noisy=copy(noisy_)
    scaleX=copy(scaleX_)
    scalefacX=copy(scalefacX_)
    CNTsin=copy(CNTsin_)
    shrink_Delta=copy(shrink_Delta_)
    scale=copy(scale_)
    shift_Y=copy(shift_Y_)
    info=copy(info_)
    options=copy(options_)
    values=copy(values_)

    old_delta=copy(delta)
    sigma=1
    rho_factor=0.3
    tau1=copy(gamma2)
    tau2=copy(gamma3)
    tau3=5
    constrained_pbl=copy(me)
    null_step=0
    merit=info.f + sigma * norm_(info.ce)
    msg='Unexpected message from ecdfo_main'
    m=size_(X,2)-1
    indfree=find_(vstatus == c.free)
    indfix=find_(vstatus >= c.fixed)
    nfix=length_(indfix)
    Y=X[indfree,ind_Y]
    x=X[indfree,i_xbest]
    if not isempty_(indfree):
        indfree=indfree.reshape(-1)
    if not isempty_(indfix):
        indfix=indfix.reshape(-1)
    n=size_(Y,1)
    fY=fX[ind_Y]
    fx=fX[i_xbest]
    itype=' '
    pc=0
    norms=0
    pred=0
    if mi > 0:
        ciY=copy(ciX[:,ind_Y])
    else:
        ciY=array([])
        gci=array([])
    if me > 0:
        ceY=copy(ceX[:,ind_Y])
    else:
        ceY=array([])
        gce=array([])
    radius_has_been_rejected=copy(False)
    while 1:

        if info.niter >= options.miter:
            info.flag=values.stop_on_max_iter
            ecdfo_iter_printout_(info,old_delta,norms,pc,itype,values,nb,mi,options,constrained_pbl,merit)
            return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
        if info.nsimul[1] >= options.msimul:
            info.flag=values.stop_on_max_simul
            ecdfo_iter_printout_(info,old_delta,norms,pc,itype,values,nb,mi,options,constrained_pbl,merit)
            return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
        xk=copy(x)
        lbounds=- inf * ones_(size_(x))
        ubounds=inf * ones_(size_(x))
        ilb=(abs(lb[indfree] - x) < 1e-05).reshape(-1)
        iub=(abs(ub[indfree] - x) < 1e-05).reshape(-1)
        lbounds[ilb]=lb[indfree[ilb]]
        ubounds[iub]=ub[indfree[iub]]
        lm,info=sqplab_lsmult_(x,lbounds,ubounds,info,options,values,nargout=2)
        feas,comp,info=ecdfo_optimality_(x,lm,lb[indfree],ub[indfree],info,options,nargout=3)
        if info.flag:
            return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
        info.glagn=norm_(info.glag,inf)
        info.feasn=norm_(feas,inf)
        info.compl=norm_(comp,inf)
        if (info.niter > 0) and (options.verbose >= 3):
            fprintf_(options.fout,'\nOptimality:\n')
            if constrained_pbl:
                fprintf_(options.fout,'  |grad Lag|      = %12.5e\n'%(info.glagn))
                fprintf_(options.fout,'  feasibility     = %12.5e\n'%(info.feasn))
            else:
                fprintf_(options.fout,' |grad f| = %12.5e\n'%(info.glagn))
        if ((info.glagn <= options.tol[0]) and (info.feasn <= options.tol[1]) and (info.compl <= options.tol[2])) or delta <= epsilon * 1e-05 or (pred == - 1.0):
            augment=rep_degree - cur_degree
            if (augment <= 0):
                poised,Y_radius=bcdfo_poisedness_Y_(QZ,RZ,Y,eps_L,x,lSolver,whichmodel,hardcons,lb,ub,indfree,stratLam,scale,shift_Y,nargout=2)
                poisedness_known=1
                errg=poised * Y_radius / factor_CV
                if options.verbose >= 3:
                    disp_('error on gradient before set improvement = ',num2str_(errg))
                if (((info.glagn <= options.tol[0]) and (info.feasn <= options.tol[1]) and (info.compl <= options.tol[2]) and errg <= epsilon) or delta <= epsilon * 1e-05) and level=='toplevel':
                    info.niter=info.niter + 1
                    itype='conv'
                    ecdfo_iter_printout_(info,old_delta,norms,pc,itype,values,nb,mi,options,constrained_pbl,merit)
                    info.flag=values.success
                    msg='Convergence in '+str(neval)+' evaluations of the objective function.'
                    return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
                else:
                    info.niter=info.niter + 1
                    ecdfo_iter_printout_(info,old_delta,norms,pc,itype,values,nb,mi,options,constrained_pbl,merit)
            if options.verbose >= 3:
                disp_('not immediately converged - improve set!')
            itype='impr'
            if (augment <= 0):
                eps_current=max_(mu * eps_current,epsilon)
            if (normgx <= epsilon):
                effective_FPR=1
            else:
                effective_FPR=copy(factor_FPR)
            if (augment > 0):
                itype='augm'
                if (info.glagn <= epsilon):
                    if (lSolver == 2):
                        delta=epsilon / sqrt_(n)
                        eps_current=epsilon / sqrt_(n)
                    else:
                        delta=copy(epsilon)
                        eps_current=copy(epsilon)
                ynew=- delta * ones_(n,1) + 2 * delta * rand_(n,1)
                cur_degree,QZ,RZ,Y,xbase,scale=bcdfo_augment_Y_(ynew,Y[:,0:cur_degree],whichmodel,shift_Y,delta,normgx,kappa_ill,nargout=6)
                ind_Y[cur_degree-1]=cur_degree-1
                if (hardcons):
                    ynew,improvement=bcdfo_find_new_yj_bc_(QZ,RZ,Y,cur_degree,delta,eps_L,xbase,lSolver,whichmodel,xl,xu,indfree,stratLam,scale,shift_Y,nargout=2)
                else:
                    ynew,improvement=bcdfo_find_new_yj_(QZ,RZ,Y,cur_degree,delta,eps_L,xbase,lSolver,whichmodel,scale,shift_Y,nargout=2)
                QZ,RZ,Y,xbase,scale=bcdfo_replace_in_Y_(QZ,RZ,ynew,Y,cur_degree,xbase,whichmodel,scale,shift_Y,delta,normgx,kappa_ill,nargout=5)
                replaced=array([cur_degree-1])
            else:
                if (info.glagn <= factor_CV * epsilon):
                    if (lSolver == 2):
                        radius=min_(delta / sqrt_(n),epsilon / sqrt_(n))
                    else:
                        radius=min_(delta,epsilon)
                else:
                    radius=max_(delta,eps_current)
                QZ,RZ,Y,replaced,poised,Y_radius,x,scale=bcdfo_repair_Y_(QZ,RZ,Y,radius,effective_FPR,Lambda_FP,Lambda_CP,eps_L,x,lSolver,whichmodel,hardcons,lb,ub,indfree,stratLam,scale,shift_Y,normgx,kappa_ill,nargout=8)
                if options.verbose >= 3:
                    disp_('improve interpolation set (in radius = ',num2str_(radius),') : replaced = ',num2str_(replaced),', poised = ',num2str_(poised),', Y_radius = ',num2str_(Y_radius))
            if (options.verbose >= 4):
                poised,Y_radius=bcdfo_poisedness_Y_(QZ,RZ,Y,eps_L,x,lSolver,whichmodel,hardcons,lb,ub,indfree,stratLam,scale,shift_Y,nargout=2)
                disp_(' poisedness(Y) = ',num2str_(poised))
            poised_model=1
            for i in range(0,length_(replaced)):
                j=replaced[i]
                m=m + 1
                xstatus[ind_Y[j]]=c.unused
                ind_Y[j]=m
                try:
                    xstatus[m]=c.inY
                except IndexError:
                    concatenate_([xstatus,[m]])
                X,fX,ciX,ceX,neval,xstatus,sstatus,dstatus,info,outdic=ecdfo_augmX_evalf_(func,Y[:,[j]],m,X,fX,ciX,ceX,nfix,xfix,indfix,indfree,fxmax,neval,xstatus,c.inY,sstatus,dstatus,scaleX,scalefacX,info,options,values,nargout=10)
                fY[j]=fX[m]
                if mi > 0:
                    ciY[:,j]=copy(info.ci.T)
                if me > 0:
                    ceY[:,j]=copy(info.ce.T)
                poised_model=0
                if msg=='Error':
                    if level=='toplevel':
                        disp_(msg)
                    return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
            i_xold=copy(i_xbest)
            x,fx,QZ,RZ,Y,fY,ciY,ceY,ind_Y,i_xbest,scale,info=ecdfo_find_smallf_(c,QZ,RZ,Y,fY,ciY,ceY,ind_Y,i_xbest,cur_degree,indfree,x,lb,ub,fx,dstatus,whichmodel,scale,shift_Y,delta,normgx,kappa_ill,sigma,info,nargout=12)
            fcmodel=bcdfo_computeP_(QZ,RZ,Y,concatenate_([fY.reshape(1,-1),ciY,ceY]),whichmodel,fcmodel,ind_Y,i_xold,m,gx,scale,shift_Y)
            gx=bcdfo_gradP_(fcmodel[[0],:],x,x,scale,shift_Y)
            normgx,_=bcdfo_projgrad_(n,x,gx,lb[indfree],ub[indfree])
            if mi > 0:
                gci=zeros_(mi,n)
                for i in range(0,mi):
                    gci[i,:]=bcdfo_gradP_(fcmodel[[1 + i],:],x,x,scale,shift_Y).T
            if me > 0:
                gce=zeros_(me,n)
                for i in range(0,me):
                    gce[i,:]=bcdfo_gradP_(fcmodel[[1 + mi + i],:],x,x,scale,shift_Y).T
            M,pc,info=ecdfo_computeHessian_(func,x,null_step,constrained_pbl,lm,M,n,me,mi,s,gx,gci,gce,info,options,values,fcmodel,Y,fY,ciY,ceY,sigma,scale,shift_Y,QZ,RZ,whichmodel,ind_Y,i_xbest,m,nargout=3)
            lbounds=- inf * ones_(size_(x))
            ubounds=inf * ones_(size_(x))
            ilb=(abs(lb[indfree] - x) < 1e-05).reshape(-1)
            iub=(abs(ub[indfree] - x) < 1e-05).reshape(-1)
            lbounds[ilb]=lb[indfree[ilb]]
            ubounds[iub]=ub[indfree[iub]]
            lm,info=sqplab_lsmult_(x,lbounds,ubounds,info,options,values,nargout=2)
            feas,comp,info=ecdfo_optimality_(x,lm,lb[indfree],ub[indfree],info,options,nargout=3)
            if info.flag:
                return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
            info.glagn=norm_(info.glag,inf)
            info.feasn=norm_(feas,inf)
            info.compl=norm_(comp,inf)
            if (info.niter > 0) and (options.verbose >= 3):
                fprintf_(options.fout,'\nOptimality:\n')
                if constrained_pbl:
                    fprintf_(options.fout,'  |grad Lag|      = %12.5e\n'%(info.glagn))
                    fprintf_(options.fout,'  feasibility     = %12.5e\n'%(info.feasn))
                else:
                    fprintf_(options.fout,' |grad f| = %12.5e\n'%(info.glagn))
            errg=poised * Y_radius / factor_CV
            if options.verbose >= 3:
                disp_('error on gradient after set improvement = ',num2str_(errg))
            if (info.glagn / factor_CV <= options.tol[0]) and (info.feasn / factor_CV <= options.tol[1]) and (info.compl / factor_CV <= options.tol[2]) and errg <= epsilon and cur_degree >= rep_degree and level=='toplevel':
                info.niter=info.niter + 1
                itype='conv'
                ecdfo_iter_printout_(info,old_delta,norms,pc,itype,values,nb,mi,options,constrained_pbl,merit)
                msg='Convergence in '+str(neval)+' evaluations of the objective function.'
                return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
            if options.verbose >= 3:
                disp_('not converged after improvement of interpolation set')
            delta=copy(radius)
        if not(radius_has_been_rejected):
            f0=info.f
            ce0=info.ce
            ce0n=norm_(ce0)
            merit0=f0 + sigma * ce0n
            prec_r=options.tol[1] / 10
            prec_t=options.tol[0] / 10
            if options.verbose >= 5:
                fprintf_(options.fout,'\nStep computation: merit = %12.5e\n'%(merit0))
            if options.verbose == 4:
                fprintf_(options.fout,'   radius     |r|      |t|      |s|     sigma     rho\n')
        info.niter=info.niter + 1
        ecdfo_iter_printout_(info,old_delta,norms,pc,itype,values,nb,mi,options,constrained_pbl,merit)
        if options.verbose >= 5:
            fprintf_(options.fout,'  Trust radius = %8.2e\n'%(delta))
        old_delta=copy(delta)
        xnew,deltaTR,rpred,active_r,active_t,lm_computed,lm,info=ecdfo_solve_TR_bc_(func,x,lb[indfree],ub[indfree],delta,mi,me,M,prec_r,prec_t,info,options,values,radius_has_been_rejected,lm,ceY,ciY,gx,nargout=8)
        if info.flag == values.fail_strange:
            return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
        s=xnew - xk
        x=copy(xk)
        norms=norm_(s)
        if options.verbose >= 3:
            fprintf_(options.fout,'  Full step:\n    |s| = %8.2e\n'%(norms))
        qcost=info.g.T.dot(s) + 0.5 * (s.T.dot(M.dot(s)))
        if rpred == 0:
            sigmab=0.0
        else:
            sigmab=qcost / ((1 - rho_factor) * rpred)
        if sigma < sigmab:
            sigma=max_(sigmab,1.5 * sigma)
            merit0=f0 + sigma * ce0n
        if options.verbose >= 4:
            fprintf_(options.fout,'  Penalty parameter = %8.2e (threshold %8.2e)\n'%(sigma,sigmab))
        if (interpol_TR == 1):
            gTs=gx.T.dot(s)
        xplus=x + s
        m=m + 1
        X,fX,ciX,ceX,neval,xstatus,sstatus,dstatus,info,outdic=ecdfo_augmX_evalf_(func,xplus,m,X,fX,ciX,ceX,nfix,xfix,indfix,indfree,fxmax,neval,xstatus,0,sstatus,dstatus,scaleX,scalefacX,info,options,values,nargout=10)
        if (info.flag):
            return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
        else:
            fxplus=copy(fX[m])
        if outdic:
            if outdic == 1:
                if options.verbose >= 5:
                    fprintf_(options.fout,'  Step rejected (out of an implicit domain)\n')
            else:
                if outdic == 2:
                    if options.verbose:
                        fprintf_(options.fout,'\n### ecdfo_main: the simulator wants to stop\n\n')
                    info.flag=values.stop_on_simul
                    return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
                else:
                    if options.verbose:
                        fprintf_(options.fout,'\n### ecdfo_main: error in the simulator (outdic = %0i)\n\n'%(outdic))
                    info.flag=values.fail_on_simul
                    return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
            itype='xfail'
            radius_has_been_rejected=copy(True)
            if options.verbose >= 5:
                fprintf_(options.fout,'  rho = %8.2e \n'%(rho))
            info.ce=ce0
            delta=tau1 * delta
            if options.verbose == 3 or options.verbose >= 5:
                fprintf_(options.fout,'  Step rejected due to failure in function evaluation\n')
        else:
            merit=info.f + sigma * norm_(info.ce)
            if options.verbose >= 3:
                fprintf_(options.fout,'  Merit function: %15.8e -> %15.8e\n'%(merit0,merit))
            ared=merit0 - merit
            pred=- qcost + sigma * rpred
            if pred < 0:
                if options.verbose:
                    fprintf_(options.fout,'\n### ecdfo_main: pred = %9.2e should be positive\n\n'%(pred))
                info.flag=values.fail_strange
                return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
            else:
                if pred == 0:
                    pred=- 1.0
                    disp_('### ecdfo_main : Warning : predicted reduction is 0 ###')
            rho=ared / pred
            if pred == - 1.0:
                rho=- 1.0
            if (rho >= eta1):
                succ=1
            else:
                succ=0
            if options.verbose == 4:   
                fprintf_(options.fout,'  %8.2e  %7.1e  %7.1e  %9.2e\n'%(delta,norms,sigma,rho))
#                fprintf_(options.fout,'  %8.2e  %7.1e  %7.1e  %7.1e  %7.1e  %9.2e\n'%(delta,norm_r,norm_(t),norms,sigma,rho))
            i_xold=copy(i_xbest)
            pos=-1
            if (rho >= eta1):
                if options.verbose >= 5:
                    fprintf_(options.fout,'  Step accepted (rho = %9.2e; ared = %9.2e, pred = %9.2e)\n'%(rho,ared,pred))
                if (merit >= merit0):
                    info.flag=values.fail_on_non_decrease
                    return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
                if (cur_degree < pquad or (whichmodel == 3 and cur_degree < pquad + pquad)):
                    cur_degree,QZ,RZ,Y,xbase,scale=bcdfo_augment_Y_(xplus,Y[:,0:cur_degree],whichmodel,shift_Y,delta,normgx,kappa_ill,nargout=6)
                    pos=copy(cur_degree)-1
                else:
                    QZ,RZ,Y,pos,x,scale=bcdfo_include_in_Y_(xplus,QZ,RZ,Y,arange(0,cur_degree),Lambda_XN,criterion_S,x,whichmodel,succ,scale,shift_Y,delta,normgx,kappa_ill,nargout=6)
                    if (pos >= 0):
                        xstatus[ind_Y[pos]]=c.unused
                if (pos >= 0):
                    itype='succ'
                    if (options.verbose >= 3):
                        disp_(' replacing/including interpolation point ',str(pos),' (successful)')
                    xstatus[m]=c.inY
                    try:
                        ind_Y[pos]=m
                        fY[pos]=copy(fxplus)
                    except IndexError:
                        ind_Y=concatenate_([ind_Y,[m]],axis=1)
                        fY=concatenate_([fY, [fxplus]],axis=1)
                    if mi > 0:
                        try:
                            ciY[:,pos]=copy(info.ci.T)
                        except IndexError:
                            ciY=concatenate_([ciY, info.ci],axis=1)
                    if me > 0:
                        try:
                            ceY[:,pos]=copy(info.ce.T)
                        except IndexError:
                            ceY=concatenate_([ceY, info.ce],axis=1)

                    QZ,RZ,Y,ind_Y,fY,ciY,ceY,x,scale=ecdfo_swap_in_Y_(0,pos,QZ,RZ,Y,ind_Y,fY,ciY,ceY,x,whichmodel,scale,shift_Y,delta,normgx,kappa_ill,nargout=9)
                    fx=copy(fxplus)
                    i_xbest=copy(m)
                    if (not shift_Y):
                        x=copy(Y[:,[0]])
                    poised_model=0
                    fcmodel=bcdfo_computeP_(QZ,RZ,Y,concatenate_([fY.reshape(1,-1),ciY,ceY]),whichmodel,fcmodel[[0],:],ind_Y,i_xold,m,gx,scale,shift_Y)
                    gx=bcdfo_gradP_(fcmodel[[0],:],x,x,scale,shift_Y)
                    normgx,_=bcdfo_projgrad_(n,x,gx,lb[indfree],ub[indfree])
                    if mi > 0:
                        gci=zeros_(mi,n)
                        for i in range(0,mi):
                            gci[i,:]=bcdfo_gradP_(fcmodel[[1 + i],:],x,x,scale,shift_Y).T
                    if me > 0:
                        gce=zeros_(me,n)
                        for i in range(0,me):
                            gce[i,:]=bcdfo_gradP_(fcmodel[[1 + mi + i],:],x,x,scale,shift_Y).T
                    if rho >= eta2:
                        if (active_r or active_t):
                            delta=delta * tau3
                        else:
                            delta=min_(max_(tau3 * norms,delta),Deltamax)
                    else:
                        if rho >= eta1:
                            if (active_r or active_t):
                                delta=delta * tau2
                            else:
                                delta=min_(max_(tau2 * norms,delta),Deltamax)
                    radius_has_been_rejected=copy(False)
                    if lm_computed == 0:
                        lbounds=- inf * ones_(size_(x))
                        ubounds=inf * ones_(size_(x))
                        ilb=(abs(lb[indfree] - x) < 1e-05).reshape(-1)
                        iub=(abs(ub[indfree] - x) < 1e-05).reshape(-1)
                        lbounds[ilb]=lb[indfree[ilb]]
                        ubounds[iub]=ub[indfree[iub]]
                        lm,info=sqplab_lsmult_(x,lbounds,ubounds,info,options,values,nargout=2)
                    M,pc,info=ecdfo_computeHessian_(func,x,null_step,constrained_pbl,lm,M,n,me,mi,s,gx,gci,gce,info,options,values,fcmodel,Y,fY,ciY,ceY,sigma,scale,shift_Y,QZ,RZ,whichmodel,ind_Y,i_xbest,m,nargout=3)
            if pred == - 1.0:
                pos=0
                rho=1
            if (rho < eta1) or (pos == -1):
                itype='repD,repF,repC,redD'
                itype='unsuc'
                radius_has_been_rejected=copy(True)
                if options.verbose == 3 or options.verbose >= 5:
                    fprintf_(options.fout,'  Step rejected (rho = %9.2e; ared = %9.2e, pred = %9.2e)\n'%(rho,ared,pred))
                if (((cur_degree < pquad) or (whichmodel == 3 and cur_degree < pquad + pquad)) and (rho < eta1)):
                    cur_degree,QZ,RZ,Y,xbase,scale=bcdfo_augment_Y_(xplus,Y[:,0:cur_degree],whichmodel,shift_Y,delta,normgx,kappa_ill,nargout=6)
                    if (options.verbose >= 3):
                        disp_(' including interpolation point ',str(cur_degree-1),' (augm)')
                    xstatus[m]=c.inY
                    try:
                        ind_Y[cur_degree-1]=m
                        fY[cur_degree-1]=copy(fxplus)
                    except IndexError:
                        ind_Y=concatenate_([ind_Y,[m]],axis=1)
                        fY=concatenate_([fY, [fxplus]],axis=1)
                    if mi > 0:
                        try:
                            ciY[:,cur_degree-1]=copy(info.ci.T)
                        except IndexError:
                            ciY=concatenate_([ciY, info.ci],axis=1)
                    if me > 0:
                        try:
                            ceY[:,cur_degree-1]=copy(info.ce.T)
                        except IndexError:
                            ceY=concatenate_([ceY, info.ce],axis=1)    
                    poised_model=0
                    fcmodel=bcdfo_computeP_(QZ,RZ,Y,concatenate_([fY.reshape(1,-1),ciY,ceY]),whichmodel,fcmodel,ind_Y,i_xold,m,gx,scale,shift_Y)
                    gx=bcdfo_gradP_(fcmodel[[0],:],x,x,scale,shift_Y)
                    normgx,_=bcdfo_projgrad_(n,x,gx,lb[indfree],ub[indfree])
                    if mi > 0:
                        gci=zeros_(mi,n)
                        for i in range(0,mi):
                            gci[i,:]=bcdfo_gradP_(fcmodel[[1 + i],:],x,x,scale,shift_Y).T
                    if me > 0:
                        gce=zeros_(me,n)
                        for i in range(0,me):
                            gce[i,:]=bcdfo_gradP_(fcmodel[[1 + mi + i],:],x,x,scale,shift_Y).T
                    itype='augm'
                    pos=copy(m)
                    if (shrink_Delta == 1 and delta > epsilon):
                        delta=gamma2 * delta
                if (cur_degree >= pquad or pos == -1):
                    if ((pos == -1) and (poised_model == 0 or delta <= eps_current)):
                        d=zeros(cur_degree)
                        if (rho >= eta1):
                            for j in range(0,cur_degree):
                                if (lSolver == 1):
                                    d[j]=norm_(Y[:,[j]] - xplus)
                                else:
                                    d[j]=norm_(Y[:,[j]] - xplus,inf)
                        else:
                            for j in range(1,cur_degree): #was a range(2,cur_degree) in matlab
                                if (lSolver == 1):
                                    d[j]=norm_(Y[:,[j]] - x)
                                else:
                                    d[j]=norm_(Y[:,[j]] - x,inf)
                        FPlength=factor_FPU * (1 + eps_TR) * delta
                        if (rho >= eta1):
                            criterion_FPn='weighted'
                        else:
                            criterion_FPn=copy(criterion_FP)
                        QZ,RZ,Y,pos,x,scale=bcdfo_include_in_Y_(xplus,QZ,RZ,Y,find_(d > FPlength),Lambda_FP,criterion_FPn,x,whichmodel,succ,scale,shift_Y,delta,normgx,kappa_ill,nargout=6)
                        if (pos >= 0):
                            itype='repF'
                            if (options.verbose >= 3):
                                disp_(' replacing interpolation point ',pos,' (far)')
                            xstatus[ind_Y[pos]]=c.unused
                            xstatus[m]=c.inY
                            ind_Y[pos]=m
                            fY[pos]=fxplus
                            if mi > 0:
                                ciY[:,pos]=info.ci.T
                            if me > 0:
                                ceY[:,pos]=info.ce.T
                            if (rho >= eta1):
                                QZ,RZ,Y,ind_Y,fY,ciY,ceY,x,scale=ecdfo_swap_in_Y_(0,pos,QZ,RZ,Y,ind_Y,fY,ciY,ceY,x,whichmodel,scale,shift_Y,delta,normgx,kappa_ill,nargout=9)
                                fx=copy(fxplus)
                                info.f=fx
                                if mi > 0:
                                    info.ci=copy(ciY[:,[0]])
                                if me > 0:
                                    info.ce=copy(ceY[:,[0]])
                                i_xbest=copy(m)
                                if (not shift_Y):
                                    x=Y[:,[0]]
                                poised_model=0
                                if (options.verbose >= 3):
                                    disp_(' swapped point to position 1')
                                itype='repFs'
                                delta=min_(max_(gamma3 * norms,delta),Deltamax)
                            else:
                                if (shrink_Delta == 1 and delta > epsilon):
                                    delta=gamma2 * delta
                            fcmodel=bcdfo_computeP_(QZ,RZ,Y,concatenate_([fY.reshape(1,-1),ciY,ceY]),whichmodel,fcmodel,ind_Y,i_xold,m,gx,scale,shift_Y)
                            gx=bcdfo_gradP_(fcmodel[[0],:],x,x,scale,shift_Y)
                            normgx,_=bcdfo_projgrad_(n,x,gx,lb[indfree],ub[indfree])
                            if mi > 0:
                                gci=zeros_(mi,n)
                                for i in range(0,mi):
                                    gci[i,:]=bcdfo_gradP_(fcmodel[[1 + i],:],x,x,scale,shift_Y).T
                            if me > 0:
                                gce=zeros_(me,n)
                                for i in range(0,me):
                                    gce[i,:]=bcdfo_gradP_(fcmodel[[1 + mi + i],:],x,x,scale,shift_Y).T
                        if (pos == -1):
                            if (rho >= eta1):
                                criterion_CPn='standard'
                            else:
                                criterion_CPn=copy(criterion_CP)
                            if (rho >= eta1):
                                Lambda_CPn=1e-15
                            else:
                                Lambda_CPn=copy(Lambda_CP)
                                d[0]=2 * FPlength
                            QZ,RZ,Y,pos,x,scale=bcdfo_include_in_Y_(xplus,QZ,RZ,Y,find_(d <= FPlength),Lambda_CPn,criterion_CPn,x,whichmodel,succ,scale,shift_Y,delta,normgx,kappa_ill,nargout=6)
                            if (pos >= 0):
                                itype='repC'
                                if (pos == 0):
                                    i_xold=ind_Y[1] #was a Y[2] in matlab
                                if (options.verbose >= 3):
                                    disp_(' replacing interpolation point ',str(pos-1),' (close)')
                                xstatus[ind_Y[pos]]=c.unused
                                xstatus[m]=c.inY
                                ind_Y[pos]=m
                                fY[pos]=fxplus
                                if mi > 0:
                                    ciY[:,pos]=info.ci.T
                                if me > 0:
                                    ceY[:,pos]=info.ce.T
                                if (rho >= eta1):
                                    QZ,RZ,Y,ind_Y,fY,ciY,ceY,x,scale=ecdfo_swap_in_Y_(0,pos,QZ,RZ,Y,ind_Y,fY,ciY,ceY,x,whichmodel,scale,shift_Y,delta,normgx,kappa_ill,nargout=9)
                                    fx=copy(fxplus)
                                    info.f=fx
                                    if mi > 0:
                                        info.ci=ciY[:,[0]]
                                    if me > 0:
                                        info.ce=ceY[:,[0]]
                                    i_xbest=copy(m)
                                    if (not shift_Y):
                                        x=copy(Y[:,[0]])
                                    poised_model=0
                                    if (options.verbose >= 3):
                                        disp_(' swapped point to position 1')
                                    itype='repCs'
                                    delta=min_(max_(gamma3 * norms,delta),Deltamax)
                                else:
                                    if (shrink_Delta == 1 and delta > epsilon):
                                        delta=gamma2 * delta
                                fcmodel=bcdfo_computeP_(QZ,RZ,Y,concatenate_([fY.reshape(1,-1),ciY,ceY]),whichmodel,fcmodel,ind_Y,i_xold,m,gx,scale,shift_Y)
                                gx=bcdfo_gradP_(fcmodel[[0],:],x,x,scale,shift_Y)
                                normgx,_=bcdfo_projgrad_(n,x,gx,lb[indfree],ub[indfree])
                                if mi > 0:
                                    gci=zeros_(mi,n)
                                    for i in range(0,mi):
                                        gci[i,:]=bcdfo_gradP_(fcmodel[[1 + i],:],x,x,scale,shift_Y).T
                                if me > 0:
                                    gce=zeros_(me,n)
                                    for i in range(0,me):
                                        gce[i,:]=bcdfo_gradP_(fcmodel[[1 + mi + i],:],x,x,scale,shift_Y).T
                    if (pos == -1):
                        if (options.verbose >= 3):
                            disp_(' decreasing the TR radius')
                        xstatus[m]=c.unused
                        if (interpol_TR == 1):
                            curvature=- pred - gTs
                            gam_inter=(eta2 - 1) * gTs / (fxplus - fx - gTs - eta2 * curvature)
                            delta=max_(gamma1,min_(gam_inter,gamma2)) * min_(delta,norms)
                        else:
                            delta=gamma2 * norms
                        itype='redD'
                        if (delta < stallfact * norm_(x) or delta < epsilon * 1e-05):
                            if options.verbose >= 2 and level=='toplevel':
                                ecdfo_iter_printout_(info,old_delta,norms,pc,itype,values,nb,mi,options,constrained_pbl,merit)
                                if (show_errg):
                                    disp_('************************************* Trust-region',' radius small *********************************')
                                else:
                                    disp_('******************************** Trust-region',' radius small ****************************')
                            msg='Algorithm stopped after '+str(neval)+' evaluations of the objective function because Delta small.'
                            info.flag=values.stop_on_small_trust_region
                            return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
                info.ce=ce0
                info.f=f0
                lbounds=- inf * ones_(size_(x))
                ubounds=inf * ones_(size_(x))
                ilb=(abs(lb[indfree] - x) < 1e-05).reshape(-1)
                iub=(abs(ub[indfree] - x) < 1e-05).reshape(-1)
                lbounds[ilb]=lb[indfree[ilb]]
                ubounds[iub]=ub[indfree[iub]]
                lm,info=sqplab_lsmult_(x,lbounds,ubounds,info,options,values,nargout=2)
                M,pc,info=ecdfo_computeHessian_(func,x,null_step,constrained_pbl,lm,M,n,me,mi,s,gx,gci,gce,info,options,values,fcmodel,Y,fY,ciY,ceY,sigma,scale,shift_Y,QZ,RZ,whichmodel,ind_Y,i_xbest,m,nargout=3)

    return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info

