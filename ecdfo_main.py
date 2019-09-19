# -*- coding: utf-8 -*-

from runtime import *
from numpy import inf, arange
from copy import copy

from ecdfo_compute_multiplier import ecdfo_compute_multiplier_
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
#from bcdfo_find_new_yj_bc import bcdfo_find_new_yj_bc_
from bcdfo_find_new_yj import bcdfo_find_new_yj_
from bcdfo_replace_in_Y import bcdfo_replace_in_Y_
from ecdfo_find_smallf import ecdfo_find_smallf_
from bcdfo_include_in_Y import bcdfo_include_in_Y_
from numpy import array, zeros, concatenate, zeros_like
import ecdfo_global_variables as glob


def ecdfo_main_(func_=None,n_=None,nb_=None,mi_=None,me_=None,lm_=None,nitold_=None,nit_=None,\
    i_xbest_=None,lb_=None,ub_=None,m_=None,X_=None,fX_=None,ciX_=None,ceX_=None,\
    ind_Y_=None,QZ_=None,RZ_=None,delta_=None,cur_degree_=None,neval_=None,\
    maxeval_=None,maxit_=None,fcmodel_=None,gx_=None,normgx_=None,show_errg_=None,\
    pquad_=None,pdiag_=None,plin_=None,stallfact_=None,eps_rho_=None,Deltamax_=None,\
    rep_degree_=None,epsilon_=None,verbose_=None,eta1_=None,eta2_=None,gamma1_=None,\
    gamma2_=None,gamma3_=None,interpol_TR_=None,factor_CV_=None,Lambda_XN_=None,\
    Lambda_CP_=None,factor_FPU_=None,factor_FPR_=None,Lambda_FP_=None,\
    criterion_S_=None,criterion_FP_=None,criterion_CP_=None,mu_=None,theta_=None,\
    eps_TR_=None,eps_L_=None,lSolver_=None,stratLam_=None,eps_current_=None,\
    vstatus_=None,xstatus_=None,sstatus_=None,dstatus_=None,ndummyY_=None,\
    sspace_save_=None,xspace_save_=None,xfix_=None,fxmax_=None,poised_model_=None,\
    M_=None,kappa_ill_=None,kappa_th_=None,eps_bnd_=None,poised_=None,Y_radius_=None,\
    c_=None,level_=None,whichmodel_=None,hardcons_=None,noisy_=None,scaleX_=None,\
    scalefacX_=None,CNTsin_=None,shrink_Delta_=None,scale_=None,shift_Y_=None,\
    info_=None,options_=None,values_=None,*args,**kwargs):
    ###############################################################################
    # Main optimization loop for ECDFO.
    ###############################################################################

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

    # Initialization
    
    old_delta  = copy(delta);     # for printing
    sigma      = 1;         # initial penalty parameter
    rho_factor = 0.3;       # it is imposed that pred/vred >= rho_factor*sigmab 
                            # (must be in (0,1))
    tau1       = copy(gamma2);    # trust radius reduction factor if out of domain
    tau2       = copy(gamma3);    # good trust radius augmentation factor 
                                  # when active and rho is >= eta1
    tau3       = 5;         # extra trust radius augmentation factor 
                            # when active and rho is >= eta2
                            
    nbr_slacks = glob.get_nbr_slacks()
    sl = glob.get_slacks()
    slplus = zeros_like(sl)

    constrained_pbl=copy(me)
    null_step=0
    ce = info.ce
    if nbr_slacks:
        merit=info.f + sigma * \
              norm_(ce - concatenate((zeros((len(ce)-nbr_slacks,1)),sl**2)))
    else:
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
        
    n=size_(x,1)
    fY=fX[ind_Y]
    fx=fX[i_xbest]
    itype=' '
    pc=0
    s=zeros((size_(x)))
    norms=0
    pred=0
    ciplus=array([])
    ceplus=array([])
    
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
        
    if options.final_degree == values.quadratic:
        pfinal = pquad
    elif options.final_degree == values.diagonal:
        pfinal = pdiag
    elif options.final_degree == values.linear:
        pfinal = plin

    ##########################################################################
    # Begin main loop
    ##########################################################################

    radius_has_been_rejected=copy(False)
    
    while 1:
        
        # Stop on counter.

        if info.niter >= options.miter:
            info.flag=values.stop_on_max_iter

            # Final printout

            ecdfo_iter_printout_(info,old_delta,norms,pc,itype,values,nb,mi,\
                                 options,constrained_pbl,merit)
            return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
            cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
            ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
            
        if info.nsimul[1] >= options.msimul:
        
            info.flag=values.stop_on_max_simul

            # Final printout

            ecdfo_iter_printout_(info,old_delta,norms,pc,itype,values,nb,mi,\
                                 options,constrained_pbl,merit)
            return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
            cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
            ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
            
        #   if radius_has_been_rejected == false:

        xk=copy(x)

       #--------------------------------------------------------------------
       # Check stopping criteria.
       #--------------------------------------------------------------------
    
       # Compute feasibility and complementarity.
       
        lbounds=- inf * ones_(size_(x))
        ubounds=inf * ones_(size_(x))
        ilb=(abs(lb[indfree] - x) < 1e-05).reshape(-1)
        iub=(abs(ub[indfree] - x) < 1e-05).reshape(-1)
        lbounds[ilb]=lb[indfree[ilb]]
        ubounds[iub]=ub[indfree[iub]]
        
        lm,info=\
        ecdfo_compute_multiplier_(x,lbounds,ubounds,info,options,values,nargout=2)
        
        feas,comp,info=\
        ecdfo_optimality_(x,lm,lb[indfree],ub[indfree],info,options,nargout=3)
        
        if info.flag:
            return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
            cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
            ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
            
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

        # Stop on convergence.

        if ((info.glagn <= options.tol_grad) and (info.feasn <= options.tol_feas) \
        and (info.compl <= options.tol_bnds)) or (delta<1e-10) or (pred == - 1.0):

            # check for accuracy and improve if necessary

            augment=rep_degree - cur_degree

            # No additional interpolation point is required for non-terminal repair.

            if (augment <= 0):

                #  Compute poisedness of the interpolation set Y

                poised,Y_radius=\
                bcdfo_poisedness_Y_(QZ,RZ,Y,eps_L,x,lSolver,whichmodel,hardcons,\
                                   lb,ub,indfree,stratLam,scale,shift_Y,nargout=2)
                                   
                poisedness_known=1

                #  Compute gradient accuracy 

                errg=poised * Y_radius / factor_CV
                
                if options.verbose >= 3:
                    disp_('error on gradient before set improvement = ',str(errg))
                    
                #  Check whether convergence tolerances are satisfied
                    
                if (((info.glagn <= options.tol_grad) \
                and (info.feasn <= options.tol_feas) \
                and (info.compl <= options.tol_bnds) and errg <= epsilon)) \
                and level=='toplevel':
                
                    info.niter=info.niter + 1
                    
                    itype='conv'

                    # Final printout

                    ecdfo_iter_printout_(info,old_delta,norms,pc,itype,values,\
                                        nb,mi,options,constrained_pbl,merit)
                    
                    info.flag=values.success
                    
                    msg='Convergence in '+str(neval)+\
                    ' evaluations of the objective function.'
                    
                    return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
                    cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
                    ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
                    
                else:
                
                    info.niter=info.niter + 1

                    #iteration printout
                    
                    ecdfo_iter_printout_(info,old_delta,norms,pc,itype,values,\
                                        nb,mi,options,constrained_pbl,merit)

            #  Not at a solution: improve the interpolation set.

            if options.verbose >= 3:
                disp_('not immediately converged - improve set!')
                
            itype='impr'

            #  Reduce eps_current if repair degree is reached.

            if (augment <= 0):
            
                if (pred == -1):
                
                    # in case of not yet converged feasibility problem
                    if (info.feasn > options.tol_feas):
                        eps_current = mu * delta
                        
                    # else improve in mu * eps_current
                    else:
                        eps_current = min_(mu * eps_current, epsilon)
                    
                else:
                    eps_current = epsilon

            #  Rebuild a poised model in the eps_current ball, and ...

            if (normgx <= epsilon or pred == -1):

                # ... remove far points (with strict criterion: farfact = 1, forcing all
                # new interpolation points to be within distance epsilon)

                effective_FPR=1
            else:

                # ... remove far points (with looser criterion)

                effective_FPR=copy(factor_FPR)

            #  One needs to add new interpolation points to reach the desired degree.
            #  (only entered if rep_degree is higher than linear!)

            if (augment > 0):
            
                itype='augm'
                
                #  If gradient small, find a new point in the epsilon-environment, 
                #  not in Delta (distinguish between infty-norm and 2-norm local solver)
         
                if (info.glagn <= epsilon):
                
                    if (lSolver == 2):
                        delta=epsilon / sqrt_(n)
                        eps_current=epsilon / sqrt_(n)
                    else:
                        delta=copy(epsilon)
                        eps_current=copy(epsilon)

                #  Pick a random interpolation point.

                ynew=- delta * ones_(n,1) + 2 * delta * rand_(n,1)
                
                # add the random point to the set Y
                
                cur_degree,QZ,RZ,Y,xbase,scale=\
                bcdfo_augment_Y_(ynew,Y[:,0:cur_degree],whichmodel,\
                                shift_Y,delta,normgx,kappa_ill,nargout=6)
                                
                ind_Y[cur_degree-1]=cur_degree-1
                
                #  Find an optimal new point to replace the random point

                if (hardcons):
                    ynew,improvement=\
                    bcdfo_find_new_yj_bc_(QZ,RZ,Y,cur_degree,delta,eps_L,xbase,\
                                         lSolver,whichmodel,xl,xu,indfree,stratLam,\
                                         scale,shift_Y,nargout=2)
                else:
                    ynew,improvement=\
                    bcdfo_find_new_yj_(QZ,RZ,Y,cur_degree,delta,eps_L,xbase,\
                                      lSolver,whichmodel,scale,shift_Y,nargout=2)
                    
                # replace the random point by the new point in Y 
                
                QZ,RZ,Y,xbase,scale=\
                bcdfo_replace_in_Y_(QZ,RZ,ynew,Y,cur_degree,xbase,whichmodel,\
                                   scale,shift_Y,delta,normgx,kappa_ill,nargout=5)
                                   
                replaced=array([cur_degree-1])

            #  The current interpolation set has the requested degree.

            else:

                #  If gradient small, repair in epsilon-radius, else in Delta
                #  (distinguish between infty-norm and 2-norm local solver)

                if pred == -1:
                    radius = min_(delta, eps_current)
                elif (info.glagn <= factor_CV * epsilon):
                    if (lSolver == 2):
                        radius=min_(delta / sqrt_(n),epsilon / sqrt_(n))
                    else:
                        radius=min_(delta, epsilon)
                else:
                    radius=max_(delta,eps_current)
                    
                #  Check that the trust-region radius has not become so small that a 
                #  repair step of this size will not be meaningful.
                
                if (radius < stallfact * norm_(x) or radius < epsilon * 1e-5):
                
                    msg='Algorithm stopped after '+str(neval)+\
                        ' evaluations of the objective function because Delta small.'
                    
                    info.flag=values.stop_on_small_trust_region
                    
                    # final printout
                    
                    ecdfo_iter_printout_(info,radius,norms,pc,itype,values,nb,mi,\
                          options,constrained_pbl,merit)
                    
                    return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
                    cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
                    ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info

                # repair Y
                
                QZ,RZ,Y,replaced,poised,Y_radius,x,scale=\
                bcdfo_repair_Y_(QZ,RZ,Y,radius,effective_FPR,Lambda_FP,Lambda_CP,\
                               eps_L,x,lSolver,whichmodel,hardcons,lb,ub,indfree,\
                               stratLam,scale,shift_Y,normgx,kappa_ill,nargout=8)
                               
                if options.verbose >= 3:
                    disp_('improve interpolation set (in radius = ',str(radius),\
                          ') : replaced = ',str(replaced),', poised = ',str(poised),\
                          ', Y_radius = ',str(Y_radius))
                    
            if (options.verbose >= 4):
                poised,Y_radius=\
                bcdfo_poisedness_Y_(QZ,RZ,Y,eps_L,x,lSolver,whichmodel,\
                hardcons,lb,ub,indfree,stratLam,scale,shift_Y,nargout=2)
                
                disp_(' poisedness(Y) = ',str(poised))
                
            poised_model=1

            #  Compute the corresponding function values.

            for i in range(0,length_(replaced)):
            
                j=int(replaced[i])

                #  Set index of new point and update status of the old point

                m=m + 1
                xstatus[ind_Y[j]]=c.unused
                ind_Y[j]=m
                
                try:
                    xstatus[m]=c.inY
                except IndexError:
                    concatenate_([xstatus,[m]],axis=1)

                #  Update X and evaluate function

                X,fX,ciX,ceX,neval,xstatus,sstatus,dstatus,info,outdic=\
                ecdfo_augmX_evalf_(func,Y[:,[j]],m,X,fX,ciX,ceX,nfix,xfix,\
                                   indfix,indfree,fxmax,neval,xstatus,c.inY,\
                                   sstatus,dstatus,scaleX,scalefacX,info,\
                                   options,values,nargout=10)
                                   
                fY[j]=fX[m]
                
                if mi > 0:
                    ciY[:,j]=copy(ciX[:,m].T)
                if me > 0:
                    ceY[:,j]=copy(ceX[:,m].T)
                    
                poised_model=0
                
                if msg=='Error':
                
                    if level=='toplevel':
                        disp_(msg)
                        
                    return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
                    cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
                    ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info

            #  Move to the best point found, if different from x.

            i_xold=copy(i_xbest)
            x,fx,QZ,RZ,Y,fY,ciY,ceY,ind_Y,i_xbest,scale,info=\
            ecdfo_find_smallf_(c,QZ,RZ,Y,fY,ciY,ceY,ind_Y,i_xbest,cur_degree,\
                              indfree,x,lb,ub,fx,dstatus,whichmodel,scale,shift_Y,\
                              delta,normgx,kappa_ill,sigma,info,nargout=12)

            #  Compute new model(s).

            fcmodel=bcdfo_computeP_(QZ,RZ,Y,concatenate_([fY.reshape(1,-1),ciY,ceY]),\
                                   whichmodel,fcmodel,ind_Y,i_xold,m,gx,scale,shift_Y)

            #  Compute the gradient(s) of the new model(s).

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

            # Update Hessian approximation (and gradients in info.g, info.ai, info.ae)
            
            M,pc,info=\
            ecdfo_computeHessian_(func,x,null_step,constrained_pbl,lm,M,n,me,mi,s,\
                                 gx,gci,gce,info,options,values,fcmodel,Y,fY,ciY,\
                                 ceY,sigma,scale,shift_Y,QZ,RZ,whichmodel,ind_Y,\
                                 i_xbest,m,nargout=3)

            #  Terminate if the solution has been found.
      
            # Compute feasibility and complementarity.
            
            lbounds=- inf * ones_(size_(x))
            ubounds=inf * ones_(size_(x))
            ilb=(abs(lb[indfree] - x) < 1e-05).reshape(-1)
            iub=(abs(ub[indfree] - x) < 1e-05).reshape(-1)
            lbounds[ilb]=lb[indfree[ilb]]
            ubounds[iub]=ub[indfree[iub]]
            
            # compute multiplier
            
            lm,info=\
            ecdfo_compute_multiplier_(x,lbounds,ubounds,info,options,values,nargout=2)
            
            # compute optimality
            
            feas,comp,info=\
            ecdfo_optimality_(x,lm,lb[indfree],ub[indfree],info,options,nargout=3)
            
            if info.flag:
                return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
                cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
                ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
                
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

            # Stop on convergence.

            errg=poised * Y_radius / factor_CV
            
            if options.verbose >= 3:
                disp_('error on gradient after set improvement = ',str(errg))
                
            # check whether convergence tolerances are satisfied   
             
            if (info.glagn / factor_CV <= options.tol_grad)\
                and (info.feasn / factor_CV <= options.tol_feas)\
                and (info.compl / factor_CV <= options.tol_bnds) and errg <= epsilon\
                and cur_degree >= rep_degree and level=='toplevel':
                
                info.niter=info.niter + 1
                
                itype='conv'                
                
                ecdfo_iter_printout_(info,old_delta,norms,pc,itype,values,nb,mi,\
                      options,constrained_pbl,merit)
                
                msg='Convergence in '+str(neval)+\
                ' evaluations of the objective function.'
                
                return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
                cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
                ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
                
            if options.verbose >= 3:
                disp_('not converged after improvement of interpolation set')

            #  Reset the radius to a multiple of ||gx||.

            delta=copy(radius)

        #if ( augment <= 0 ):
        #    delta = min( min( theta * normgx, epsilon ), Deltamax )

        #-----------------------------------------------------------------------
        # Globalization
        #-----------------------------------------------------------------------
  
        if not(radius_has_been_rejected):
            f0      = info.f;           # memorize f at the current iterate, 
                                        # useful in the sufficient decrease condition
            if nbr_slacks:     
                ce0 = copy(info.ce) - \
                      concatenate((zeros((len(info.ce)-nbr_slacks,1)),sl**2))
            else:
                ce0 = info.ce;          # memorize ce at the current iterate in
                                         # case of step rejection
            ce0n    = norm_(ce0);
            merit0  = f0 + sigma * ce0n;  # inital value of the merit function, 
                                          # for the given sigma
            prec_r  = options.tol_feas/10;  # initial precision for the restoration
                                            # problem
            prec_t  = options.tol_grad/10;  # initial precision for the tangent 
                                            # problem

            if options.verbose >= 5:
                fprintf_(options.fout,'\nStep computation: merit = %12.5e\n'%(merit0))
                
            if options.verbose == 4:
                fprintf_(options.fout,'   radius     |r|      |t|      |s|',\
                '     sigma     rho\n')
                
        info.niter=info.niter + 1

        #-----------------------------------------------------------------------
        # Iteration printout
        #-----------------------------------------------------------------------
        
        ecdfo_iter_printout_(info,old_delta,norms,pc,itype,values,nb,mi,\
              options,constrained_pbl,merit)
        
        if options.verbose >= 5:
            fprintf_(options.fout,'  Trust radius = %8.2e\n'%(delta))

        # -----------------------------------------------
        # Compute new trial point in TR and inside bounds
        # -----------------------------------------------

        # Compute Hessian of the interpolation model
  
        #M = bcdfo_hessP( fcodel(1,:), x, x, scale, shift_Y );
   
        # solve trust-region subproblem in delta
        
        old_delta=copy(delta)

        xnew,deltaTR,rpred,active_r,active_t,lm_computed,lm,info,slplus=\
        ecdfo_solve_TR_bc_(func,x,lb,ub,delta,mi,me,M,\
                           prec_r,prec_t,info,options,values,\
                           radius_has_been_rejected,lm,ceY,ciY,gx,indfree,\
                           nargout=9)
        
        # check for error

        if info.flag == values.fail_unexpected:
            return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
            cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
            ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info

        # ---------------------
        # Full step computation
        # ---------------------

        #s=xnew - xk
        #x=copy(xk)
        s=xnew-x
        
        if nbr_slacks:
            norms=norm_(concatenate((s,slplus-sl)))
        else:
            norms=norm_(s)
        
        if options.verbose >= 3:
            fprintf_(options.fout,'  Full step:\n    |s| = %8.2e\n'%(norms))

        # -----------------------------
        # Compute new penalty parameter
        # -----------------------------

        qcost=info.g.T.dot(s) + 0.5 * (s.T.dot(M.dot(s)))
        
        if rpred == 0:
            sigmab=0.0
        else:
            sigmab=qcost / ((1 - rho_factor) * rpred)
            
        if sigma < sigmab:        
            sigma=max_(sigmab,1.5 * sigma)
            
            if sigma > 1e+250:
                fprintf_(options.fout,'\n### ecdfo_main: Penalty parameter (sigma): %15.8e '\
                     'is too big\n\n'%(sigma))
                     
                info.flag=values.fail_unexpected
                
                return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
                cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
                ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
            
            # re-evaluate the merit function at x, since sigma has changed
            merit0=f0 + sigma * ce0n
            
        if options.verbose >= 4:
            fprintf_(options.fout,'  Penalty parameter = %8.2e (threshold %8.2e)\n'\
            %(sigma,sigmab))

        # -----------------------------------------
        # Evaluate function values at the new point
        # -----------------------------------------

        if (interpol_TR == 1):   # save for computing backtracking interpolation of TR
            gTs=gx.T.dot(s)
            
        # compute new trial point 
            
        xplus=x + s

        #  Set index of new point

        m=m + 1

        # Include point in X and evaluate f
        # (xstatus(m) is set to 0 but is updated later on)
        
        X,fX,ciX,ceX,neval,xstatus,sstatus,dstatus,info,retval=\
        ecdfo_augmX_evalf_(func,xplus,m,X,fX,ciX,ceX,nfix,xfix,indfix,indfree,\
                          fxmax,neval,xstatus,0,sstatus,dstatus,scaleX,\
                          scalefacX,info,options,values,nargout=10)
                         
        if (info.flag):
        
            return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
            cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
            ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
            
        else:
        
            fxplus=copy(fX[m])
            
            if ceX.any():
                ceX[ceX>=1e25] = 10*max_(ceX[ceX<1e25])
                ceX[ceX<=-1e25] = 10*min_(ceX[ceX>-1e25])
                ceplus = copy(array([ceX[:,m]]).T)

        # ---------------
        # Step validation
        # ---------------
        
        if retval:
        
            if retval == 1:
            
                if options.verbose >= 5:
                    fprintf_(options.fout,'  Step rejected (out of an implicit domain)\n')
                    
            elif retval == 2:
            
                if options.verbose > 1:
                    fprintf_(options.fout,'\n### ecdfo_main: evaluation of the function stopped\n\n')
                    
                info.flag=values.stop_on_simul
                
                return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
                cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
                ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
                
            else: # unexpected
            
                if options.verbose > 1:
                    fprintf_(options.fout,'\n### ecdfo_main: error during evaluation of the function')
                    
                info.flag=values.fail_on_simul
                
                return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
                cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
                ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
                
            itype='xfail'
            
            radius_has_been_rejected=copy(True)

            # Recover ce at the current iterate

            info.ce=ce0

            # Shrink trust region radius

            #delta = tau1*norms;
            delta=tau1 * delta
            
            if options.verbose == 3 or options.verbose >= 5:
                fprintf_(options.fout,'  Step rejected due to failure in function'\
                                ' evaluation\n')
                
        else:
        
          # -------------------------------------
          # Compute merit function and ratio rho
          # -------------------------------------
      
            if nbr_slacks:
                merit=fxplus + sigma * norm_(ceplus - \
                      concatenate((zeros((len(ceplus)-nbr_slacks,1)),slplus**2))) 
            else:
                merit=fxplus + sigma * norm_(ceplus)
                
            if np.isinf(merit):
                fprintf_(options.fout,'  Merit function: %15.8e -> %15.8e\n'\
                     %(merit0,merit))
                     
                info.flag=values.fail_unexpected
                
                return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
                cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
                ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info

            if options.verbose >= 3:
                fprintf_(options.fout,'  Merit function: %15.8e -> %15.8e\n'\
                %(merit0,merit))
                
            ared=merit0 - merit
            pred=- qcost + sigma * rpred
            
            if rpred < 0:
                # inaccurate model is assumed, thus an improvement step is performed
                # at the next iteration
                pred = -1.0
            
            if pred < 0:  
                # should not occur, since here s ~= 0
                # an improvement step is performed at the next iteration
                pred=-1.0
                            
                if options.verbose >= 3:
                    fprintf_(options.fout,'\n### ecdfo_main: pred = %9.2e should be positive\n\n'%(pred))
                
            elif pred == 0:            
                # here, stationarity is assumed but model maybe inaccurate to state 
                # convergence thus, an improvement step is performed next iteration
                pred=- 1.0
                
                if options.verbose >= 3:
                    disp_('### ecdfo_main : Warning : predicted reduction is 0 ###')
                    
            rho=ared / pred
            
            if pred == - 1.0:
                rho=- 1.0
                
            if (rho >= eta1):
                succ=1
            else:
                succ=0
                
            if options.verbose == 4:   
                fprintf_(options.fout,'  %8.2e  %7.1e  %7.1e  %9.2e\n'\
                %(delta,norms,sigma,rho))
#               fprintf_(options.fout,'  %8.2e  %7.1e  %7.1e  %7.1e  %7.1e  %9.2e\n'\
#               %(delta,norm_r,norm_(t),norms,sigma,rho))

            ###################################################################
            # Include the new point in the interpolation set Y 
            ###################################################################

            i_xold=copy(i_xbest)   # save to compute min-frob-norm model
            pos=-1
            
            # --------------------------------------
            # Successful iteration (accept the step)
            # --------------------------------------
            
            if (rho >= eta1):
            
                if options.verbose >= 3:
                    fprintf_(options.fout,'  Step accepted (rho = %9.2e;'\
                    ' ared = %9.2e, pred = %9.2e)\n'%(rho,ared,pred))
                    
                if (merit >= merit0):
                    info.flag=values.fail_on_non_decrease
                    return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
                    cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
                    ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info
         
                #  Augment interpolation set if not fully quadratic yet

                if (cur_degree < pfinal or (whichmodel == 3 \
                and cur_degree < pfinal + pfinal)):
                    cur_degree,QZ,RZ,Y,xbase,scale=\
                    bcdfo_augment_Y_(xplus,Y[:,0:cur_degree],whichmodel,shift_Y,\
                                    delta,normgx,kappa_ill,nargout=6)
                    
                    pos=copy(cur_degree)-1
                    
                #  Include xplus in the interpolation set, by replacing another point if
                #  the model is already fully quadratic.
                
                else:
                    QZ,RZ,Y,pos,x,scale=\
                    bcdfo_include_in_Y_(xplus,QZ,RZ,Y,arange(0,cur_degree),Lambda_XN,\
                                       criterion_S,x,whichmodel,succ,scale,shift_Y,\
                                       delta,normgx,kappa_ill,nargout=6)
                   
                    if (pos >= 0):
                        xstatus[ind_Y[pos]]=c.unused
                        
                #  If xplus could/should be included in the interpolation set

                if (pos >= 0):
                
                    itype='succ'
                    
                    if (options.verbose >= 3):
                        disp_(' replacing/including interpolation point ',\
                        str(pos),' (successful)')
                        
                    xstatus[m]=c.inY
                    
                    try:
                        ind_Y[pos]=m
                        fY[pos]=copy(fxplus)
                    except IndexError:
                        ind_Y=concatenate_([ind_Y,[m]],axis=1)
                        fY=concatenate_([fY, [fxplus]],axis=1)
                                                    
                    if me > 0:
                        try:
                            info.ce=copy(ceplus)
                            ceY[:,pos]=copy(info.ce.T)
                        except IndexError:
                            info.ce=copy(ceplus)
                            ceY=concatenate_([ceY, info.ce],axis=1)
  
                        if nbr_slacks:
                            # move slack variable away from zero if inequality
                            # value gets above zero (here trial use of 0.1)
                            
                            for i in range(0,nbr_slacks):
                                if slplus[i]==0 and info.ce[me-nbr_slacks+i]>0.01:
                                   slplus[i] = sqrt_(info.ce[me-nbr_slacks+i])
                                    
                            sl = slplus
                            glob.set_slacks(slplus)
                            
                    #  Move it in the first position, redefining the base point.

                    QZ,RZ,Y,ind_Y,fY,ciY,ceY,x,scale=\
                    ecdfo_swap_in_Y_(0,pos,QZ,RZ,Y,ind_Y,fY,ciY,ceY,x,whichmodel,\
                                    scale,shift_Y,delta,normgx,kappa_ill,nargout=9)
                                    
                    fx=copy(fxplus)
                    i_xbest=copy(m)
                    
                    if (not shift_Y):
                        x=copy(Y[:,[0]])
                        
                    poised_model=0
                    
                    #  Compute the associated polynomial interpolation models.

                    fcmodel=\
                    bcdfo_computeP_(QZ,RZ,Y,concatenate_(\
                                 [fY.reshape(1,-1),ciY,ceY]),whichmodel,\
                                 fcmodel[[0],:],ind_Y,i_xold,m,gx,scale,shift_Y)
                  
                    #  Compute model gradients for objective and constraint
                    #  functions.
                    
                    gx=bcdfo_gradP_(fcmodel[[0],:],x,x,scale,shift_Y)
                    normgx,_=bcdfo_projgrad_(n,x,gx,lb[indfree],ub[indfree])
                    
                    if mi > 0:
                        gci=zeros_(mi,n)
                        for i in range(0,mi):
                            gci[i,:]=\
                            bcdfo_gradP_(fcmodel[[1 + i],:],x,x,scale,shift_Y).T
                    if me > 0:
                        gce=zeros_(me,n)
                        for i in range(0,me):
                            gce[i,:]=\
                            bcdfo_gradP_(fcmodel[[1 + mi + i],:],x,x,scale,shift_Y).T
                            
                    #  Update the trust-region radius.
 
                    if rho >= eta2:
                        if (active_r or active_t) and delta<5.0:
                            delta=delta * tau3
                        else:
                            delta=min_(max_(tau2 * norms,delta),Deltamax)
                    else:
                        if rho >= eta1:
                            if (active_r or active_t):
                                delta=delta * tau2
                            else:
                                delta=min_(max_(tau2 * norms,delta),Deltamax)
                                
                    radius_has_been_rejected=copy(False)
                    
                    # Re-compute Lagrange multipliers (if not already done)

                    if lm_computed == 0:
                        lbounds=- inf * ones_(size_(x))
                        ubounds=inf * ones_(size_(x))
                        ilb=(abs(lb[indfree] - x) < 1e-05).reshape(-1)
                        iub=(abs(ub[indfree] - x) < 1e-05).reshape(-1)
                        lbounds[ilb]=lb[indfree[ilb]]
                        ubounds[iub]=ub[indfree[iub]]
                        
                        lm,info=\
                        ecdfo_compute_multiplier_(x,lbounds,ubounds,info,options,\
                                                 values,nargout=2)
                        
                    # Update Hessian approximation
                    
                    M,pc,info=\
                    ecdfo_computeHessian_(func,x,null_step,constrained_pbl,lm,M,\
                                         n,me,mi,s,gx,gci,gce,info,options,values,\
                                         fcmodel,Y,fY,ciY,ceY,sigma,scale,shift_Y,\
                                         QZ,RZ,whichmodel,ind_Y,i_xbest,m,nargout=3)
                                        
            # ---------------------------------------
            # Unsuccessful iteration - step rejection
            # ---------------------------------------
    
            # model accuracy is questionable - go to next iteration 
            # and choose new interpolation set
            
            if pred == - 1.0:
                pos=0
                rho=1
                
            if (rho < eta1) or (pos == -1):
            
                itype='repD,repF,repC,redD'
                itype='unsuc'
                
                radius_has_been_rejected=copy(True)
                
                if options.verbose == 3 or options.verbose >= 5:
                    fprintf_(options.fout,'  Step rejected (rho = %9.2e;'\
                    ' ared = %9.2e, pred = %9.2e)\n'%(rho,ared,pred))
                    
                #  The model is not fully quadratic yet: add (if possible)
                #  the new point to the interpolation set and recompute the model.
                
                if (((cur_degree < pfinal) or (whichmodel == 3 \
                    and cur_degree < pfinal + pfinal)) and (rho < eta1)):
                    
                    cur_degree,QZ,RZ,Y,xbase,scale=\
                    bcdfo_augment_Y_(xplus,Y[:,0:cur_degree],whichmodel,shift_Y,\
                                    delta,normgx,kappa_ill,nargout=6)
                                    
                    if (options.verbose >= 3):
                        disp_(' including interpolation point ',\
                        str(cur_degree-1),' (augm)')
                        
                    # Update status and position of the new point

                    xstatus[m]=c.inY
                    
                    try:
                        ind_Y[cur_degree-1]=m
                        fY[cur_degree-1]=copy(fxplus)
                    except IndexError:
                        ind_Y=concatenate_([ind_Y,[m]],axis=1)
                        fY=concatenate_([fY, [fxplus]],axis=1)
                        
                    if mi > 0:
                        try:
                            ciY[:,cur_degree-1]=copy(ciplus.T)
                        except IndexError:
                            ciY=concatenate_([ciY, ciplus],axis=1)
                    if me > 0:
                        try:
                            ceY[:,cur_degree-1]=copy(ceplus.T)
                        except IndexError:
                            ceY=concatenate_([ceY, ceplus],axis=1)  
                              
                    poised_model=0
                    
                    #  Compute new model(s).

                    fcmodel=bcdfo_computeP_(QZ,RZ,Y,concatenate_([fY.reshape(1,-1),\
                                           ciY,ceY]),whichmodel,fcmodel,ind_Y,i_xold,\
                                           m,gx,scale,shift_Y)
                    
                    #  Compute the gradient(s) of the new model(s).
                    
                    gx=bcdfo_gradP_(fcmodel[[0],:],x,x,scale,shift_Y)
                    normgx,_=bcdfo_projgrad_(n,x,gx,lb[indfree],ub[indfree])

                    if mi > 0:
                        gci=zeros_(mi,n)
                        for i in range(0,mi):
                            gci[i,:]=\
                            bcdfo_gradP_(fcmodel[[1 + i],:],x,x,scale,shift_Y).T
                    if me > 0:
                        gce=zeros_(me,n)
                        for i in range(0,me):
                            gce[i,:]=\
                            bcdfo_gradP_(fcmodel[[1 + mi + i],:],x,x,scale,shift_Y).T
                    itype='augm'
                    pos=copy(m)
                    
                    #  Shrink trust region in unsuccessful iteration

                    if (shrink_Delta == 1 and delta > epsilon):
                        #delta = gamma2 * norms;

                        delta=gamma2 * delta

                #  Enter if the model is already fully quadratic *or* 
                #  xplus could not yet be included in the set. 
                #  The decision to include xplus here depends on possibly eliminating 
                #  another point.
                
                if (cur_degree >= pfinal or pos == -1):
                
                    if ((pos == -1) and (poised_model == 0 or delta <= eps_current)):
                    
                        #  Compute the distance of the interpolation points to the
                        #  current iterate. (Distinguish between the badly conditioned
                        #  successful and the unsuccessful case!)
                        
                        d=zeros(cur_degree)
                        
                        if (rho >= eta1):
                            for j in range(0,cur_degree):
                                if (lSolver == 1):
                                    d[j]=norm_(Y[:,[j]] - xplus)
                                else:
                                    d[j]=norm_(Y[:,[j]] - xplus,inf)
                        else:
                            for j in range(1,cur_degree): 
                                if (lSolver == 1):
                                    d[j]=norm_(Y[:,[j]] - x)
                                else:
                                    d[j]=norm_(Y[:,[j]] - x,inf)
                                    
                        #  Compute the basic distance used to define far/close points.

                        FPlength=factor_FPU * (1 + eps_TR) * delta
                        
                        #  Replace a far interpolation point.

                        if (rho >= eta1):
                            criterion_FPn='weighted'   # use weighted measure, 
                                                       # not furthest point
                        else:
                            criterion_FPn=copy(criterion_FP)
                            
                        QZ,RZ,Y,pos,x,scale=\
                        bcdfo_include_in_Y_(xplus,QZ,RZ,Y,find_(d > FPlength),\
                                           Lambda_FP,criterion_FPn,x,whichmodel,\
                                           succ,scale,shift_Y,delta,normgx,\
                                           kappa_ill,nargout=6)
                                          
                        if (pos >= 0):
                        
                            itype='repF'
                            
                            if (options.verbose >= 3):
                                disp_(' replacing interpolation point ',pos,' (far)')
                                
                            #  Update status and position of the new point

                            xstatus[ind_Y[pos]]=c.unused
                            xstatus[m]=c.inY
                            ind_Y[pos]=copy(m)
                            fY[pos]=copy(fxplus)
                            
                            if mi > 0:
                                ciY[:,pos]=copy(ciplus.T)
                            if me > 0:
                                ceY[:,pos]=copy(ceplus.T)
                                
                            #  Swap points if included a successful point

                            if (rho >= eta1):
                            
                                QZ,RZ,Y,ind_Y,fY,ciY,ceY,x,scale=\
                                ecdfo_swap_in_Y_(0,pos,QZ,RZ,Y,ind_Y,fY,ciY,ceY,\
                                                x,whichmodel,scale,shift_Y,delta,\
                                                normgx,kappa_ill,nargout=9)
                                
                                fx=copy(fxplus)
                                info.f=fx
                                
                                if mi > 0:
                                    info.ci=copy(ciY[:,[0]])
                                if me > 0:
                                    info.ce=copy(ceY[:,[0]])
                                    
                                    if nbr_slacks:
                                        # set slacks if any
                                        for i in range(0,nbr_slacks):
                                            if slplus[i]==0 \
                                            and info.ce[me-nbr_slacks+i]>0.01:
                                                slplus[i] = sqrt_(info.ce[\
                                                            me-nbr_slacks+i])
                                        sl = slplus
                                        glob.set_slacks(slplus)
                                    
                                i_xbest=copy(m)
                                
                                if (not shift_Y):
                                    x=Y[:,[0]]
                                    
                                poised_model=0
                                
                                if (options.verbose >= 3):
                                    disp_(' swapped point to position 1')
                                    
                                itype='repFs'
                                
                                #  Update the trust-region radius.

                                delta=min_(max_(gamma3 * norms,delta),Deltamax)
                                
                            else:
                                #  Shrink trust region in unsuccessful iteration

                                if (shrink_Delta == 1 and delta > epsilon):
                                
                                    #delta = gamma2 * norms;
                                    delta=gamma2 * delta
                                    
                            #  Compute the associated polynomial interpolation model.

                            fcmodel=bcdfo_computeP_(QZ,RZ,Y,concatenate_\
                                                   ([fY.reshape(1,-1),ciY,ceY]),\
                                                   whichmodel,fcmodel,ind_Y,i_xold,\
                                                   m,gx,scale,shift_Y)
                                                   
                            #  Compute the gradient(s) of the new model(s).
                            
                            gx=bcdfo_gradP_(fcmodel[[0],:],x,x,scale,shift_Y)
                            normgx,_=bcdfo_projgrad_(n,x,gx,lb[indfree],ub[indfree])
                            
                            if mi > 0:
                                gci=zeros_(mi,n)
                                for i in range(0,mi):
                                    gci[i,:]=\
                                    bcdfo_gradP_(fcmodel[[1 + i],:],x,x,scale,shift_Y).T
                            if me > 0:
                                gce=zeros_(me,n)
                                for i in range(0,me):
                                    gce[i,:]=\
                                    bcdfo_gradP_(fcmodel[[1 + mi +\
                                               i],:],x,x,scale,shift_Y).T
                                    
                        # Replace a close interpolation point.

                        if (pos == -1):
                        
                            if (rho >= eta1):
                                criterion_CPn='standard'   # find best improvement
                            else:
                                criterion_CPn=copy(criterion_CP)
                                
                            if (rho >= eta1):
                                Lambda_CPn=1e-15  # try hard to include a successful 
                                                  # point
                            else:
                                Lambda_CPn=copy(Lambda_CP)
                                d[0]=2 * FPlength   # excludes the current iterate
                                
                            QZ,RZ,Y,pos,x,scale=\
                            bcdfo_include_in_Y_(xplus,QZ,RZ,Y,find_(d <= FPlength),\
                                               Lambda_CPn,criterion_CPn,x,whichmodel,\
                                               succ,scale,shift_Y,delta,normgx,\
                                               kappa_ill,nargout=6)
                            
                            if (pos >= 0):
                            
                                itype='repC'      
                                
                                #  Safeguard i_x for frobenius model type 4 when 
                                #  replacing point 1

                                if (pos == 0):
                                    i_xold=ind_Y[1] 
                                    
                                if (options.verbose >= 3):
                                    disp_(' replacing interpolation point ',\
                                    str(pos-1),' (close)')
                                    
                                #  Update status and position of the new point

                                xstatus[ind_Y[pos]]=c.unused
                                xstatus[m]=c.inY
                                ind_Y[pos]=copy(m)
                                fY[pos]=copy(fxplus)
                                
                                if mi > 0:
                                    ciY[:,pos]=copy(ciplus.T)
                                if me > 0:
                                    ceY[:,pos]=copy(ceplus.T)
                                    
                                #  Swap points if included a successful point

                                if (rho >= eta1):
                                    QZ,RZ,Y,ind_Y,fY,ciY,ceY,x,scale=\
                                    ecdfo_swap_in_Y_(0,pos,QZ,RZ,Y,ind_Y,fY,ciY,\
                                                    ceY,x,whichmodel,scale,shift_Y,\
                                                    delta,normgx,kappa_ill,nargout=9)
                                                    
                                    fx=copy(fxplus)
                                    info.f=fx
                                    
                                    if mi > 0:
                                        info.ci=ciY[:,[0]]
                                        
                                    if me > 0:
                                        info.ce=ceY[:,[0]]
                                        
                                        if nbr_slacks:
                                            # set slacks if any
                                            for i in range(0,nbr_slacks):
                                                if slplus[i]==0 \
                                                and info.ce[me-nbr_slacks+i]>0.01:
                                                    slplus[i] = sqrt_(info.ce[\
                                                                me-nbr_slacks+i])
                                            sl = slplus
                                            glob.set_slacks(slplus)
                                        
                                    i_xbest=copy(m)
                                    
                                    if (not shift_Y):
                                        x=copy(Y[:,[0]])
                                        
                                    poised_model=0
                                    
                                    if (options.verbose >= 3):
                                        disp_(' swapped point to position 1')
                                        
                                    itype='repCs'
                                    
                                    #  Update the trust-region radius.

                                    delta=min_(max_(gamma3 * norms,delta),Deltamax)
                                    
                                else:
                                    #  Shrink trust region in unsuccessful iteration

                                    if (shrink_Delta == 1 and delta > epsilon):
                                        #delta = gamma2 * norms;
                                        delta=gamma2 * delta
                                        
                                #  Compute associated polynomial interpolation model.

                                fcmodel=bcdfo_computeP_(QZ,RZ,Y,concatenate_\
                                                       ([fY.reshape(1,-1),ciY,ceY]),\
                                                       whichmodel,fcmodel,ind_Y,\
                                                       i_xold,m,gx,scale,shift_Y)
                                                       
                                #  Compute the gradient(s) of the new model(s).
                                                       
                                gx=bcdfo_gradP_(fcmodel[[0],:],x,x,scale,shift_Y)
                                normgx,_=bcdfo_projgrad_(n,x,gx,lb[indfree],ub[indfree])
                                
                                if mi > 0:
                                    gci=zeros_(mi,n)
                                    for i in range(0,mi):
                                        gci[i,:]=\
                                        bcdfo_gradP_(fcmodel[[1 +\
                                                    i],:],x,x,scale,shift_Y).T
                                if me > 0:
                                    gce=zeros_(me,n)
                                    for i in range(0,me):
                                        gce[i,:]=\
                                        bcdfo_gradP_(fcmodel[[1 + mi +\
                                                    i],:],x,x,scale,shift_Y).T
                                        
                    # Decrease the radius.

                    if (pos == -1):
                    
                        if (options.verbose >= 3):
                            disp_(' decreasing the TR radius')
                            
                        #  Set status of the new point

                        xstatus[m]=c.unused
                        
                        #  Compute new trust-region radius

                        if (interpol_TR == 1):
                            curvature=- pred - gTs
                            gam_inter=(eta2 - 1) * gTs / (fxplus - fx - gTs - \
                                      eta2 * curvature)
                            delta=max_(gamma1,min_(gam_inter,gamma2)) * \
                                  min_(delta,norms)
                        else:
#                           delta     = gamma2 * delta;
                            delta=gamma2 * norms
                            
                        itype='redD'
                        
                        #  Check that the trust-region radius has not become so small  
                        #  that a step of this size will not be significant.
                        
                        if (delta < stallfact * norm_(x) or delta < epsilon * 1e-5):
                            msg='Algorithm stopped after '+str(neval)+\
                           ' evaluations of the objective function because Delta small.'
                                
                            info.flag=values.stop_on_small_trust_region
                            
                            # final printout
                            
                            ecdfo_iter_printout_(info,delta,norms,pc,itype,\
                                      values,nb,mi,options,constrained_pbl,merit)
                                      
                            return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,\
                            eps_current,cur_degree,fcmodel,gx,normgx,vstatus,\
                            xstatus,sstatus,dstatus,M,ndummyY,sspace_save,\
                            xspace_save,msg,CNTsin,neval,lm,info
                            
                # Recover ce at the current iterate

                #info.ce=ce0
                info.f=f0

                # Recompute Lagrange multiplier

                lbounds=- inf * ones_(size_(x))
                ubounds=inf * ones_(size_(x))
                ilb=(abs(lb[indfree] - x) < 1e-05).reshape(-1)
                iub=(abs(ub[indfree] - x) < 1e-05).reshape(-1)
                lbounds[ilb]=lb[indfree[ilb]]
                ubounds[iub]=ub[indfree[iub]]
                lm,info=\
                ecdfo_compute_multiplier_(x,lbounds,ubounds,info,options,values,\
                                         nargout=2)
                
                #  Compute / update Hessian
                
                M,pc,info=\
                ecdfo_computeHessian_(func,x,null_step,constrained_pbl,lm,M,n,me,mi,\
                                     s,gx,gci,gce,info,options,values,fcmodel,Y,fY,\
                                     ciY,ceY,sigma,scale,shift_Y,QZ,RZ,whichmodel,\
                                     ind_Y,i_xbest,m,nargout=3)
                
    return nit,i_xbest,x,fx,m,X,fX,ciX,ceX,ind_Y,delta,eps_current,\
    cur_degree,fcmodel,gx,normgx,vstatus,xstatus,sstatus,dstatus,M,\
    ndummyY,sspace_save,xspace_save,msg,CNTsin,neval,lm,info

