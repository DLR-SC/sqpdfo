# -*- coding: utf-8 -*-
import sys
#try:
from runtime import *
#except ImportError:
#from smop.runtime import *

import helper
from bcdfo_solve_TR_MS_bc import *
#from bcdfo_solve_TR_MS_bc_decorator import *
#from blls import *
#from ecdfo_check_cond import *
from ecdfo_check_convex import *
from sqplab_lsmult import *
from sqplab_tcg import *
from copy import copy
from ecdfo_global_variables import get_threshold, get_check_condition

def ecdfo_solve_TR_bc_(simul=None,x_=None,lb=None,ub=None,delta_=None,mi=None,me=None,M=None,prec_r=None,prec_t=None,info_=None,options=None,values=None,radius_has_been_rejected=None,lm_=None,ceY=None,ciY=None,gx=None,*args,**kwargs):
    """
%  function to compute the new SQP-trust-region step inside delta and subject
%  to simple bounds
%
%  inputs:
%
% 
%  outputs:
%
    """
#    varargin = cellarray(args)
#    nargin = 18-[simul,x,lb,ub,delta,mi,me,M,prec_r,prec_t,info,options,values,radius_has_been_rejected,lm,ceY,ciY,gx].count(None)+len(args)

#    violated=None
    info_r=helper.dummyUnionStruct()
    
    x=copy(x_)
    delta=copy(delta_)
    info=copy(info_)
    lm=copy(lm_)

    #If somewhere we have done set_check_condition (in the tests for instance we have set_check_condition(0)), then
    #we get this value, otherwise we take '1' by default.
    try:
        check_condition=get_check_condition()
    except:
        check_condition=1      
    threshold=get_threshold()
    lm_computed=0
    n=length_(x)
    I=eye_(n)
    xi=1.0
    xnew=copy(x)
    active_r=0
    active_t=0
    rpred=0
    norm_r=0
    x_fix= []
    glocal=copy(gx)
    delta_min=1e-08
    plevel_r=0
    if options.verbose >= 5:
        plevel_r=1
    plevel_t=0
    if options.verbose >= 5:
        plevel_t=1
    if me + mi == 0:
        lb_r=lb[0:n] - x
        ub_r=ub[0:n] - x
        stratLam=1
        s,_lambda,norms,value,gplus,nfact,neigd,msg=bcdfo_solve_TR_MS_bc_(glocal,M,lb_r,ub_r,delta,1e-07,stratLam,nargout=8)
        xnew=x + s
        rpred=0
        active_r=0
        if norm_(s) < delta:
            active_t=0
        else:
            active_t=1
        return xnew,delta,rpred,active_r,active_t,lm_computed,lm,info
    constraints=info.ce
    gconstraints=info.ae
    x_active=zeros_(size_(x))
    look_for_active_bounds=1
    if look_for_active_bounds == 1:
        gradlag=copy(glocal)
        bounds=(abs(lb - x) < 1e-05) + (abs(ub - x) < 1e-05)
        A=find_(bounds[0:n])
        gradlag[A]=gradlag[A] + lm[A]
        if me > 0:
            gradlag=gradlag + info.ae.T.dot(lm[n + mi:n + mi + me])
        for i in range(0,n):
            if (x[i] - gradlag[i] <= lb[i]) and (abs(x[i] - lb[i]) < 1e-05):
                x_active[i]=1
                if options.verbose >= 3:
                    disp_('lb ',num2str_(i),' is initially active')
#                constraints=array([[constraints],[0]])
#                gconstraints=array([[gconstraints],[I[i,:]]])
                constraints=concatenate_([constraints,array([[0]])])
                gconstraints=concatenate_([gconstraints,I[[i],:]])
            if (x[i] - gradlag[i] >= ub[i]) and (abs(x[i] - ub[i]) < 1e-05):
                x_active[i]=1
                if options.verbose >= 3:
                    disp_('ub ',num2str_(i),' is initially active')
#                constraints=array([[constraints],[0]])
#                gconstraints=array([[gconstraints],[I[i,:]]])
                constraints=concatenate_([constraints,array([[0]])])
                gconstraints=concatenate_([gconstraints,I[[i],:]])
        if options.verbose >= 3 and sum_(x_active) > 0:
            xactive=x_active.T
    finished=0
    iter_active=0
    while finished == 0:

        iter_active=iter_active + 1
        if options.verbose >= 3:
            disp_('********* Loop over active set ********* iteration ',num2str_(iter_active),' *********')
        if options.verbose >= 3:
            fprintf_(options.fout,'  Restoration step:\n')
        delta_r=xi * delta
        sol=gconstraints.T.dot(constraints)
        if iter_active == 1 and norm_(sol) > 1e-14:
            g1=gconstraints.T.dot(constraints)
            H1=gconstraints.T.dot(gconstraints)
            if check_condition:
                cthreshold=1e+16
                H1,badcond=ecdfo_check_cond_(H1,cthreshold,options,nargout=2)
            check_convex=1
            if check_convex:
                H1=ecdfo_check_convex_(H1,options)
            lb_r=lb[0:n] - x
            ub_r=ub[0:n] - x
            stratLam=1
            r,_lambda,norms,value,gplus,nfact,neigd,msg=bcdfo_solve_TR_MS_bc_(g1,H1,lb_r,ub_r,delta_r,prec_r,stratLam,nargout=8)
            info_r.prec=sqrt_(gplus.T.dot(gplus))
            xr=x + r
            norm2_r=r.T.dot(r)
            norm_r=sqrt_(norm2_r)
            if msg=='error' or msg=='limit':
                info_r.flag=- 1
            elif msg=='boundary' or (delta_r - norm_r < 1e-08):
                info_r.flag=1
            else:
                info_r.flag=0
            rpred_m=norm_(constraints) - norm_(constraints + gconstraints.dot(r))
            rpred=rpred + rpred_m
            if rpred < 0:
                if options.verbose >= 3:
                    fprintf_(options.fout,'\n### ecdfo_solve_TR_bc: rpred = %9.2e should not be negative\n\n'%(rpred))
            active_r=(info_r.flag == 1) or (info_r.flag == 2)
            if options.verbose >= 5:
                if - 1 == info_r.flag:
                    fprintf_(options.fout,'    max of %0i iterations reached\n'%(20 * me))
                elif 0 == info_r.flag:
                    fprintf_(options.fout,'    precision is less than required tolerance %8.2e\n'%(prec_r))
                elif 1 == info_r.flag:
                    fprintf_(options.fout,'    TR boundary is reached\n')
                elif 2 == info_r.flag:
                    fprintf_(options.fout,'    negative curvature direction encountered\n')
                fprintf_(options.fout,'    |r|   = %8.2e\n'%(norm_r))
                fprintf_(options.fout,'    rpred = %8.2e\n'%(rpred))
        else:
            if options.verbose >= 5:
                fprintf_(options.fout,'    unchanged\n')
            r=zeros_(size_(x))
            active_r=copy(False)
            xr=copy(x)
            norm2_r=0.0
            norm_r=0.0
        if options.verbose == 3:
            disp_('r = (',num2str_(r.T),')')
            disp_('delta_r = ',num2str_(delta_r),', norm_r = ',num2str_(norm_r))
        if options.verbose >= 3:
            fprintf_(options.fout,'  Tangent step:\n')
        delta_t=copy(delta)
        deg_freedom=n - length_(constraints)
        if deg_freedom > 0:
            Z_=null_(gconstraints)
            M_t=Z_.T.dot(M.dot( Z_))
            g_t=Z_.T.dot((glocal + M.dot(r)))
            u,info_t=sqplab_tcg_(M_t,- g_t,delta_t,20 * (n - me),prec_t,plevel_t,options.fout,nargout=2)
            t=Z_.dot(u)
            active_t=(info_t.flag == 1) or (info_t.flag == 2)
            if options.verbose >= 5:
                if - 1 == info_t.flag:
                    fprintf_(options.fout,'    max of %0i iterations reached\n'%(20 * (n - me)))
                elif 0 == info_t.flag:
                    fprintf_(options.fout,'    precision is less than required tolerance %8.2e\n'%(prec_t))
                elif 1 == info_t.flag:
                    fprintf_(options.fout,'    TR boundary is reached\n')
                elif 2 == info_t.flag:
                    fprintf_(options.fout,'    negative curvature direction encountered\n')
                fprintf_(options.fout,'    |t| = %8.2e\n'%(norm_(t)))
        else:
            t=zeros_(1,n).T
            active_t=0
        if options.verbose == 3:
            disp_('t = (',num2str_(t.T),')')
            disp_('delta_t = ',num2str_(delta_t),', norm_t = ',num2str_(norm_(t)))
            disp_('delta ',num2str_(delta),', norm_s = ',num2str_(norm_(r + t)))
        xnew=x + r + t
        x_active=zeros_(size_(xnew))
        x_viol=zeros_(size_(xnew))
        for i in range(0,n):
            if (xnew[i] - lb[i] < - threshold):
                x_viol[i]=-(i+1) #we have to be careful there because we have (i+1) in python which equals i in matlab ! It is tricky to use it both as assignment value and index.
#                x_fix=array([x_fix,i])
                x_fix=concatenate_([x_fix, array([i])], axis=1)   #here we do not have to increment i because x_fix is an index array
                if options.verbose >= 3:
                    disp_('lb ',int2str_(i),' is violated')
#                constraints=array([[constraints],[0]])
#                gconstraints=array([[gconstraints],[I[i,:]]])
                constraints=concatenate_([constraints,array([[0]])])
                gconstraints=concatenate_([gconstraints,I[[i],:]])
                violated=1
                break
            elif (abs(xnew[i] - lb[i]) < 1e-07):
                x_active[i]=-(i+1)#we have to be careful there because we have (i+1) in python which equals i in matlab ! It is tricky to use it both as assignment value and index.
                if options.verbose >= 3:
                    disp_('lb ',int2str_(i),' is active')
            else:
                if (xnew[i] - ub[i] > threshold):
                    x_viol[i]=i+1 #we have to be careful there because we have (i+1) in python which equals i in matlab ! It is tricky to use it both as assignment value and index.
#                        x_fix=array([x_fix,i])
                    x_fix=concatenate_([x_fix, array([i])], axis=1) #here we do not have to increment i because x_fix is an index array
                    if options.verbose >= 3:
                        disp_('ub ',int2str_(i),' is violated')
#                   constraints=array([[constraints],[0]])
#                   gconstraints=array([[gconstraints],[I[i,:]]])
                    constraints=concatenate_([constraints,array([[0]])])
                    gconstraints=concatenate_([gconstraints,I[[i],:]])
                    violated=1
                    break
                elif (abs(xnew[i] - ub[i]) < 1e-07):
                    x_active[i]=i+1 #we have to be careful there because we have (i+1) in python which equals i in matlab ! It is tricky to use it both as assignment value and index.
                    if options.verbose >= 3:
                        disp_('ub ',int2str_(i),' is active')
        if sum_(x_viol) == 0:
            violated=0
            if options.verbose >= 3:
                disp_('no new bound violated')
        else:
            if options.verbose >= 3:
                disp_(x_fix)
        if norm_(r + t) <= 1e-16:
            if (iter_active >= 10 * n or delta < delta_min):
                if options.verbose >= 3:
                    disp_('### ecdfo_solve_TR_bc: active-set iteration limit exceeded ###')
                return xnew,delta,rpred,active_r,active_t,lm_computed,lm,info
            lbounds=- inf * ones_(size_(x))
            ubounds=inf * ones_(size_(x))
            ilb=abs(lb - xnew) < 1e-05
            iub=abs(ub - xnew) < 1e-05
            lbounds[ilb]=lb[ilb]
            ubounds[iub]=ub[iub]
            lm,info=sqplab_lsmult_(xnew,lbounds,ubounds,info,options,values,nargout=2)
            min_lm,ind_min_lm=min_(lm[x_fix],nargout=2)
            if options.verbose >= 3:
                disp_('smallest Lagrange multiplier (for the bounds) = ',num2str_(min_lm))
            if min_lm < 0:
                if options.verbose >= 3:
                    disp_('Zero step but not converged - release one bound!!')
                constraints=constraints[0:me]
                gconstraints=gconstraints[0:me,:]
                xa=find_(x_active < 0).reshape(-1)
                if length_(xa) == length_(x_fix):
                    if (np.sort(xa) == np.sort(x_fix)).all():
                        for i in range(0,length_(x_fix)):
                            if i != ind_min_lm:
#                                constraints=array([[constraints],[0]])
#                                gconstraints=array([[gconstraints],[I[x_fix[i],:]]])
                                constraints=concatenate_([constraints,array([[0]])])
                                gconstraints=concatenate_([gconstraints,I[[x_fix[i]],:]])
#                        x_fix[ind_min_lm]=[]
                        x_fix = np.delete(x_fix, ind_min_lm) 
                        finished=1
                    else:
                        x_fix=array([])
                        for i in range(0,length_(find_(x_active < 0))):
#                            constraints=array([[constraints],[0]])
#                            gconstraints=array([[gconstraints],[I[xa[i],:]]])
                            constraints=concatenate_([constraints,array([[0]])])
                            gconstraints=concatenate_([gconstraints,I[[xa[i]],:]])
#                            x_fix=array([x_fix,xa[i]])
                            x_fix=concatenate_([x_fix,[xa[i]]],axis=1)
                else:
                    x_fix=array([])
                    for i in range(0,length_(find_(x_active < 0))):
#                        constraints=array([[constraints],[0]])
#                        gconstraints=array([[gconstraints],[I[xa[i],:]]])
                        constraints=concatenate_([constraints,array([[0]])])
                        gconstraints=concatenate_([gconstraints,I[[xa[i]],:]])
#                        x_fix=array([x_fix,xa[i]])
                        x_fix=concatenate_([x_fix,[xa[i]]],axis=1)
                if options.verbose >= 3:
                    x_fix
            else:
                if options.verbose >= 3:
                    disp_('Zero step and converged - go back to TR-loop...')
                finished=1
        else:
            if violated == 0:
                if options.verbose >= 3:
                    disp_('non zero feasible step - go back to TR-loop...')
                finished=1
            else:
                if options.verbose >= 3:
                    disp_('non zero infeasible step - continue finding correct active set...')
        if violated == 1:
            tstep=copy(t)
            if options.verbose >= 3:
                disp_('shorten tangential step')
#            aT=array([[eye_(n)],[- eye_(n)]])
            aT=concatenate_([eye_(n),- eye_(n)])
            aTx=aT.dot(xr)
#            alpha=array([[ub[1:n]],[- lb[1:n]]])
            alpha=concatenate_([ub[0:n],- lb[0:n]])
            divisor=aT.dot(tstep)
            ratio=(alpha - aTx) / divisor
            minratio=min_(ratio[divisor > 0])
            if (minratio < 0):
                minratio=0.
            else:
                minratio=float(minratio)
            tstep=minratio*tstep
            x=xr + tstep
            step=r + tstep
            if options.verbose == 3:
                disp_('t = (',num2str_(tstep.T),')')
                disp_('delta_t = ',num2str_(delta_t),', norm_t = ',num2str_(norm_(tstep)))
                disp_('delta ',num2str_(delta),', norm_s = ',num2str_(norm_(r + tstep)))
            glocal=glocal + M .dot(step)
            for i in range(0,me):
                gconstraints[i,:]=gconstraints[[i],:] + (M .dot( step)).T
            delta=delta - norm_(step)

    return xnew,delta,rpred,active_r,active_t,lm_computed,lm,info

