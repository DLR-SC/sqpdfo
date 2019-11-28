# -*- coding: utf-8 -*-
from runtime import *
import helper
from bcdfo_solve_TR_MS_bc import *
from minq import Minq

from sqpdfo_check_convex import *
from sqpdfo_compute_multiplier import *
from sqpdfo_truncated_cg import *
from copy import copy
from numpy import zeros_like
import sqpdfo_global_variables as glob

def sqpdfo_solve_TR_bc_(simul=None,x=None,lb=None,ub=None,delta=None,mi=None,me=None,\
    M=None,prec_r=None,prec_t=None,info_=None,options=None,values=None,\
    radius_has_been_rejected=None,lm=None,ceY=None,ciY=None,gx=None,indfree=None,\
    *args,**kwargs):
#############################################################################
#
#  function to compute the new SQP-trust-region step inside delta and subject
#  to simple bounds
#
############################################################################
    
    info_r=helper.dummyUnionStruct()
    
    info=copy(info_)

    nbr_slacks = glob.get_nbr_slacks()
    sl = glob.get_slacks()
 
    #If somewhere we have done set_check_condition (in the tests for instance 
    #we have set_check_condition(0)), then
    #we get this value, otherwise we take '1' by default.
    try:
        check_condition=glob.get_check_condition()
    except:
        check_condition=1    
          
    threshold=glob.get_threshold()
    lm_computed=0
    n=length_(x)
    I=eye_(n+nbr_slacks)
    xi=1.0  # trust radius contrictor (must be in (0,1))
    xnew=copy(x)
    active_r=0
    active_t=0
    rpred=0
    norm_r=0
    x_fix= []
    glocal=copy(gx)
    delta_min=1e-08
    slnew=zeros_like(sl)

    plevel_r=0 # printing level for the restoration step (0; nothing, >0: something)
    if options.verbose >= 5: 
        plevel_r=1 
    plevel_t=0  # printing level for the tangent step (0; nothing, >0: something)
    if options.verbose >= 5:
        plevel_t=1

#####################################################################

# in the unconstrained case (only simple bounds)

#####################################################################

    if me + mi == 0:
    
        lb_r=lb[0:n] - x
        ub_r=ub[0:n] - x
        
        stratLam=1
        
        s,_lambda,norms,value,gplus,nfact,neigd,msg=\
        bcdfo_solve_TR_MS_bc_(glocal,M,lb_r,ub_r,delta,1e-07,stratLam,options,nargout=8)
        
        xnew=x + s
        
        rpred=0
        active_r=0
        
        if norm_(s) < delta:
            active_t=0
        else:
            active_t=1
            
        return xnew,delta,rpred,active_r,active_t,lm_computed,lm,info,slnew


######################################################################

# in the constrained case

#####################################################################

    constraints=info.ce
    gconstraints=info.ae
    
    # modify constraints in case of slack variables due to inequalities
    
    if nbr_slacks:
        constraints = constraints - concatenate((zeros((me-nbr_slacks,1)),sl**2))
        hh1 = concatenate((zeros((me-nbr_slacks,nbr_slacks)),-2*diag(sl.T[0])))
        gconstraints = concatenate((gconstraints,hh1),axis=1)

    # Compute initial active bounds

    x_active=zeros_(size_(x))
    look_for_active_bounds=1
    
    if look_for_active_bounds == 1:
        gradlag=copy(glocal)
        bounds=(abs(lb[indfree] - x) < 1e-05) + (abs(ub[indfree] - x) < 1e-05)
        A=find_(bounds[0:n])
        gradlag[A]=gradlag[A] + lm[A]
        
        if me > 0:
            if nbr_slacks:
                gradlag = concatenate((gradlag,zeros((nbr_slacks,1))))
                gradlag = gradlag + gconstraints.T.dot(lm[n:n+me])
            else:
                gradlag=gradlag + gconstraints.T.dot(lm[n:n+me])
            
        for i in range(0,n):
            if (x[i] - gradlag[i] <= lb[i]) and (abs(x[i] - lb[i]) < 1e-05):
                x_active[i]=1
                if options.verbose >= 3:
                    disp_('lb ',str(i),' is initially active')
                    
                # append active bound to the equality constraints

                constraints=concatenate_([constraints,array([[0]])])
                gconstraints=concatenate_([gconstraints,I[[i],:]])
                
            if (x[i] - gradlag[i] >= ub[i]) and (abs(x[i] - ub[i]) < 1e-05):
                x_active[i]=1
                if options.verbose >= 3:
                    disp_('ub ',str(i),' is initially active')
                    
                # append active bound to the equality constraints

                constraints=concatenate_([constraints,array([[0]])])
                gconstraints=concatenate_([gconstraints,I[[i],:]])
                
        if options.verbose >= 3 and sum_(x_active) > 0:
            print("xactive="+str(x_active.T))
            
    # prepare g and M for slack handling
    
    if nbr_slacks:
        glocal = concatenate((glocal,zeros((nbr_slacks,1))))
        MM1 = concatenate((M,zeros((n,nbr_slacks))),axis=1)
        MM2 = concatenate((zeros((nbr_slacks,n)),0*eye(nbr_slacks)),axis=1)
        M = concatenate((MM1,MM2))
    
    #-------------------------------------------------
    # Start loop over active set
    #-------------------------------------------------
            
    finished=0
    iter_active=0
    
    while finished == 0:
    
        #----------------------------------------------------------
        # Restoration step computation 
        #----------------------------------------------------------

        iter_active=iter_active + 1
        
        if options.verbose >= 3:
            disp_('********* Loop over active set ********* iteration ',\
            str(iter_active),' *********')
        if options.verbose >= 3:
            fprintf_(options.fout,'  Restoration step:\n')
            
        delta_r=xi * delta

        g1=gconstraints.T.dot(constraints)
        #g1=numpy.dot(gconstraints.T,constraints)
        #if iter_active == 1 and info.feasn > 1e-7: 
        # no need to recompute the restoration step otherwise

        if iter_active == 1 and norm_(g1) > 1e-14: 
        # no need to recompute the restoration step otherwise
        
            # form the linear least squares system AT*A

            H1=gconstraints.T.dot(gconstraints)
            #H1=numpy.dot(gconstraints.T,gconstraints)
            # check condition of matrix H1
            
            check_condition=0

            if check_condition:
                cthreshold=1e+16
                H1,badcond=sqpdfo_check_cond_(H1,cthreshold,options,nargout=2)
                
            check_convex=1
            
            #  check that matrix H1 is convex (otherwise convexify)

            if check_convex:
                H1=sqpdfo_check_convex_(H1,options)
                
            # solve the LLS problem inside the variables bounds
            
            if nbr_slacks:
                lb_r=lb[0:n+nbr_slacks] - concatenate((x,sl))
                ub_r=ub[0:n+nbr_slacks] - concatenate((x,sl))
            else:
                lb_r=lb[0:n] - x
                ub_r=ub[0:n] - x
            
            lls_solver = 'bcdfo'

            if lls_solver == 'bcdfo':
                stratLam=1
                r_sl,_lambda,norms,value,gplus,nfact,neigd,msg=\
                bcdfo_solve_TR_MS_bc_(g1,H1,lb_r,ub_r,delta_r,prec_r,\
                                      stratLam,options,nargout=8)
                                      
            elif lls_solver == 'minq':
                cc = np.array(g1.T[0])
                G = np.matrix(H1)
                yu = np.array(np.maximum(lb_r, -delta_r).T[0])
                yo = np.array(np.minimum(ub_r, delta_r).T[0])
                prt = 0
                m = Minq(0, cc, G, yu, yo, prt)
                xy, value, msg, nsub = m.search()
                if msg > 0: msg = 'error'
                gplus = g1
                r_sl = numpy.array([xy]).T
                
            r = r_sl[0:n]

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
            
            if nbr_slacks:
                rpred_m=norm_(constraints) - norm_(constraints + gconstraints.dot(r_sl))
            else:   
                #rpred_m = norm_(info.ce)-norm_(info.ce+info.ae.dot(r)) 
                rpred_m=norm_(constraints) - norm_(constraints + gconstraints.dot(r))
            
            rpred=rpred + rpred_m
            
            if rpred < 0:
                if options.verbose >= 3:
                    fprintf_(options.fout,'\n### sqpdfo_solve_TR_bc:',\
                    ' rpred = %9.2e should not be negative\n\n'%(rpred))
                    
            active_r=(info_r.flag == 1) or (info_r.flag == 2)
            
            if options.verbose >= 5:
                if - 1 == info_r.flag:
                    fprintf_(options.fout,'    max of %0i iterations reached\n'\
                    %(20 * me))
                elif 0 == info_r.flag:
                    fprintf_(options.fout,'    precision is less than',\
                    ' required tolerance %8.2e\n'%(prec_r))
                elif 1 == info_r.flag:
                    fprintf_(options.fout,'    TR boundary is reached\n')
                elif 2 == info_r.flag:
                    fprintf_(options.fout,'    negative curvature direction',\
                    ' encountered\n')
                fprintf_(options.fout,'    |r|   = %8.2e\n'%(norm_r))
                fprintf_(options.fout,'    rpred = %8.2e\n'%(rpred))
                
        else:
            if options.verbose >= 3:
                fprintf_(options.fout,'    unchanged\n')
                
            r=zeros_(size_(x))
            if nbr_slacks:
                r_sl=zeros_like(concatenate((x,sl)))
            else:
                r_sl=zeros_like(x)
            active_r=copy(False)
            xr=copy(x)
            norm2_r=0.0
            norm_r=0.0
            
        if options.verbose == 3:
            disp_('r = (',str(r.T),')')
            disp_('delta_r = ',str(delta_r),', norm_r = ',str(norm_r))
            
        #-------------------------------------------------------------
        # Tangent step computation
        #-------------------------------------------------------------

        if options.verbose >= 3:
            fprintf_(options.fout,'  Tangent step:\n')
            
        delta_t=copy(delta)
        deg_freedom=n - length_(constraints) + nbr_slacks
        
        if deg_freedom > 0:
        
            # compute nullspace of the Jacobian and project M and g

            Z_=null_(gconstraints)
            M_t=Z_.T.dot(M.dot( Z_))
            g_t=Z_.T.dot((glocal + M.dot(r_sl)))
         
            #  Compute minimizer of the model in the tangent space

            u_sl,info_t=\
            sqpdfo_truncated_cg_(M_t,- g_t,delta_t,20 * (n-me+nbr_slacks),\
                               prec_t,nargout=2)
            
            t_sl=Z_.dot(u_sl)
            t = t_sl[0:n]
            
            active_t=(info_t.flag == 1) or (info_t.flag == 2)
            
            if options.verbose >= 5:
                if - 1 == info_t.flag:
                    fprintf_(options.fout,'    max of %0i iterations reached\n'\
                    %(20 * (n - me)))
                elif 0 == info_t.flag:
                    fprintf_(options.fout,'    precision is less than',\
                    ' required tolerance %8.2e\n'%(prec_t))
                elif 1 == info_t.flag:
                    fprintf_(options.fout,'    TR boundary is reached\n')
                elif 2 == info_t.flag:
                    fprintf_(options.fout,'    negative curvature direction',\
                    ' encountered\n')
                fprintf_(options.fout,'    |t| = %8.2e\n'%(norm_(t)))
                
        else:
            t_sl=zeros_(1,n+nbr_slacks).T
            t=zeros_(1,n).T
            active_t=0
            
        if options.verbose == 3:
            disp_('t = (',str(t_sl.T),')')
            disp_('delta_t = ',str(delta_t),', norm_t = ',str(norm_(t)))
            disp_('delta ',str(delta),', norm_s = ',str(norm_(r + t)))
            
        #-----------------------------------------------------
        # Compute trial values of x and slack variables
        #-----------------------------------------------------
        
        xnew=x + r + t
        
        if nbr_slacks:
            r_slsl = r_sl[n:n+nbr_slacks]
            t_slsl = t_sl[n:n+nbr_slacks]
            
            if iter_active == 1:
                slnew = sl + r_slsl + t_slsl
            else:
                slnew = slnew + t_slsl
                
            # check slacks to stay in their bounds
            
            slnew[slnew > ub[n:n+nbr_slacks]] = slnew[slnew > ub[n:n+nbr_slacks]] -\
                                                t_slsl[slnew > ub[n:n+nbr_slacks]]
            slnew[slnew < lb[n:n+nbr_slacks]] = slnew[slnew < lb[n:n+nbr_slacks]] -\
                                                t_slsl[slnew < lb[n:n+nbr_slacks]]
            
        #-----------------------------------------------------------------------
        # Compute violated bounds and append these to the equality constraints.
        #-----------------------------------------------------------------------

        x_active=zeros_(size_(xnew))
        x_viol=zeros_(size_(xnew))
        
        for i in range(0,n):
            if (xnew[i] - lb[i] < - threshold):
                x_viol[i]=-(i+1) 
                x_fix=concatenate_([x_fix, array([i])], axis=1)  
                
                if options.verbose >= 3:
                    disp_('lb ',str(i),' is violated')
                constraints=concatenate_([constraints,array([[0]])])
                gconstraints=concatenate_([gconstraints,I[[i],:]])
                violated=1
                break
                
            elif (abs(xnew[i] - lb[i]) < 1e-07):
                x_active[i]=-(i+1)
                if options.verbose >= 3:
                    disp_('lb ',str(i),' is active')
                    
            elif (xnew[i] - ub[i] > threshold):
                x_viol[i]=-(i+1)
                x_fix=concatenate_([x_fix, array([i])], axis=1) 
                if options.verbose >= 3:
                    disp_('ub ',str(i),' is violated')
                constraints=concatenate_([constraints,array([[0]])])
                gconstraints=concatenate_([gconstraints,I[[i],:]])
                violated=1
                break
                
            elif (abs(xnew[i] - ub[i]) < 1e-07):
                x_active[i]=-(i+1) 
                if options.verbose >= 3:
                    disp_('ub ',str(i),' is active')
                    
        if sum_(x_viol) == 0:
            violated=0
            if options.verbose >= 3:
                disp_('no new bound violated')
                
        else:
            if options.verbose >= 3:
                disp_(x_fix)
                
        #-----------------------------------------------------------------------
        # Check optimality of x (if step is zero)
        #-----------------------------------------------------------------------

        if norm_(r + t) <= 1e-16:
            if (iter_active >= 10 * n or delta < delta_min):
                if options.verbose >= 3:
                    disp_('### sqpdfo_solve_TR_bc: active-set iteration',\
                    ' limit exceeded ###')
                return xnew,delta,rpred,active_r,active_t,lm_computed,lm,info,slnew
                
            if len(x_fix) > 0:
            
              # compute new Lagrange multipliers

              lbounds=- inf * ones_(size_(x))
              ubounds=inf * ones_(size_(x))
              ilb=(abs(lb[indfree] - xnew) < 1e-05).reshape(-1)
              iub=(abs(ub[indfree] - xnew) < 1e-05).reshape(-1)
              lbounds[ilb]=lb[indfree[ilb]]
              ubounds[iub]=ub[indfree[iub]]
            
              lm,info=\
              sqpdfo_compute_multiplier_(xnew,lbounds,ubounds,info,options,values,\
                                       nargout=2)
            
              # compute smallest LM of the active bounds
                       
              min_lm,ind_min_lm=min_(lm[x_fix],nargout=2)

              if options.verbose >= 3:
                disp_('smallest Lagrange multiplier (for the bounds) = ',str(min_lm))
                
              if min_lm < 0:
                if options.verbose >= 3:
                    disp_('Zero step but not converged - release one bound!!')
                    
                constraints=constraints[0:me]
                gconstraints=gconstraints[0:me,:]
                
                if not isempty_(find_(x_active < 0)):
                    xa = find_(x_active < 0).reshape(-1)
                else:
                    xa = []
                
                if length_(xa) == length_(x_fix):
                    if (np.sort(xa) == np.sort(x_fix)).all():
                    
                        # check the Lagrange multipliers to release one bound

                        for i in range(0,length_(x_fix)):
                            if i != ind_min_lm:
                                constraints=concatenate_([constraints,array([[0]])])
                                gconstraints=concatenate_([gconstraints,I[[x_fix[i]],:]])
                        x_fix = np.delete(x_fix, ind_min_lm) 
                        finished=1
                    else:
                        # the active are not the fixed variables - just fix the
                        # currently active vars and release all other

                        x_fix=array([])
                        
                        for i in range(0,length_(find_(x_active < 0))):
                            constraints=concatenate_([constraints,array([[0]])])
                            gconstraints=concatenate_([gconstraints,I[[xa[i]],:]])
                            x_fix=concatenate_([x_fix,[xa[i]]],axis=1)
                else:
                    # number of active and fixed variables is not the same - 
                    # just fix the currently active vars and release all other

                    x_fix=array([])
                    
                    for i in range(0,length_(find_(x_active < 0))):
                        constraints=concatenate_([constraints,array([[0]])])
                        gconstraints=concatenate_([gconstraints,I[[xa[i]],:]])
                        x_fix=concatenate_([x_fix,[xa[i]]],axis=1)
                        
                if options.verbose >= 3:
                    print("x_fix="+str(x_fix))
                    
              else:
                if options.verbose >= 3:
                    disp_('Zero step and converged - go back to TR-loop...')
                finished=1
                
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
                    disp_('non zero infeasible step - continue finding',\
                    ' correct active set...')
                    
        #-----------------------------------------------------------------------
        # Update x, g, delta and ce (if violated bound)
        #-----------------------------------------------------------------------
   
        if violated == 1:
        
            tstep=copy(t_sl)
            
            if options.verbose >= 3:
                disp_('shorten tangential step')
                
            # compute ratio of current step inside the bounds
            
            if nbr_slacks:
                aT=concatenate_([eye_(n+nbr_slacks),-eye_(n+nbr_slacks)])
                aTx=aT.dot(concatenate_([xr,sl+r_sl[n:n+nbr_slacks]]))
                alpha=concatenate_([ub[0:n+nbr_slacks],-lb[0:n+nbr_slacks]])
            else:
                aT=concatenate_([eye_(n),- eye_(n)])
                aTx=aT.dot(xr)
                alpha=concatenate_([ub[0:n],- lb[0:n]])
                        
            divisor=aT.dot(tstep)
            indices1 = find_(divisor < 1e-15)
            indices2 = find_(divisor > -1e-15)
            iii = [i in indices2 for i in indices1]
            indices = indices1[iii]
            divisor[indices] = 1e-15
            ratio=(alpha - aTx) / divisor
            minratio=min_(ratio[divisor > 0])
            
            if (minratio < 0):
                # in case xk is slightly in the infeasible region, this should not
                # be corrected here (tstep points in ascending direction if minratio<0)
                minratio=0.
                
            else:
                minratio=float(minratio)
                
            # compute step inside the bounds

            tstep=minratio*tstep
            x=xr + tstep[0:n]
            step=r_sl + tstep
            
            if nbr_slacks:
                if iter_active == 1:
                    slnew = sl + r_sl[n:n+nbr_slacks] + tstep[n:n+nbr_slacks]
                else:
                    slnew = slnew + tstep[n:n+nbr_slacks]
            
            if options.verbose == 3:
                disp_('t = (',str(tstep[0:n].T),')')
                disp_('delta_t = ',str(delta_t),', norm_t = ',str(norm_(tstep[0:n])))
                disp_('delta ',str(delta),', norm_s = ',str(norm_(r + tstep[0:n])))
                
            # update gradient and Jacobian 

            glocal=glocal + M .dot(step)
            
            for i in range(0,me):
                gconstraints[i,:]=gconstraints[[i],:] + (M .dot( step)).T
                
            delta=delta - 0.5*norm_(step[0:n])

    return xnew,delta,rpred,active_r,active_t,lm_computed,lm,info,slnew

