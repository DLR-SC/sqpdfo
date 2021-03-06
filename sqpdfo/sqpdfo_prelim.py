# -*- coding: utf-8 -*-

import types
from sqpdfo.sqpdfo_options import *
from sqpdfo.bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from sqpdfo.bcdfo_poisedness_Y import bcdfo_poisedness_Y_
from sqpdfo.bcdfo_computeP import bcdfo_computeP_
from sqpdfo.bcdfo_gradP import bcdfo_gradP_
from sqpdfo.bcdfo_projgrad import bcdfo_projgrad_
from sqpdfo.bcdfo_repair_Y import bcdfo_repair_Y_
from sqpdfo.sqpdfo_augmX_evalf import sqpdfo_augmX_evalf_
from sqpdfo.sqpdfo_compute_multiplier import sqpdfo_compute_multiplier_
from sqpdfo.sqpdfo_optimality import sqpdfo_optimality_
import sqpdfo.sqpdfo_global_variables as glob
from copy import copy
from numpy import array, arange
from sqpdfo.helper import *


def sqpdfo_prelim_(func_=None,x0_=None,lm0_=None,Delta0_=None,lb_=None,ub_=None,\
   scaleX_=None,scalefacX_=None,cur_degree_=None,rep_degree_=None,plin_=None,\
   pdiag_=None,pquad_=None,c_=None,initial_Y_=None,kappa_ill_=None,\
   factor_FPR_=None,Lambda_FP_=None,Lambda_CP_=None,eps_L_=None,lSolver_=None,\
   hardcons_=None,stratLam_=None,xstatus_=None,sstatus_=None,dstatus_=None,\
   options_=None,*args,**kwargs):

###############################################################################
# This function realizes the following preliminary jobs:
# - build initial poised interpolation set
# - check the bounds and that initial point inside the bounds
# - check the given options
# - compute function and constraint values
# - compute initial multipliers (if not given)
# - initial printings
###############################################################################

    func=copy(func_)
    x0=copy(x0_)
    lm0=copy(lm0_)
    Delta0=copy(Delta0_)
    lb=copy(lb_)
    ub=copy(ub_)
    scaleX=copy(scaleX_)
    scalefacX=copy(scalefacX_)
    cur_degree=copy(cur_degree_)
    rep_degree=copy(rep_degree_)
    plin=copy(plin_)
    pdiag=copy(pdiag_)
    pquad=copy(pquad_)
    c=copy(c_)
    initial_Y=copy(initial_Y_)
    kappa_ill=copy(kappa_ill_)
    factor_FPR=copy(factor_FPR_)
    Lambda_FP=copy(Lambda_FP_)
    Lambda_CP=copy(Lambda_CP_)
    eps_L=copy(eps_L_)
    lSolver=copy(lSolver_)
    hardcons=copy(hardcons_)
    stratLam=copy(stratLam_)
    xstatus=copy(xstatus_)
    sstatus=copy(sstatus_)
    dstatus=copy(dstatus_)
    options=copy(options_)

    info = dummyUnionStruct()
    info.nsimul = array([])
    
    nbr_slacks = glob.get_nbr_slacks()
    sl = glob.get_slacks()
    
    Y = array([])
    gamma1 = 0.010000000000000
    eps =  2.220446049250313e-16
    stallfact=10 * eps
    
    nfix = None				
    indfix = None				
    xfix = None
    vstatus = None				
    QZ = None
    RZ = None
    scale = None
    poised = None
    Y_radius = None
    poised_model = None
    X = None
    fX = None
    #Y = None
    fY = None
    ciX = None
    ciY = None
    ceX = None
    ceY = None
    poisedness_known = None
    m = None
    normgx = None
    fcmodel = None
    ind_Y = None
    i_xbest = None
    indfree = None
    
    # Set output arguments
    
    n=0
    nb=0
    mi=0
    me=0
    lm=array([])
    info.f=np.NaN
    info.ce=array([])
    info.g=[]
    info.ai=[]
    info.ae=[]
    info.hl=[]
    info.niter=0
    info.glagn = float('NaN')
    info.feasn = float('NaN')
    info.compl = float('NaN')
    shift_Y=1
    x=copy(np.NaN)
    fx=copy(np.NaN)
    gx=copy(np.NaN)

    # check user options or get default values if no values given by the user

    info,options,values=sqpdfo_options_(info,options,nargout=3)
    if info.flag:
        return n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta0,nfix,indfix,xfix,vstatus,xstatus,\
           sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,poised_model,X,fX,Y,fY,\
           ciX,ciY,ceX,ceY,poisedness_known,m,gx,normgx,fcmodel,ind_Y,i_xbest,\
           cur_degree,rep_degree,plin,pdiag,pquad,indfree,info,options,values
    info.nsimul=np.zeros(values.nsimultype)
    whichmodel = options.whichmodel

    # Check the argument x0; deduce n

    n=size_(x0,1)
    if size_(x0,2) != 1:
        if options.verbose:
            fprintf_(options.fout,'### sqpdfo_prelim: the initial x must be an n-vector\n\n')
        info.flag=values.fail_on_argument
        return n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta0,nfix,indfix,xfix,vstatus,xstatus,\
           sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,poised_model,X,fX,Y,fY,\
           ciX,ciY,ceX,ceY,poisedness_known,m,gx,normgx,fcmodel,ind_Y,i_xbest,\
           cur_degree,rep_degree,plin,pdiag,pquad,indfree,info,options,values
    if n < 1:
        if options.verbose:
            fprintf_(options.fout,'### sqpdfo_prelim: the initial x must be an n-vector with n > 0\n\n')
        info.flag=values.fail_on_argument
        return n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta0,nfix,indfix,xfix,vstatus,xstatus,\
           sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,poised_model,X,fX,Y,fY,\
           ciX,ciY,ceX,ceY,poisedness_known,m,gx,normgx,fcmodel,ind_Y,i_xbest,\
           cur_degree,rep_degree,plin,pdiag,pquad,indfree,info,options,values
           
    # Compute the number of bounds

    if isempty_(lb):
        lb=- options.inf * ones_(n + mi,1)
    if isempty_(ub):
        ub=options.inf * ones_(n + mi,1)
    nb_lo=sum_(lb[0:n] > - options.inf)
    nb_up=sum_(ub[0:n] < options.inf)
    nb=sum_(min_((lb[0:n] > - options.inf) + (ub[0:n] < options.inf),1))
    
    #  Checking the bounds and correct Delta0 if there is insufficient space 
    #  between the bounds. Modification of x0 if construction of first
    #  interpolation model would interfere with given bounds.
    
    zero=0.0
    nfix=0
    indfix=array([])
    xfix=zeros_(n,1)
    vstatus=zeros_(n,1)
    temp=zeros_(n,1)
    
    for j in range(0,n):
    
        #  Check lower and upper bounds.

        if (lb[j] > ub[j]):
            disp_('Error: Lower bound of component ',str(j),' exceeds upper bound !!')
            info.flag=2
            return n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta0,nfix,indfix,xfix,vstatus,\
            xstatus,sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,poised_model,X,fX,Y,fY,\
            ciX,ciY,ceX,ceY,poisedness_known,m,gx,normgx,fcmodel,ind_Y,i_xbest,\
            cur_degree,rep_degree,plin,pdiag,pquad,indfree,info,options,values
            
        #  Check difference between bounds.

        temp[j]=ub[j] - lb[j]
        if (temp[j] < Delta0 + Delta0):
            if (temp[j] == zero):
                nfix=nfix + 1
                indfix=concatenate_([indfix,j],axis=1)
                vstatus[j]=c.alwaysfixed
                xfix[j]=lb[j]
                continue
            else:
                Delta0=0.5 * temp[j]
                disp_(' Diff. between lower and upper bound of component ',str(j),\
                ' is less than 2*Delta0 !! New Delta0=',str(Delta0))
                
        #  Move the starting point inside the bounds if necessary

        templ=lb[j] - x0[j]
        tempu=ub[j] - x0[j]
        if (templ >= - Delta0):
            x0[j]=lb[j] + Delta0
        else:
            if (tempu <= Delta0):
                x0[j]=ub[j] - Delta0
                
    #  Scale x0 and bounds if user-defined.

    if (scaleX):
        for i in range(0,n):
            if (scalefacX[i] > 0):
                x0[i]=x0[i] * scalefacX[i]
                lb[i]=lb[i] * scalefacX[i]
                ub[i]=ub[i] * scalefacX[i]
            else:
                scalefacX[i]=1
                
    #  Reset constants if fixed some variables.

    if (nfix > 0):
        nfree=n - nfix
        if (nfree <= 0):
            disp_('No free variables. Please, enlarge search space!')
            info.flag=2
            return n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta0,nfix,indfix,xfix,vstatus,\
            xstatus,sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,poised_model,X,fX,Y,fY,\
            ciX,ciY,ceX,ceY,poisedness_known,m,gx,normgx,fcmodel,ind_Y,i_xbest,\
            cur_degree,rep_degree,plin,pdiag,pquad,indfree,info,options,values
            
        indfree=setdiff_(arange(0,n),indfix)
        x0=x0[indfree]
        
        if (cur_degree == plin):
            cur_degree=nfree + 1
        else:
            if (cur_degree == pdiag):
                cur_degree=2 * nfree + 1
            else:
                if (cur_degree == pquad):
                    cur_degree=((nfree + 1) * (nfree + 2)) / 2
        if (rep_degree == plin):
            rep_degree=nfree + 1
        else:
            if (rep_degree == pdiag):
                rep_degree=2 * nfree + 1
            else:
                if (rep_degree == pquad):
                    rep_degree=((nfree + 1) * (nfree + 2)) / 2
                    
        plin=nfree + 1
        pdiag=2 * nfree + 1
        pquad=((nfree + 1) * (nfree + 2)) / 2
        n=copy(nfree)
        
    else:
        indfree=arange(0,n)
        
    x=copy(x0)
    
    # Compute an interpolation set around x0 and compute f, ci, ce, g, ai, ae

    getfY=1
    while (getfY):

        if (options.verbose > 2):
            disp_(' Degree of the initial  model = ',str(cur_degree))
            
        #  Compute an initial poised interpolation set around the starting point.
   
        if initial_Y=='random':
            Y[:,0]=x0
            
            #  Loop in case of an accidentally ill-conditioned initial system

            ill_init=1
            while (ill_init):
                Y[:,1:cur_degree]=-ones_(n,cur_degree-1) + 2 * rand_(n,cur_degree-1) 
                for j in range(1,cur_degree): 
                    Y[:,j]=Y[:,0] + Y[:,j]*(Delta0 / norm_(Y[:,j]))
                    
                QZ,RZ,x,scale=\
                bcdfo_build_QR_of_Y_(Y,whichmodel,shift_Y,Delta0,1,kappa_ill,nargout=4)
                
                # check the condition

                if (cond_(RZ) < kappa_ill):
                    ill_init=0

            #  Make the set poised.

            QZ,RZ,Y,replaced,poised,Y_radius,x,scale=\
            bcdfo_repair_Y_(QZ,RZ,Y,Delta0,factor_FPR,Lambda_FP,Lambda_CP,eps_L,x,\
                lSolver,whichmodel,hardcons,lb,ub,indfree,stratLam,scale,shift_Y,\
                1,kappa_ill,nargout=8)
               
            poisedness_known=1
            
        elif initial_Y=='simplx':

            #  Compute the initial interpolation set (simplex plus midpoints).

            I=eye_(n)
            Y=zeros_(n,n+1)
            Y[:,0]=x0.reshape(-1)

            for j in range(0,n):
                # initial degree is linear
                step1=- Delta0
                Y[:,j + 1]=x0.reshape(-1) + step1 * I[:,j]

                if (cur_degree >= pdiag):
                    # initial degree is diagonal
                    step2=copy(Delta0)
                    Y[:,j + 1 + n]=x0.reshape(-1) + step2 * I[:,j]

        #    if (cur_degree == pquad):
                 # initial degree is quadratic
        #        k=2 * n + 2
        #        for j in range(0,n):
        #            for jj in range(j + 1,n+1):
        #                Y[:,k-1]=0.5 * (Y[:,j + 1] + Y[:,jj + 1])
        #                k=k + 1

            #  Build the initial factorization.

            QZ,RZ,x,scale=\
            bcdfo_build_QR_of_Y_(Y,whichmodel,shift_Y,Delta0,1,kappa_ill,nargout=4)
            
            poised,Y_radius=\
            bcdfo_poisedness_Y_(QZ,RZ,Y,eps_L,x,1,whichmodel,hardcons,lb,ub,indfree,\
                stratLam,scale,shift_Y,nargout=2)
                
            poisedness_known=1
        
        poised_model=1 # The initial interpolation set is poised.

       #  Compute the associated function values, possibly reducing Delta0 to
       #  ensure that the objective function remains finite at all interpolation
       #  points.
       
        X=array([])
        fX=array([])
        ciX=array([])
        ceX=array([])
        ind_Y=array([])
        info.ci=array([])
        info.ce=array([])
        
        for i in range(0,cur_degree):
            X,fX,ciX,ceX,neval,xstatus,sstatus,dstatus,info,outdic=\
            sqpdfo_augmX_evalf_(func,Y[:,[i]],i,X,fX,ciX,ceX,nfix,xfix,indfix,\
                indfree,1e+25,info.nsimul[1],xstatus,c.inY,sstatus,\
                dstatus,scaleX,scalefacX,info,options,values,nargout=9)
                
            if info.flag:
                return n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta0,nfix,indfix,xfix,\
                   vstatus,xstatus,sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,\
                   poised_model,X,fX,Y,fY,ciX,ciY,ceX,ceY,poisedness_known,m,\
                   gx,normgx,fcmodel,ind_Y,i_xbest,cur_degree,rep_degree,plin,pdiag,\
                   pquad,indfree,info,options,values

            #  If the computed function value is infinite, restart with a smaller Delta0.

            if (abs(fX[i]) > 1e+25):
                break
                
            #  All functions values at points of Y are finite.  No need for
            #  another pass.
            
            if (i == cur_degree-1):
                getfY=0
            ind_Y=concatenate_([ind_Y,array([i])],axis=1)
            
        # Check for infinity constraint values
        
        ceX[ceX>=1e25] = 10 * max_(ceX[ceX<1e25])
        ceX[ceX<=-1e25] = 10 * min_(ceX[ceX>-1e25])
                    
        fY=copy(fX)
        ciY=copy(ciX)
        ceY=copy(ceX)
        xstatus=xstatus.T
        dstatus=dstatus.T

        #  Another pass is needed with a smaller Delta0 (at least one function
        #  value is infinite). 
        
        if (getfY):
        
            #  Decrease Delta0.

            Delta0=gamma1 * Delta0

            #  Terminate with an error message if Delta0 becomes negligible wrt x0.

            if (Delta0 < stallfact * norm_(x0)):
                disp_('Error: cannot find enough finite objective function values',\
                ' in the neighbourhood of the starting point! Terminating.')

                #  including fixed variables at return

                if (nfix > 0):
                    I=eye_(n + nfix)
                    x=I[:,indfix] * zeros_(nfix,1) + I[:,indfree].dot( x )
                    gx=I[:,indfix] * zeros_(nfix,1) + I[:,indfree].dot( gx )
                    
                return n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta0,nfix,indfix,xfix,\
                vstatus,xstatus,sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,\
                poised_model,X,fX,Y,fY,ciX,ciY,ceX,ceY,poisedness_known,m,\
                gx,normgx,fcmodel,ind_Y,i_xbest,cur_degree,rep_degree,plin,\
                pdiag,pquad,indfree,info,options,values

    fx0=copy(fY[0])
    info.f=fx0
    m=copy(cur_degree)-1
    i_xbest=0
    
    #  Move to the best point in the interpolation set, if different from x0.
    #  ATTENTION for constrained problems...! Check merit function value!
    
    #[x, fx, QZ, RZ, Y, fY, ciY, ceY, ind_Y, i_xbest, scale] = ...
    #   sqpdfo_find_smallf(c, QZ, RZ, Y, fY, ciY, ceY,...
    #   1:cur_degree, 1, cur_degree, indfree, x, lb, ub, fx0, ...
    #   zeros(cur_degree,1), whichmodel, scale, shift_Y, Delta0, 1, kappa_ill);
    
    #  Compute the associated polynomial interpolation model(s) 
    #  for objective and possible constraints.
    
    initmodel=zeros_(1,pquad)
    rhsY=concatenate_([fY.reshape(1,-1),ciY,ceY], axis=0)
    
    fcmodel=\
    bcdfo_computeP_(QZ,RZ,Y,rhsY,whichmodel,initmodel,ind_Y,0,0,gx,scale,shift_Y,Delta0)
    gx=bcdfo_gradP_(fcmodel[[0],:],x,x,scale,shift_Y)
    normgx,_=bcdfo_projgrad_(n,x,gx,lb[indfree],ub[indfree])
    
    if any_(size_(gx) != [n,1]):
        if options.verbose:
            fprintf_(options.fout,'### sqpdfo: the computed gradient g has a wrong ',\
            'size, (%0i,%0i) instead of (%0i,1)\n\n'%(size_(gx),n))
            
        info.flag=values.fail_on_simul
        
        return n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta0,nfix,indfix,xfix,vstatus,\
        xstatus,sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,poised_model,X,fX,\
        Y,fY,ciX,ciY,ceX,ceY,poisedness_known,m,gx,normgx,fcmodel,ind_Y,\
        i_xbest,cur_degree,rep_degree,plin,pdiag,pquad,indfree,info,options,values
        
    info.g=gx

    #  Check constraints and deduce mi and me

    mi=size_(ciY,1)
    if mi > 0:
        info.ci=copy(ciY[:,[0]])

        # compute model gradient associated with the inequality constraints

        gci=zeros_(mi,n)
        for i in range(0,mi):
            gci[i,:]=bcdfo_gradP_(fcmodel[[1 + i],:],x,x,scale,shift_Y).T
        info.ai=copy(gci)
    else:
        info.ci=array([])
        info.ai=array([])
        
    me=size_(ceY,1)
    if me > 0:
        info.ce=copy(ceY[:,[0]])

        # compute model gradient associated with the equality constraints

        gce=zeros_(me,n)
        for i in range(0,me):
            gce[i,:]=bcdfo_gradP_(fcmodel[[1 + mi + i],:],x,x,scale,shift_Y).T
        info.ae=gce
    else:
        info.ce=array([])
        info.ae=array([])
        
    # Initialize slack variables if necessary
    
    if nbr_slacks > 0:
    
        ce = info.ce

        for i in range(0,nbr_slacks):
            if ce[-nbr_slacks+i] > 0:
                sl[i] = min_(sqrt_(ce[-nbr_slacks+i]), sqrt_(ub[-nbr_slacks+i]))
            else:
                sl[i] = 0
                
        glob.set_slacks(sl)
        
    # Initial printing (1)

    fprintf_('\n')
    fprintf_('**************************************************************************************\n')
    fprintf_('*                                                                                    *\n')
    fprintf_('*       SQPDFO: Sequential-Quadratic-Programming Derivative-Free Optimization        *\n')
    fprintf_('*                                                                                    *\n')
    fprintf_('*                        (c)  A. Troeltzsch, 2013-2019                               *\n')
    fprintf_('*                                                                                    *\n')
    fprintf_('**************************************************************************************\n')
    fprintf_('\n')
    values.dline='-------------------------------------------'
    values.dline=strcat_(values.dline,values.dline) # 80 dashes
    values.eline='==========================================='
    values.eline=strcat_(values.eline,values.eline)  # 80 '=' characters
    values.sline='*******************************************'
    values.sline=strcat_(values.sline,values.sline) # 80 '*' characters
    if options.verbose > 0 and options.verbose < 4:
        fprintf_(options.fout,' iter  neval     fvalue            merit      ')
        if (me > 0):
            fprintf_(options.fout,' |grad Lag|  feasibility')
        else:
            fprintf_(options.fout,'  gradient  ')
        fprintf_(options.fout,'  delta     stepsize')
        if options.hess_approx == values.bfgs:
            fprintf_(options.fout,'  BFGS\n')
        else:
            fprintf_(options.fout,'  \n')
        fprintf_(options.fout,'  \n')
    if options.verbose >= 4:
        fprintf_(options.fout,'\nSQPDFO optimization solver\n\n')
        if isinstance(func, types.FunctionType):
            func_name=str(func)
        else:
            func_name=copy(func)
        fprintf_(options.fout,'  function name: "%s"\n'%(func_name))
        fprintf_(options.fout,'  dimensions:\n')
        fprintf_(options.fout,'  . variables (n):               %4i\n'%(n))
        if nb > 0:
            fprintf_(options.fout,'  . bounds on variables (nb):    (%0i lower, %0i double, %0i upper)\n'%(nb,nb_lo,nb_up))
        if mi > 0:
            fprintf_(options.fout,'  . inequality constraints (mi): %4i\n'%(mi))
        if me > 0:
            fprintf_(options.fout,'  . equality constraints (me):   %4i\n'%(me))
        fprintf_(options.fout,'  required tolerances for optimality:\n')
        if nb + mi + me > 0:
            fprintf_(options.fout,'  . gradient of the Lagrangian      %8.2e\n'%(options.tol_grad))
            fprintf_(options.fout,'  . feasibility                     %8.2e\n'%(options.tol_feas))
        else:
            fprintf_(options.fout,'  . gradient of the cost function   %8.2e\n'%(options.tol[0]))
        fprintf_(options.fout,'  counters:\n')
        fprintf_(options.fout,'  . max iterations                  %4i\n'%(options.miter))
        fprintf_(options.fout,'  . max function evaluations        %4i\n'%(options.msimul))
        fprintf_(options.fout,'  algorithm:\n')
        fprintf_(options.fout,'  . globalization by trust regions\n')
        fprintf_(options.fout,'  various input/initial values:\n')
        if  (nb + mi + me == 0) and (options.df1 > 0) and (info.f > 0):
            fprintf_(options.fout,'  . expected initial decrease       %8.2e\n'%(options.df1 * info.f))
        if nb + mi > 0:
            fprintf_(options.fout,'  . infinite bound threshold        %8.2e\n'%(options.inf))
        fprintf_(options.fout,'  . |x|_2                           %8.2e\n'%(norm_(x)))

    # Compute an initial multiplier if not given (takes a while); info.ci must be known

    if (nb + mi + me > 0):
        if isempty_(lm0):
            lm,info=sqpdfo_compute_multiplier_(x,[],[],info,options,values,nargout=2)
            
            if info.flag:
                return n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta0,nfix,indfix,xfix,vstatus,\
                xstatus,sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,poised_model,X,fX,\
                Y,fY,ciX,ciY,ceX,ceY,poisedness_known,m,gx,normgx,fcmodel,ind_Y,\
                i_xbest,cur_degree,rep_degree,plin,pdiag,pquad,indfree,\
                info,options,values
                
            if options.verbose >= 4:
                fprintf_(options.fout,'  . |lm|_2                          %8.2e (default: least-squares value)\n'%(norm_(lm)))
                
        else:
            lm=copy(lm0)
            if options.verbose >= 4:
                fprintf_(options.fout,'  . |lm|_2                          %8.2e\n'%(norm_(lm)))

    # Initial optimality

    feas,compl,info=\
    sqpdfo_optimality_(x,lm,lb[indfree],ub[indfree],info,options,nargout=3)
    
    if info.flag:
        return n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta0,nfix,indfix,xfix,vstatus,\
        xstatus,sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,poised_model,\
        X,fX,Y,fY,ciX,ciY,ceX,ceY,poisedness_known,m,gx,normgx,fcmodel,\
        ind_Y,i_xbest,cur_degree,rep_degree,plin,pdiag,pquad,indfree,info,options,values

    # Initial printing (2)

    if options.verbose >= 4:
        fprintf_(options.fout,'  . |g|_inf                         %8.2e\n'\
        %(norm_(info.g,inf)))
        if nb + mi + me > 0:
            fprintf_(options.fout,'  . |glag|_inf                      %8.2e\n'\
            %(norm_(info.glag,inf)))
        if nb:
            fprintf_(options.fout,'  . |x^#|_inf                       %8.2e\n'\
            %(norm_(feas[0:n],inf)))
        if mi:
            fprintf_(options.fout,'  . |ci^#|_inf                      %8.2e\n'\
            %(norm_(feas[n:n + mi],inf)))
        if me:
            fprintf_(options.fout,'  . |ce|_inf                        %8.2e\n'\
            %(norm_(feas[n + mi:n + mi + me],inf)))
        fprintf_(options.fout,'  tunings:\n')
        fprintf_(options.fout,'  . printing level                  %0i\n'\
        %(options.verbose))

    if info.flag:
        return n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta0,nfix,indfix,xfix,vstatus,\
        xstatus,sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,poised_model,X,fX,\
        Y,fY,ciX,ciY,ceX,ceY,poisedness_known,m,gx,normgx,fcmodel,ind_Y,\
        i_xbest,cur_degree,rep_degree,plin,pdiag,pquad,indfree,info,options,values
        
    return n,nb,mi,me,x,lm,lb,ub,scalefacX,Delta0,nfix,indfix,xfix,vstatus,\
    xstatus,sstatus,dstatus,QZ,RZ,scale,poised,Y_radius,poised_model,X,fX,\
    Y,fY,ciX,ciY,ceX,ceY,poisedness_known,m,gx,normgx,fcmodel,ind_Y,\
    i_xbest,cur_degree,rep_degree,plin,pdiag,pquad,indfree,info,options,values
