# -*- coding: utf-8 -*-
from __future__ import division
#try:
from runtime import *
from numpy import arange,append
#except ImportError:
    #from smop.runtime import *
def blls_(A=None,b=None,lb=None,ub=None,*args,**kwargs):
    """
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  blls is a solver for bound-constrained linear least-squares problems, that
%  is problem where, given 
%  - a m x n matrix A,
%  - a m x 1 right-hand side vector b,
%  - a n x 1 vector of (possibly -Inf) lower bounds lb,
%  - a n x 1 vector of (possibly +Inf) upper bounds ub,
%  one seeks a solution vector s solving
%
%           min || As - b ||^2     subject to    lb <= s <= ub
%
%  where ||.|| is the Euclidean norm and the contraint is understood 
%  componentwise.  No restriction is made on the dimension of A or its
%  rank, except that n>0 and m>0.
%
%  On exit, s contains the problem's solution if exitc = 0. If exitc = 1, 
%  s contains the best approximate solution found in maxiter = 2*n iterations.
%  The scalar resn contains the associated residual's norm, ie. ||As-b||.
%  If exitc = -1, the algorithm failed in the sense that the Cauchy step (see
%  below) could not be computed in maxbck backtackings (this error sould never
%  happen).
%
%  The maximum number of iterations, as other algorithmic parameters
%  (including tolerances for primal and dual feasibility and verbosity level),
%  may be modified in the "Algorithmic parameters" section in the beginning of
%  the code.
%
%  The method is intended for small-dimensional problems.  It is an active-set
%  algorithm where the unconstrained problem is solved at each iteration in
%  the subspace defined by the currently active bounds, themselves being
%  determined by a projected Cauchy step. Each subspace solution is computed
%  using a SVD decomposition of the reduced matrix.
%
%  Programming : Ph. Sampaio, Ph. Toint, A. Troeltzsch, April 2014.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    """
#    varargin = cellarray(args)
#    nargin = 4-[A,b,lb,ub].count(None)+len(args)

    m,n=size_(A,nargout=2)
    verbose=0
    epsfeas=1e-12
    epsconv=1e-10
    epsres=1e-11
    maxiter=2 * n
    epsdzer=1e-14
    armijob=0.5
    armijor=0.01
    maxback=15
    inds=arange(0,n)
    nit=0
    nuns=0
    ssub=zeros_(n,1)
    exitc=0
    s=min_(max_(pinv_(A).dot(b),lb),ub)
    res=A .dot(s) - b
    resn=norm_(res)
    g=A.T.dot(res)
    stry=min_(max_(s - g,lb),ub) - s
    opt=norm_(stry)
    free=find_(stry)
    atlb=find_(abs(max_(s - g,lb) - s) <= epsfeas)
    atub=find_(abs(min_(s - g,ub) - s) <= epsfeas)
    latlb=length_(atlb)
    latub=length_(atub)
    lfree=length_(free)
    if (verbose > 0):
        disp_(' ')
        disp_('   **************************************************************')
        disp_('   *                                                            *')
        disp_('   *                          BLLS                              *')
        disp_('   *                                                            *')
        disp_('   *   a direct bound-constrained linear least-squares solver   *')
        disp_('   *                                                            *')
        disp_('   *                                                            *')
        disp_('   *     (c) Ph. Sampaio, Ph. L. Toint, A. Troeltzsch, 2014     *')
        disp_('   *                                                            *')
        disp_('   **************************************************************')
        disp_(' ')
        disp_('     The problem has ',int2str_(n),' variables and ',int2str_(m),' rows.')
        disp_(' ')
        if (verbose > 2):
            problem_matrix=copy(A)
            right_hand_side=b.T
            lower_bounds=lb.T
            upper_bounds=ub.T
            disp_(' ')
        fprintf_('     nit     ||r||    optimality')
        fprintf_('                 nfr nlow nupp\n\n')
        fprintf_('   %5d  %.4e  %.4e                %4d %4d %4d\n' % (nit,resn,opt,lfree,latlb,latub))
        if (verbose > 1):
            if (verbose > 2):
                unconstrained_solution=s.T
            disp_(' ')
            disp_('   --------------------------------------------------------------')
            disp_(' ')
    if (opt <= epsconv or resn <= epsres):
        maxit=0
    else:
        maxit=copy(maxiter)
    for i in arange(1,maxit+1):
        nit=nit + 1
        if (verbose > 1):
            disp_('   Iteration ',int2str_(i))
            disp_(' ')
            fprintf_('   Cauchy point projected search\n')
            fprintf_('       k     ||r||     stepsize')
            fprintf_('                  nfr nlow nupp\n')
            fprintf_('   %5d  %.4e                            %4d %4d %4d'%(0,resn,lfree,latlb,latub))
        g=A.T.dot(res)
        alpha=(norm_(g) / norm_(A.dot( g))) ** 2
        for kc in arange(1,maxback+1):
            stry=min_(max_(s - alpha*g,lb),ub)
            dtry=stry - s
            ltry=g.T.dot(dtry)
            qtry=ltry + 0.5 * norm_(A .dot(dtry)) ** 2
            if (verbose > 1):
                fprintf_('\n   %5d  %.4e  %.4e'&(kc,norm_(A * stry - b),alpha))
            if (qtry <= armijor * ltry):
                break
            else:
                alpha=armijob * alpha
        if (kc >= maxback):
            exitc=- 1
            break
        atlb=inds[find_(abs(stry - lb) <= epsfeas)].reshape(-1)
        atub=inds[find_(abs(stry - ub) <= epsfeas)].reshape(-1)
        atb=append(atlb,atub)
        latlb=length_(atlb)
        latub=length_(atub)
        free=copy(inds)
        free = np.delete(free, atb)
        lfree=length_(free)
        s[atlb]=lb[atlb]
        s[atub]=ub[atub]
        s[free]=stry[free]
        res=A.dot(s) - b
        resn=norm_(res)
        if (verbose > 1):
            fprintf_('                %4d %4d %4d\n'%(lfree,latlb,latub))
            if (verbose > 2):
                Cauchy_point=s.T
                indices_of_free_variables=copy(free)
                indices_of_variables_at_their_lower_bound=copy(atlb)
                indices_of_variables_at_their_upper_bound=copy(atub)
        if (lfree == 0 or resn <= epsres):
            if (verbose > 1):
                fprintf_('   No nested subspace search\n')
            opt=0
        else:
            if (verbose > 1):
                fprintf_('   Nested subspace search\n')
                fprintf_('       k     ||r||     stepsize      ||r*||')
                fprintf_('      nfr nlow nupp\n')
                fprintf_('   %5d  %.4e                            %4d %4d %4d\n'&(0,resn,lfree,latlb,latub))
            for k in arange(0,n):
                if (verbose > 2):
                    disp_('    > Solving in subspace ',int2str_(k))
                    indices_of_free_variables=copy(free)
                    indices_of_variables_at_their_lower_bound=copy(atlb)
                    indices_of_variables_at_their_upper_bound=copy(atub)
                rhs=copy(b)
                if not isempty_(atlb):
                    rhs=rhs - A[0:m,atlb].dot(lb[atlb])
                if not isempty_(atub):
                    rhs=rhs - A[0:m,atub].dot(ub[atub])
                ssub[free]=pinv_(A[0:m,free]).dot(rhs)
                ssub[atlb]=lb[atlb]
                ssub[atub]=ub[atub]
                rsubo=A .dot(ssub) - b
                rsubon=norm_(rsubo)
                natlb=find_(ssub[free] < lb[free])
                natub=find_(ssub[free] > ub[free])
                lnatb=length_(natlb) + length_(natub)
                if (lnatb > 0):
                    alpha=1
                    dtry=ssub - s
                    rred=rsubon - resn
                    found=0
                    nback=4 * (1 - nuns)
                    for kb in arange(1,nback+1):
                        stry=(1 - alpha) * s+ alpha * ssub
                        natlbt=free[find_(stry[free] < lb[free])]
                        natubt=free[find_(stry[free] > ub[free])]
                        lnatbt=length_(natlbt) + length_(natubt)
                        stry=min_(max_(stry,lb),ub)
                        if (verbose >= 1):
                            rtry=A.dot(stry) - b
                            rtryn=norm_(rtry)
                            atlb=append(atlb,natlbt)
                            atub=append(atub,natubt)
                            atb=append(atlbt,atubt)
                            freet=copy(inds)
#                            freet[atbt]=[]
                            freet = np.delete(free, atbt)
                            latlbt=length_(atlbt)
                            latubt=length_(atubt)
                            lfreet=length_(freet)
                            fprintf_('   %5dp %.4e  %.4e   %.4e   %4d %4d %4d\n'%(kb,rtryn,alpha,rsubon,lfreet,latlbt,latubt))
                        if (lnatbt == 0):
                            break
                        if (verbose == 0):
                            rtry=A.dot(stry) - b
                            rtryn=norm_(rtry)
                        if (rtryn <= resn - armijor * alpha * rred):
                            s=copy(stry)
                            res=copy(rtry)
                            resn=copy(rtryn)
                            if (verbose == 0):
                                atlb=append(atlb,natlbt)
                                atub=append(atub,natubt)
                                atb=append(atlb,atub)
                                free=copy(inds)
#                                free[atb]=[]
                                free = np.delete(free, atb)
                                latlb=length_(atlb)
                                latub=length_(atub)
                                lfree=length_(free)
                            else:
                                atlb=copy(atlbt)
                                atub=copy(atubt)
                                free=copy(freet)
                                latlb=copy(latlbt)
                                latub=copy(latubt)
                                lfree=copy(lfreet)
                            found=1
                            break
                        alpha=armijob * alpha
                    if (found):
                        break
                    else:
                        if (kb >= nback):
                            nuns=nuns + 1
                        alpha=1
                        for kf in arange(0,length_(free)):
                            kk=free[kf]
                            if (dtry[kk] >= epsdzer):
                                alpha=min_(alpha,(ub[kk] - s[kk]) / dtry[kk])
                            else:
                                if (dtry[kk] <= - epsdzer):
                                    alpha=min_(alpha,(lb[kk] - s[kk]) / dtry[kk])
                        ssub=s + alpha * dtry
                        rsub=(1 - alpha) * res + alpha * rsubo
                        rsubn=norm_(rsub)
                        if (verbose > 1):
                            fprintf_('   %5ds %.4e  %.4e   %.4e   %4d %4d %4d\n'%(k,rsubn,alpha,rsubon,lfree,latlb,latub))
                        natlb=free[find_(abs(ssub[free] - lb[free]) <= epsfeas)].reshape(-1)
                        natub=free[find_(abs(ssub[free] - ub[free]) <= epsfeas)].reshape(-1)
                        atlb=append(atlb,natlbt)
                        atub=append(atub,natubt)
                        atb=append(atlb,atub)
                        free=copy(inds)
#                        free[atb]=[]
                        free = np.delete(free, atb)
                        latlb=length_(atlb)
                        latub=length_(atub)
                        lfree=length_(free)
                        if (verbose > 2):
                            current_subspace_solution=ssub.T
                            indices_of_free_variables=copy(free)
                            indices_of_variables_at_their_lower_bound=copy(atlb)
                            indices_of_variables_at_their_upper_bound=copy(atub)
                        s=copy(ssub)
                        res=copy(rsub)
                        resn=copy(rsubn)
                else:
                    s=copy(ssub)
                    res=copy(rsubo)
                    resn=copy(rsubon)
                    if (verbose > 1):
                        fprintf_('   %5df %.4e  %.4e   %.4e   %4d %4d %4d\n'%(k,resn,1,resn,lfree,latlb,latub))
                        if (verbose > 2):
                            current_subspace_solution=ssub.T
                            indices_of_variables_at_their_lower_bound=copy(atlb)
                            indices_of_variables_at_their_upper_bound=copy(atub)
                            free
                    break
            opt=norm_(min_(max_(s - A.T.dot(res),lb),ub) - s)
        if (verbose == 1):
            fprintf_('   %5d  %.4e  %.4e                %4d %4d %4d\n'%(nit,resn,opt,lfree,latlb,latub))
        else:
            if (verbose > 1):
                disp_(' ')
                fprintf_('     nit    ||r||     optimality')
                fprintf_('                 nfr nlow nupp\n\n')
                fprintf_('   %5d  %.4e  %.4e                %4d %4d %4d\n'%(nit,resn,opt,lfree,latlb,latub))
                if (verbose > 2):
                    current_solution=s.T
                    indices_of_free_variables=copy(free)
                    indices_of_variables_at_their_lower_bound=copy(atlb)
                    indices_of_variables_at_their_upper_bound=copy(atub)
                disp_('   --------------------------------------------------------------')
                disp_(' ')
        if (opt <= epsconv or resn <= epsres):
            break
    if (exitc == 0 and nit >= maxiter):
        exitc=1
    if (verbose > 0):
        disp_(' ')
        if (verbose > 2):
            indices_of_free_variables=copy(free)
            indices_of_variables_at_their_lower_bound=copy(atlb)
            indices_of_variables_at_their_upper_bound=copy(atub)
            final_solution=s.T
            final_residual=res.T
        if (exitc == 1):
            disp_('   !!! maxit reached !!!')
            keyboard
        else:
            if (exitc == - 1):
                disp_('   !!! Cauchy point calculation failure :-(  !!!')
            else:
                disp_('   ---> Solved.')
        disp_(' ')
    return s,resn,opt,exitc
