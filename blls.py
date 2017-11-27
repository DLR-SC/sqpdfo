# -*- coding: utf-8 -*-

#try:
from runtime import *
from numpy import arange,append
#except ImportError:
    #from smop.runtime import *
def blls_(A=None,b=None,lb=None,ub=None,*args,**kwargs):
    """

#  blls is a solver for bound-constrained linear least-squares problems, that
#  is problem where, given 
#  - a m x n matrix A,
#  - a m x 1 right-hand side vector b,
#  - a n x 1 vector of (possibly -Inf) lower bounds lb,
#  - a n x 1 vector of (possibly +Inf) upper bounds ub,
#  one seeks a solution vector s solving
#
#           min || As - b ||^2     subject to    lb <= s <= ub
#
#  where ||.|| is the Euclidean norm and the contraint is understood 
#  componentwise.  No restriction is made on the dimension of A or its
#  rank, except that n>0 and m>0.
#
#  On exit, s contains the problem's solution if exitc = 0. If exitc = 1, 
#  s contains the best approximate solution found in maxiter = 2*n iterations.
#  The scalar resn contains the associated residual's norm, ie. ||As-b||.
#  If exitc = -1, the algorithm failed in the sense that the Cauchy step (see
#  below) could not be computed in maxbck backtackings (this error sould never
#  happen).
#
#  The maximum number of iterations, as other algorithmic parameters
#  (including tolerances for primal and dual feasibility and verbosity level),
#  may be modified in the "Algorithmic parameters" section in the beginning of
#  the code.
#
#  The method is intended for small-dimensional problems.  It is an active-set
#  algorithm where the unconstrained problem is solved at each iteration in
#  the subspace defined by the currently active bounds, themselves being
#  determined by a projected Cauchy step. Each subspace solution is computed
#  using a SVD decomposition of the reduced matrix.
#
#  Programming : Ph. Sampaio, Ph. Toint, A. Troeltzsch, April 2014.


    """
#    varargin = cellarray(args)
#    nargin = 4-[A,b,lb,ub].count(None)+len(args)

    [m,n]   = size_(A,nargout=2);        # problem size

    # Algorithmic parameters

    verbose = 0;              # verbosity flag (0: silent, larger -> more output)
    epsfeas = 1.0e-12;        # the tolerance for bound feasibility (primal)
    epsconv = 1.0e-10;        # the tolerance for optimality        (dual)
    epsres  = 1.0e-11;        # The tolerance for zero residual
    maxiter = 2*n;            # the maximum number of iterations
    epsdzer = 1.0e-14;        # the tolerance for considering a component of
                              # a search direction to be zero
    armijob = 0.5;            # the backtracking factor in the Armijo search
    armijor = 0.01;           # the linear reduction factor in the Armijo search
    maxback = 15;             # the maximum number of backtracks in the Armijo
                              # searches for and beyond the Cauchy point. This
                              # value and armijob should be chosen such that 
                              # armijob^maxback is of the order of the desired 
                              # precision on the variables.

    # Initialization

    inds  = arange(0,n);
    nit   = 0;                # the iteration counter
    nuns  = 0;                # the number of unsuccessful bactrackings (if
                              # set to 1, disables backtracking beyond the
                              # Cauchy point)

    ssub  = zeros_(n,1);       # prepare the shape of the subspace solution
    exitc = 0;                # initialize exit condition to successful

    # Start from the projection of the unconstrained minimizer 
    # on the feasible set.    
    
    s=min_(max_(pinv_(A).dot(b),lb),ub)
    # Compute the associated residual and optimality measure.

    res   = A.dot(s) - b;                                 # initial residual    
    resn  = norm_(res);                               # initial residual norm
    g     = A.T.dot(res);                                  # initial gradient
    stry  = min_(max_(s - g,lb),ub) - s;
    opt   = norm_( stry );                            # initial optimality measure
    
    # Compute the activity status of each variables.
    
    free  = find_( stry );                           # indices of free variables 
    atlb  = find_( abs( max_ (s-g,lb)-s) <= epsfeas); # inds of vars at lower bound
    atub  = find_( abs( min_ (s-g,ub)-s) <= epsfeas); # inds of vars at upper bound
    latlb = length_(atlb);     # the number of variables at their lower bound
    latub = length_(atub);     # the number of variables at their upper bound
    lfree = length_(free);     # the number of free variables
    
    # Print output banner and information on the initial iterate.
    
    
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
        disp_('     The problem has ',str(n),' variables and ',str(m),' rows.')
        disp_(' ')
        if (verbose > 2):
            print("problem_matrix="+str(A))
            print("right_hand_side="+str(b.T))
            print("lower_bounds="+str(lb.T))
            print("upper_bounds="+str(ub.T))
            disp_(' ')
        fprintf_('     nit     ||r||    optimality')
        fprintf_('                 nfr nlow nupp\n\n')
        fprintf_('   %5d  %.4e  %.4e                %4d %4d %4d\n' % (nit,resn,opt,lfree,latlb,latub))
        if (verbose > 1):
            if (verbose > 2):
                print("unconstrained_solution="+str(s.T))
            disp_(' ')
            disp_('   --------------------------------------------------------------')
            disp_(' ')
   # Terminate if the true solution is the initial one.

    if (opt <= epsconv or resn <= epsres):
        maxit=0 # Skips useless iterations and branches to final printout.
    else:
        maxit=copy(maxiter)
 # Main iteration loop
    for i in arange(1,maxit+1): #was indexed 1:maxit in matlab
        nit=nit + 1
   ##########################################################################
   #        Compute a Cauchy point using a projected Armijo search.         #
   ##########################################################################
        if (verbose > 1):
            disp_('   Iteration ',str(i))
            disp_(' ')
            fprintf_('   Cauchy point projected search\n')
            fprintf_('       k     ||r||     stepsize')
            fprintf_('                  nfr nlow nupp\n')
            fprintf_('   %5d  %.4e                            %4d %4d %4d'%(0,resn,lfree,latlb,latub))
        g=A.T.dot(res) # the quadratic's gradient
        alpha=(norm_(g) / norm_(A.dot( g))) ** 2   # the minimizing stepsize
        for kc in arange(1,maxback+1):  # Cauchy point loop
            stry=min_(max_(s - alpha*g,lb),ub)
            dtry=stry - s
            ltry=g.T.dot(dtry)
            qtry=ltry + 0.5 * norm_(A .dot(dtry)) ** 2   # the new quadratic value
            if (verbose > 1):
                fprintf_('\n   %5d  %.4e  %.4e'&(kc,norm_(A * stry - b),alpha))
            if (qtry <= armijor * ltry):  # the Armijo condition
                break
            else:
                alpha=armijob * alpha # backtracking
      # Exit if the Cauchy point calculation fails (this should not happen!
        if (kc >= maxback):
            exitc=- 1
            break
      # Determine the bounds which are active at the Cauchy point.

        atlb=inds[find_(abs(stry - lb) <= epsfeas)].reshape(-1)
        atub=inds[find_(abs(stry - ub) <= epsfeas)].reshape(-1)
        atb=append(atlb,atub)
        latlb=length_(atlb)
        latub=length_(atub)
        free=copy(inds)
        free = np.delete(free, atb)
        lfree=length_(free)
      # Fix the active bound and compute the residual 
      # at the Cauchy point.
        s[atlb]=lb[atlb]
        s[atub]=ub[atub]
        s[free]=stry[free]
        res=A.dot(s) - b
        resn=norm_(res)
        if (verbose > 1):
            fprintf_('                %4d %4d %4d\n'%(lfree,latlb,latub))
            if (verbose > 2):
                print("Cauchy_point="+str(s.T))
                print("indices_of_free_variables="+str(free))
                print("indices_of_variables_at_their_lower_bound="+str(atlb))
                print("indices_of_variables_at_their_upper_bound="+str(atub))
   ##########################################################################
   #       Beyond the Cauchy point: nested subspace minimization            #
   ##########################################################################
   # The Cauchy point is optimal for the problem because either there is no 
   # free variable left (which implies opt = 0) or because the residual is 
   # zero: terminate as s is the problem's solution.
        if (lfree == 0 or resn <= epsres):
            if (verbose > 1):
                fprintf_('   No nested subspace search\n')
            opt=0  #  Forces termination
      # There are free variables at the Cauchy point: start a nested 
      # subspace search.
        else: # ( lfree > 0 && resn > epsres )
           # Print initial information for nested subspace search, if requested.
 
           if (verbose > 1):
                fprintf_('   Nested subspace search\n')
                fprintf_('       k     ||r||     stepsize      ||r*||')
                fprintf_('      nfr nlow nupp\n')
                fprintf_('   %5d  %.4e                            %4d %4d %4d\n'&(0,resn,lfree,latlb,latub))
            # Loop on the successive nested subspaces.

           for k in arange(0,n):
             # Solve in the subspace of currently free variables.

                if (verbose > 2):
                    disp_('    > Solving in subspace ',str(k))
                    print("indices_of_free_variables="+str(free))
                    print("indices_of_variables_at_their_lower_bound="+str(atlb))
                    print("indices_of_variables_at_their_upper_bound="+str(atub))
                rhs=copy(b)
                if not isempty_(atlb):
                    rhs=rhs - A[0:m,atlb].dot(lb[atlb])
                if not isempty_(atub):
                    rhs=rhs - A[0:m,atub].dot(ub[atub])
              # Build the full-space version of the subspace minimizer
              # and compute the associated residual and its norm.
                ssub[free]=pinv_(A[0:m,free]).dot(rhs)
                ssub[atlb]=lb[atlb]
                ssub[atub]=ub[atub]
                rsubo=A .dot(ssub) - b
                rsubon=norm_(rsubo)
             # Check feasibility of the subspace solution wrt to inactive bounds.

                natlb=find_(ssub[free] < lb[free])
                natub=find_(ssub[free] > ub[free])
                lnatb=length_(natlb) + length_(natub)
            # If the subspace minimizer is unfeasible, attempt a projected
            # backtracking Armijo search until (at worst) back in the 
            # current subspace.
                if (lnatb > 0):
                    alpha=1
                    dtry=ssub - s  # the direction to the free minimizer
                    rred=rsubon - resn # associated residual reduction
                    found=0  # backtracking success indicator


                  # Backtracking loop
                    nback=4 * (1 - nuns)   # completely heuristic here! # any better idea?
                    for kb in arange(1,nback+1):
                        # Compute the unprojected trial point.

                        stry=(1 - alpha) * s+ alpha * ssub
                        # See which of the inactive bound are violated 
                        # at the unprojected trial point.
                        natlbt=free[find_(stry[free] < lb[free])]
                        natubt=free[find_(stry[free] > ub[free])]
                        lnatbt=length_(natlbt) + length_(natubt)
                        # Project the trial point.

                        stry=min_(max_(stry,lb),ub)
                        # Print information on the new point, if requested.

                        if (verbose >= 1):
                            rtry=A.dot(stry) - b
                            rtryn=norm_(rtry)
                            atlb=append(atlb,natlbt)
                            atub=append(atub,natubt)
                            atb=append(atlbt,atubt)
                            freet=copy(inds)
                            freet = np.delete(free, atbt)
                            latlbt=length_(atlbt)
                            latubt=length_(atubt)
                            lfreet=length_(freet)
                            fprintf_('   %5dp %.4e  %.4e   %.4e   %4d %4d %4d\n'%(kb,rtryn,alpha,rsubon,lfreet,latlbt,latubt))
                            # If the unprojected trial point is in the current subspace, 
                            # terminate backtracking. 
                        if (lnatbt == 0):
                            break
                        # If sufficient reduction in residual norm is obtained,
                        # the computed point is a suitable next iterate.  Branch out.
 
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
                            found=1  # Current point ok as next iterate
                            break
                        # Prepare for one more step of backtracking.
 
                        alpha=armijob * alpha
                # A suitable iterate has been found by backtracking: no further
                # action is needed.

                    if (found):
                        break

                # The backtracking search failed to return a suitable iterate:
                # compute the maximum feasible step toward the minimizer in the
                # subspace determined by the free variables at the Cauchy point.
                    else:
                        if (kb >= nback):
                            nuns=nuns + 1 #Increment failure counter
                        alpha=1
                        for kf in arange(0,length_(free)):   # Compute max step along dtry
                            kk=free[kf]
                            if (dtry[kk] >= epsdzer):
                                alpha=min_(alpha,(ub[kk] - s[kk]) / dtry[kk])
                            else:
                                if (dtry[kk] <= - epsdzer):
                                    alpha=min_(alpha,(lb[kk] - s[kk]) / dtry[kk])
                        ssub=s + alpha * dtry   # Take max step along dtry
                        rsub=(1 - alpha) * res + alpha * rsubo # residual
                        rsubn=norm_(rsub)
                        # Print information on the new point, if requested.

                        if (verbose > 1):
                            fprintf_('   %5ds %.4e  %.4e   %.4e   %4d %4d %4d\n'%(k,rsubn,alpha,rsubon,lfree,latlb,latub))

                       # Determine which new bounds are active and add them to 
                       # the current active set.
                        natlb=free[find_(abs(ssub[free] - lb[free]) <= epsfeas)].reshape(-1)
                        natub=free[find_(abs(ssub[free] - ub[free]) <= epsfeas)].reshape(-1)
                        atlb=append(atlb,natlbt)
                        atub=append(atub,natubt)
                        atb=append(atlb,atub)
                        free=copy(inds)
                        free = np.delete(free, atb)
                        latlb=length_(atlb)
                        latub=length_(atub)
                        lfree=length_(free)
                        # Print detailed activity information for the new subspace 
                        # iterate, if requested.
                        if (verbose > 2):
                            print("current_subspace_solution="+str(ssub.T))
                            print("indices_of_variables_at_their_lower_bound="+str(atlb))
                            print("indices_of_variables_at_their_upper_bound="+str(atub))
                            print("indices_of_free_variables="+str(free))

                       #  Prepare the next subspace minimization by memorizing the
                       #  current solution.
                        s=copy(ssub)
                        res=copy(rsub)
                        resn=copy(rsubn)
                # If the subspace minimizer is feasible, accept it as the next iterate.

                else:
                    s=copy(ssub)
                    res=copy(rsubo)
                    resn=copy(rsubon)
                    if (verbose > 1):
                        fprintf_('   %5df %.4e  %.4e   %.4e   %4d %4d %4d\n'%(k,resn,1,resn,lfree,latlb,latub))
                        if (verbose > 2):
                            print("current_subspace_solution="+str(ssub.T))
                            print("indices_of_variables_at_their_lower_bound="+str(atlb))
                            print("indices_of_variables_at_their_upper_bound="+str(atub))
                            print("indices_of_free_variables="+str(free))
                    break
            # Compute the optimality measure.

           opt=norm_(min_(max_(s - A.T.dot(res),lb),ub) - s)
       # Iteration printout

        if (verbose == 1):
            fprintf_('   %5d  %.4e  %.4e                %4d %4d %4d\n'%(nit,resn,opt,lfree,latlb,latub))
        else:
            if (verbose > 1):
                disp_(' ')
                fprintf_('     nit    ||r||     optimality')
                fprintf_('                 nfr nlow nupp\n\n')
                fprintf_('   %5d  %.4e  %.4e                %4d %4d %4d\n'%(nit,resn,opt,lfree,latlb,latub))
                if (verbose > 2):
                    print("current_solution="+str(s.T))
                    print("indices_of_free_variables="+str(free))
                    print("indices_of_variables_at_their_lower_bound="+str(atlb))
                    print("indices_of_variables_at_their_upper_bound="+str(atub))
                disp_('   --------------------------------------------------------------')
                disp_(' ')
        # Optimality termination test
        if (opt <= epsconv or resn <= epsres):
            break
############################################################################
#                         Finish off before returning                      #
############################################################################

# See if the maximum number of iterations has been reached without other
# problem. If yes, set the exit condition accordingly.
    if (exitc == 0 and nit >= maxiter):
        exitc=1
# Termination printout
    if (verbose > 0):
        disp_(' ')
        if (verbose > 2):
            print("indices_of_free_variables="+str(free))
            print("indices_of_variables_at_their_lower_bound="+str(atlb))
            print("indices_of_variables_at_their_upper_bound="+str(atub))
            print("final_solution="+str(s.T))
            print("final_residual="+str(res.T))
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
