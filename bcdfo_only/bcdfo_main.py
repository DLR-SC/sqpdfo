#! /usr/bin/env python
from numpy import *
from bcdfo_hessP import *
from bcdfo_solve_TR_MS  import *
from bcdfo_augmX_evalf import *
from bcdfo_augment_Y import *
from bcdfo_swap_in_Y import *
from bcdfo_computeP import *
from bcdfo_gradP import *
from bcdfo_projgrad import *
from bcdfo_include_in_Y import *
from bcdfo_poisedness_Y import *
from bcdfo_repair_Y import *
from bcdfo_find_smallf import *
from bcdfo_solve_TR_MS_bc import *
from bcdfo_projected_tcg import *
from bcdfo_ident_active_bnds import *

def bcdfo_main(f, nitold, nit, i_xbest, xl, xu, m, X, fX, ind_Y, QZ, RZ, Delta, 
   cur_degree, neval, maxeval, maxit, model, gx, normgx, show_errg, pquad, pdiag, 
   plin, stallfact, eps_rho, Deltamax, rep_degree, epsilon, verbose, eta1, eta2, 
   gamma1, gamma2, gamma3, interpol_TR, factor_CV, Lambda_XN, Lambda_CP,
   factor_FPU, factor_FPR, Lambda_FP, criterion_S, criterion_FP, criterion_CP, 
   mu, theta, eps_TR, eps_L, lSolver, stratLam, eps_current, vstatus, xstatus, 
   sstatus, dstatus, ndummyY, sspace_save, xspace_save, xfix, fxmax, 
   poised_model, kappa_ill, kappa_th, eps_bnd, poised, Y_radius, c, level, 
   whichmodel, hardcons, noisy, scaleX, scalefacX, CNTsin, shrink_Delta, gmat, 
   hessian, hessian_bfgs, g_bfgs, useMethod, stage, Deltamin, scale, shift_Y):
###############################################################################
#
#  function bcdfo_main
#
#  Main routine of the algorithm BC-DFO
#
#  INPUTS:
#
#  see INPUTS of bcdfo() and in addition:
#
#  nitold       : iterations already done before restarting BCDFO
#  nit          : current number of iterations
#  i_xbest      : index of current best iterate
#  m            : cardinality of X
#  X            : set of all points
#  fX           : function values of X
#  ind_Y        : indices of the interpolation set Y out of the set X
#  neval        : number of function evaluations
#  model        : coefficients of the polynomial model
#  gx           : gradient at the current iterate
#  normgx       : infinity norm of the projected gradient
#  stallfact    : termination (stalled) when ||s|| <= stallfact * norm( x )
#  eps_rho      : stabilization constant for rho
#  eps_current  : current threshold for a criticality step
#  vstatus      : status of the variables (free/fixed)
#  xstatus      : points in Y or not
#  sstatus      : points in subspace or not
#  dstatus      : points are dummy points or not
#  ndummyY      : nbr of dummy points in the set Y
#  sspace_save  : all subspaces already explored
#  xspace_save  : indices of entering and exiting points corresponding to each
#                 subspace already explored
#  xfix         : n-dim. vector that stores the values of fixed variables
#  poised_model : it is currently known that model is well-poised
#  poised       : poisedness of the set
#  Y_radius     : current radius of the set Y
#  c            : contains a bunch of constants
#  level        : "toplevel" or "sublevel"
#  CNTsin       : variable to produce a quasi-random vector in bcdfo_choose_lin
#  scale        : the model diagonal scaling
#
#  OUTPUTS:
#
#  updated inputs and
#  msg          : (error) message
#
#  DEPENDENCIES: bcdfo_ident_active_bnds, bcdfo_gradP, bcdfo_hessP, bcdfo_projgrad,
#                bcdfo_poisedness_Y, bcdfo_augmX_evalf, bcdfo_find_smallf,
#                bcdfo_repair_Y, bcdfo_augment_Y, bcdfo_include_in_Y,
#                bcdfo_find_new_yj, bcdfo_replace_in_Y, bcdfo_print_summary_vector,
#                bcdfo_print_vector, bcdfo_solve_TR_MS_bc, bcdfo_projected_tcg,
#                bcdfo_swap_in_Y
#
#  PROGRAMMING: Ph. Toint, A. Troeltzsch and S. Gratton, 2009-2011.
#               (This version 18 V 2011)
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
###############################################################################

   #  Initialization

   msg               = 'Unexpected message from bcdfo_main';
   full_n, nbr_elem  = shape( X )
   m                 = nbr_elem - 1         # m is the last index of X, not the number of elements !!
   indfree           = where( vstatus == c.free() )[0]
#   print 'indfree:',indfree
   indfix            = where( vstatus >= c.fixed() )[0]
#   print 'indfix:',indfix
   nfix              = len(indfix);
#   print 'ind_Y:',ind_Y
   Y                 = X[indfree,:]
   Y                 = Y[:,ind_Y]
#   print 'Y:',Y
   x                 = X[indfree,i_xbest]
   x                 = x[:,newaxis]
#   print 'x:',x
   n                 = size(Y,0)
#   print 'n: ',n
   fY                = fX[ind_Y]
   fx                = fX[i_xbest]
   I                 = eye(full_n);
   poisedness_known  = 0;
   yk                = zeros((1,n));
   s                 = zeros((1,n));

   #  Initialize to zero for printing reasons after possible immediate convergence

   norms             = 0;
   oldDelta          = 0;
   rho               = 0;

   ###############################################################################
   #######################  MAIN ITERATION #######################################
   ###############################################################################

   for k in range(nitold,maxit):

      nit = nit + 1;
   
      ####################################
      #### Check for active bounds #######
      ####################################

      check_bnds = 0
 
      if ( normgx > epsilon and Delta > Deltamin and ndummyY == 0 and check_bnds):
      
         [nit, i_xbest, m, X, fX, QZ, RZ, Delta, eps_current, cur_degree, sstatus,
            xstatus, dstatus, ndummyY, sspace_save, xspace_save, ind_Y, Y, model,
            gx, normgx, neval, scale, CNTsin, gmat, hessian, hessian_bfgs, stage,
            msg] = bcdfo_ident_active_bnds(
         f, i_xbest, xl, xu, m, X, fX, ind_Y, QZ, RZ, Delta,
         cur_degree, neval, maxeval, nit, maxit, model, gx, normgx, show_errg, 
         pquad, pdiag, plin, stallfact, eps_rho, Deltamax, rep_degree, epsilon, 
         verbose, eta1, eta2, gamma1, gamma2, gamma3, interpol_TR, factor_CV, 
         Lambda_XN, Lambda_CP, factor_FPU, factor_FPR, Lambda_FP, criterion_S, 
         criterion_FP, criterion_CP, mu, theta, eps_TR, eps_L, lSolver, stratLam, 
         eps_current, vstatus, xstatus, sstatus, dstatus, sspace_save, xspace_save, 
         xfix, fxmax, poised_model, kappa_ill, kappa_th, eps_bnd, poised, Y_radius, 
         c, level, whichmodel, hardcons, noisy, scaleX, scalefacX, CNTsin, 
         shrink_Delta, gmat, hessian, hessian_bfgs, g_bfgs, useMethod, stage, Deltamin, scale, 
         shift_Y);

         x        = X(indfree,i_xbest);
         fx       = fX(i_xbest);
         fY       = fX(ind_Y);
         norms    = 0;
         oldDelta = Delta;

         if ( msg[0:5] == 'Error' ):
            if ( level == 'toplevel' ):
               disp(msg)

            return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval

         #  Check if current subspace has already been explored

         if ( not(isempty(xspace_save)) ):
            if ( len(where(xspace_save[2,:] > 0)[0]) > 0 and current_method == 0 ):
               if ( level == 'toplevel' ) :
                  xspace_save[2,:]=0;

                  #  If current subspace has already been explored,
                  #  repair interpolation set in smaller radius

                  maxd = 0;
                  for i in range(0,cur_degree):
                     maxd = max( linalg.norm(Y[:,i]-Y[:,0]), maxd )
               
                  radius = maxd * 0.1;

                  if ( verbose > 1 ):
                     print 'repair set and reduce TR because subspace already explored'
                     print 'new repair radius: '+str(radius)
               

                  #  Check if TR has not become unreasonably small

                  if ( Delta < stallfact * linalg.norm( x ) or Delta < Deltamin ):
                     return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval
               
                  effective_FPR = 1;
                  [ QZ, RZ, Y, replaced, poised, Y_radius, x, scale ] = bcdfo_repair_Y(
                     QZ, RZ, Y, radius, effective_FPR, Lambda_FP,
                     Lambda_CP, eps_L, x, lSolver, whichmodel, hardcons, xl, xu, 
                     indfree, stratLam, scale, shift_Y, normgx, kappa_ill );

                  poised_model = 1;

                  #  Compute the corresponding function values.

                  for i in range(0,len( replaced )):
                     j       = replaced( i );

                     # Set index of new point and update status of replaced point

                     m                 = m + 1
                     xstatus[ind_Y[j]] = c.unused()
                     ind_Y[j]          = m

                     #  Update X and evaluate function at Y(:,j)

                     if ( whichmodel >= 8 ):
                        compute_g = 1;
                     else:
                        compute_g = 0;
                  
                     [X, fX, neval, xstatus, sstatus, dstatus, gmat, msg] = bcdfo_augmX_evalf(
                     f, Y[:,j], m, X, fX, nfix, xfix, indfix,
                     indfree, fxmax, neval, maxeval, verbose, show_errg, xstatus, 
                     c.inY, sstatus, dstatus, scaleX, scalefacX, compute_g, gmat );
                     fY[j] = fX[m]   #  ??? not sure yet
                     if ( msg[0:5] == 'Error' ):
                        if ( level == 'toplevel' ):
                           disp( msg )
                        
                        return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval
             
                  #  Move to the best point found, if different from x.

                  i_xold = i_xbest;
                  [x, fx, QZ, RZ, Y, fY, ind_Y, i_xbest, scale] = bcdfo_find_smallf(
                     c, QZ, RZ, Y, fY, ind_Y, i_xbest, cur_degree,
                     indfree, x, xl, xu, fx, dstatus, whichmodel, scale, shift_Y, 
                     Delta, normgx, kappa_ill);
     
                  #  Compute the associated polynomial interpolation model.

                  model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, model, 
                        ind_Y, i_xold, m, gx, scale, shift_Y, Delta, 
                        indfree, gmat, eye(n), epsilon, noisy );
                  gx     = bcdfo_gradP( model, x, x, scale, shift_Y );
                  normgx,pgx = bcdfo_projgrad(n,x,gx,xl[0,indfree],xu[0,indfree]);

                  itype = 'repa';

                  if ( verbose >= 1 ):

                     if ( show_errg ):

                        #  Compute poisedness and gradient error, if unavailable so far.

                        if ( not poisedness_known ):
                           [ poised, Y_radius ] = bcdfo_poisedness_Y( QZ, RZ, Y, eps_L,
                              x, lSolver, whichmodel, hardcons, xl, xu, indfree, 
                              stratLam, scale, shift_Y );
                           poisedness_known     = 1;
                     
                        errg = poised * Y_radius / factor_CV;

                        #  Iteration summary

                        print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % (nit,
                        neval, fx, normgx, errg, 0.0, Delta, 0.0, n, itype )

                     else:

                        print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % (nit,
                        neval, fx, normgx, 0.0, Delta, 0.0, n, itype )

                  

                     # printout to file convhist.m

#                     if ( verbose ):
#                        fid = fopen('convhist.m','a');
#                        fprintf(fid,'#6d  #+.14e #.2e \n', neval, fx, normgx);
#                        fclose(fid);
                  
                     if ( verbose == 3 ):
                        bcdfo_print_summary_vector( 'x ', x )
                        bcdfo_print_summary_vector( 'gx', gx )
                     elif ( verbose > 3 ):
                        bcdfo_print_vector( 'x ', x )
                        bcdfo_print_vector( 'gx', gx )
                        bcdfo_print_vector( 's', s )
                  
                     poisedness_known = 0;           # reset for the next iteration

                     if ( verbose >= 2 ):
                        print ' cond(Z) = '+str( linalg.cond(RZ) )

                        if ( verbose >= 10 ):
                           bcdfo_print_vector( 'model', model );            

                  #  Set Delta

                  Delta    = radius;
                  oldDelta = Delta;

               else:

                  #  If in the subspace, return to the toplevel to repair
                  
                  return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval
            
            elif ( len(where(xspace_save[2,:] > 0)[0]) > 0 and useMethod == 2 ):
   
               #  Switch to BCDFO because of troubles in BFGS (tries to enter the same subspace)

               print ' HYBRID: wants to enter the same subspace.'  
   
      #  Return to toplevel if gradient or trust-region radius small in the sublevel
      #  (no serious check of gradient accuracy in the sublevel)

      if ( ( normgx <= epsilon or Delta <= Deltamin ) and level == 'sublevel' ):

         #  If dummy points in Y, recompute their real function values before returning

         if ( ndummyY > 0 ):

            dummiesY = where( (dstatus == c.dummy) and (xstatus == c.inY) );
            for i in range(0,len( dummiesY )):

               k = dummiesY[i]
               j = where( k == ind_Y );

               #  Update X and evaluate function

               if ( whichmodel >= 8 ):
                  compute_g = 1;
               else:
                  compute_g = 0;
            
               [X, fX, neval, xstatus, sstatus, dstatus, gmat, msg] = bcdfo_augmX_evalf(
                  f, Y[:,j], k, X, fX, nfix, xfix, indfix, indfree,
                  fxmax, neval, maxeval, verbose, show_errg, xstatus, c.inY, sstatus,
                  dstatus, scaleX, scalefacX, compute_g, gmat );
               fY[j] = fX[k]    # ??? not yet sure
               if ( msg[0:5] == 'Error' ):
                  if ( level == 'toplevel' ):
                     disp( msg )
                  
                  return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval
                     
            ndummyY  = 0;

            #  Move to the best point found, if different from x.

            i_xold = i_xbest;
            [x, fx, QZ, RZ, Y, fY, ind_Y, i_xbest, scale] = bcdfo_find_smallf(
               c, QZ, RZ, Y, fY, ind_Y, i_xbest, cur_degree, 
               indfree, x, xl, xu, fx, dstatus, whichmodel, scale, shift_Y, 
               Delta, normgx, kappa_ill);

            #  Compute the associated polynomial interpolation model.

            model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, model, ind_Y, 
                  i_xold, m, gx, scale, shift_Y, Delta, indfree, gmat, 
                  hessian, epsilon, noisy );
            gx     = bcdfo_gradP( model, x, x, scale, shift_Y );
            normgx,pgx = bcdfo_projgrad(n,x,gx,xl[0,indfree],xu[0,indfree]);

            if ( normgx <= epsilon or Delta <= Deltamin ):
               #  Return if gradient is still small
               return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval
         
         else:
            #  Return if no dummy in Y
            return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval
      
      
 
      ####################################
      ############ Noise test ############
      ####################################
   
      # test for termination due to noise (if noise in function expected)

      if ( noisy ):
      
         ndummies = len(where( (dstatus == c.dummy) and (xstatus == c.inY()) ));
         indY_wo_dummies = where( (dstatus == c.nodummy) and (xstatus == c.inY()) );

         minf = min(fX[indY_wo_dummies].any())
         maxf = max(fX[indY_wo_dummies].any())
         noise_in_f = abs( maxf - minf )
   
         #   disp(['actual  diff in fvalues: |maxf-minf|=', str(noise_in_f)]);
   
         Delta_Y = 0;
         print 'Y'
         print Y
         for i in range(0,len( ind_Y )):
            Delta_Y = max( Delta_Y, linalg.norm( Y[:,i] - x ));
         print 'Y'
         print Y
   
         noise_tol = normgx * Delta_Y;
   
         #   disp(['allowed diff in fvalues: ||g||*Delta(Y)=', str(noise_tol), ', Delta(Y)=',str(Delta_Y)]);
   
         if ( isempty(hessian) ):
            if ( current_method == 0 ):
               hessian = bcdfo_hessP( model, x, x, scale, shift_Y );
            else:
               hessian = bcdfo_hessbfgs( s, yk, hessian );
         
      
   
         noise_tol2 = normgx * Delta_Y + Delta_Y**2 / 2 * linalg.norm(hessian);
   
         #   disp(['allowed diff in fvalues: ||g||*Delta(Y)+Delta^2/2||H||=', str(noise_tol2), ', Delta(Y)=',str(Delta_Y)]);
   
         if ( (noise_in_f > noise_tol2) and (Delta_Y < 10 * epsilon) and (cur_degree > n+1) ):         
            print 'Warning: Function is probably noisy! Calculation should be stopped.'
   
      ####################################
      #### Test for repair/convergence ###
      ####################################

      if ( normgx <= eps_current or Delta <= Deltamin ):

         # If the model is at its requested degree for convergence,
         # either we have converged or the interpolation set should be repaired
         # for a smaller eps_current.

         augment = rep_degree - cur_degree;

         if ( verbose >= 10 ):
            print ' cur_degree ='+str( cur_degree )+' rep_degree = '+str( rep_degree )      

         # No additional interpolation point is required for non-terminal repair.

         if ( augment <= 0 ):

            #  Terminate if the solution has been found.

            if ( normgx <= epsilon or Delta <= Deltamin ):
            
#               if ( useMethod >= 1 ):
               
                  # Terminate without checking gradient accuracy
                  
               poised = 0
               Y_radius = 0
                  
#               else:

                  #  Compute poisedness of the interpolation set Y

#                  [ poised, Y_radius ] = bcdfo_poisedness_Y( QZ, RZ, Y, eps_L, x, 
#                     lSolver, whichmodel, hardcons, xl, xu, indfree, stratLam, 
#                     scale, shift_Y );
#                  poisedness_known     = 1;

               #  Compute gradient accuracy

               errg = poised * Y_radius / factor_CV
            
               if ( logical_and( logical_or( errg <= epsilon, Delta <= Deltamin ), level == 'toplevel') ):
                  itype = 'conv';
              
                  if ( verbose >= 1 ):
                     if ( show_errg ):
                        print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % (nit,
                        neval, fx, normgx, errg, norms, Delta, rho, n, itype )
                        if ( normgx <= epsilon and errg <= epsilon ):
                           print '******************************************** Converged ******************************************'
                        else:
                           print '************************************* Trust-region radius small *********************************'
                     
                     else:
                        print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % (nit,
                        neval, fx, normgx, norms, Delta, rho, n, itype )
                        if ( normgx <= epsilon and errg <= epsilon ):
                           print '*************************************** Converged *************************************'
                        else:
                           print '******************************** Trust-region radius small ****************************'
                                  
                     #  Printout to file convhist.m

#                     if ( verbose ):
#                        fid = fopen('convhist.m','a');
#                        fprintf(fid,'#6d  #+.14e #.2e \n', neval, fx, normgx);
#                        fclose(fid);
                  

                     if ( verbose == 3 ):
                        bcdfo_print_summary_vector( 'x ', x )
                        bcdfo_print_summary_vector( 'gx', gx )
                     elif ( verbose > 3 ):
                        bcdfo_print_vector( 'x ', x )
                        bcdfo_print_vector( 'gx', gx )
                  
               
                  msg = ' Convergence in '+str( neval )+' evaluations of the objective function.'
                  return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval
              
         #  Not at a solution: improve the interpolation set.

         itype = 'impr';

         #  Reduce eps_current if repair degree is reached.

         if ( augment <= 0 ):
            eps_current = max( mu*eps_current, epsilon );
      

         #  Rebuild a poised model in the eps_current ball, and ...

         if ( normgx <= epsilon ) :

            # ... remove far points (with strict criterion: farfact = 1, forcing all
            # new interpolation points to be within distance epsilon)

            effective_FPR = 1;

         else:

            # ... remove far points (with looser criterion)

            effective_FPR = factor_FPR;

      

         #  One needs to add new interpolation points to reach the desired degree.
         #  (only entered if rep_degree is higher than linear!)

         if ( augment > 0 ):

            itype = 'augm';

            #  If gradient small, find a new point in the epsilon-environment,
            #  not in Delta (distinguish between infty-norm and 2-norm local solver)
         
            if ( normgx <= epsilon ):
               if ( lSolver == 2 ):
                  Delta = epsilon/sqrt(n);
                  eps_current = epsilon/sqrt(n);
               else:
                  Delta = epsilon;
                  eps_current = epsilon;
                     
            #  Pick a random interpolation point.

            ynew = -Delta * ones((n,1)) + 2 * Delta * random.rand(n,1)
            [ cur_degree, QZ, RZ, Y, xbase, scale ] = bcdfo_augment_Y(
               ynew, Y[:,0:cur_degree], whichmodel, shift_Y,
               Delta, normgx, kappa_ill );
            ind_Y[cur_degree-1] = cur_degree-1

            #  Optimally replace it.

            if ( hardcons ):
               [ ynew, improvement ] = bcdfo_find_new_yj_bc( QZ, RZ, Y, cur_degree,
                  Delta, eps_L, xbase, lSolver, whichmodel, xl, xu, indfree, 
                  stratLam, scale, shift_Y );
            else:
               [ ynew, improvement ] = bcdfo_find_new_yj( QZ, RZ, Y, cur_degree,
                  Delta, eps_L, xbase, lSolver, whichmodel, scale, shift_Y );
            

            [ QZ, RZ, Y, xbase, scale ] = bcdfo_replace_in_Y( QZ, RZ, ynew, Y, 
               cur_degree, xbase, whichmodel, scale, shift_Y, Delta, normgx, 
               kappa_ill );
            replaced = copy(cur_degree)
  
         #  The current interpolation set has the requested degree.

         else:
          
            #  If gradient small, repair in epsilon-radius, else in Delta
            #  (distinguish between infty-norm and 2-norm local solver)

            if (normgx <= factor_CV * epsilon ):
               if ( lSolver == 2 ):
                  radius = min( Delta/sqrt(n), epsilon/sqrt(n) );
               else:
                  radius = min( Delta, epsilon );
            
            else:
               radius = max( Delta, eps_current );
                     
            [ QZ, RZ, Y, replaced, poised, Y_radius, x, scale ] = bcdfo_repair_Y(
               QZ, RZ, Y, radius, effective_FPR, Lambda_FP, 
               Lambda_CP, eps_L, x, lSolver, whichmodel, hardcons, xl, xu, 
               indfree, stratLam, scale, shift_Y, normgx, kappa_ill );

         if ( verbose >= 3 ):
            [ poised, Y_radius ] = bcdfo_poisedness_Y( QZ, RZ, Y, eps_L, x, lSolver, 
                  whichmodel, hardcons, xl, xu, indfree, stratLam, scale, shift_Y );
            poisedness_known     = 1;
            print ' poisedness(Y) = '+str(poised)
      
         poised_model = 1;

         #  Compute the corresponding function values.

         for i in range(0,len( replaced )):
            j       = replaced[i]
         
            #  Set index of new point and update status of the old point
         
            m                 = m + 1;
            xstatus[ind_Y[j]] = c.unused()
            ind_Y[j]          = copy(m)
         
            #  Update X and evaluate function

            if ( whichmodel >= 8 ):
               compute_g = 1;
            else:
               compute_g = 0;
         
            [X, fX, neval, xstatus, sstatus, dstatus, gmat, msg] = bcdfo_augmX_evalf(
               f, Y[:,j], m, X, fX, nfix, xfix, indfix, indfree,
               fxmax, neval, maxeval, verbose, show_errg, xstatus, c.inY(), sstatus,
               dstatus, scaleX, scalefacX, compute_g, gmat );
            fY[j] = fX[m]       # ??? not yet sure
            if ( msg[0:5] == 'Error' ):
               if ( level == 'toplevel' ):
                  disp( msg )
               
               return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval
         
         #  Move to the best point found, if different from x.

         i_xold = i_xbest;
         [x, fx, QZ, RZ, Y, fY, ind_Y, i_xbest, scale] = bcdfo_find_smallf(
            c, QZ, RZ, Y, fY, ind_Y, i_xbest, cur_degree, 
            indfree, x, xl, xu, fx, dstatus, whichmodel, scale, shift_Y, 
            Delta, normgx, kappa_ill);
      
         #  Compute the associated polynomial interpolation model.

         model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, model, 
               ind_Y, i_xold, m, gx, scale, shift_Y, Delta, indfree, 
               gmat, hessian, epsilon, noisy );
         gx     = bcdfo_gradP( model, x, x, scale, shift_Y );
         normgx,pgx = bcdfo_projgrad(n,x,gx,xl[0,indfree],xu[0,indfree]);

         #  Terminate if the solution has been found.

         errg = poised * Y_radius / factor_CV;
         if ( normgx / factor_CV <= epsilon and errg  <= epsilon and cur_degree >= rep_degree and level == 'toplevel' ):
            if ( verbose >= 1 ):
               if ( show_errg ):
                  print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % ( nit,
                  neval, fx, normgx, errg, norms, oldDelta, rho, n, itype )
                  print '******************************************** Converged ******************************************'
               else:
                  print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % ( nit,
                  neval, fx, normgx, norms, oldDelta, rho, n, itype )
                  print '*************************************** Converged *************************************'

               # printout to file convhist.m

#            if ( verbose ):
#               fid = fopen('convhist.m','a');
#               fprintf(fid,'#6d  #+.14e #.2e \n', neval, fx, normgx);
#               fclose(fid);            

               if ( verbose == 3 ):
                  bcdfo_print_summary_vector( 'x ', x )
                  bcdfo_print_summary_vector( 'gx', gx )
               elif ( verbose > 3 ):
                  bcdfo_print_vector( 'x ', x )
                  bcdfo_print_vector( 'gx', gx )            
         
            msg = ' Convergence in '+str( neval )+' evaluations of the objective function.'

            #  Return to calling routine

            return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval

         #  Reset the radius to a multiple of ||gx||.
      
         if ( augment <= 0 ):
            Delta = min( min( theta * normgx, epsilon ), Deltamax );      

         #  Print the details.

         if ( verbose >= 1 ):
            if ( show_errg ):
               print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % ( nit,
               neval, fx, normgx, errg, norms, Delta, rho, n, itype )
            else:
               print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % ( nit,
               neval, fx, normgx, norms, Delta, rho, n, itype )
         
            # printout to file convhist.m

#         if ( verbose ):
#            fid = fopen('convhist.m','a');
#            fprintf(fid,'#6d  #+.14e #.2e \n', neval, fx, normgx);
#            fclose(fid);
         

            if ( verbose == 3 ):
               bcdfo_print_summary_vector( 'x ', x )
               bcdfo_print_summary_vector( 'gx', gx )
            elif ( verbose > 3 ):
               bcdfo_print_vector( 'x ', x )
               bcdfo_print_vector( 'gx', gx )

         #  Start a new iteration.
     
         continue   

      ##############################################
      ####### Compute the trial point by ###########
      ####### solving the TR subproblem  ###########
      ##############################################

      #  Compute Hessian

      if ( useMethod == 0 ):
         hessian = bcdfo_hessP( model, x, x, scale, shift_Y );
      else:
         hessian = bcdfo_hessbfgs( s, yk, hessian );

      lb = xl[:,indfree] - x.transpose()
      ub = xu[:,indfree] - x.transpose()

      #  Solve the trust-region subproblem

#      print 'gx='
#      print gx
#      print 'hessian='
#      print hessian



      if ( lSolver == 1 ):       
         [ s, lambd, norms, value, gplus, nfact, neigd, msgTR ] = bcdfo_solve_TR_MS_bc(
            gx, hessian, lb, ub, Delta, eps_TR, stratLam );         
#         [ s, lambd, norms, value, gplus, nfact, neigd, msgTR, hardcase ] = bcdfo_solve_TR_MS(
#            gx, hessian, Delta, eps_TR );
         
      elif ( lSolver == 2 ):
         lb = maximum( lb, -Delta )
         ub = minimum( ub, Delta )
       
         [ s, msgTR, cgits, value, gplus, interior ] = bcdfo_projected_tcg(
            1, gx, hessian, lb, ub, 10, 1e-9, 1e-9, 0 )

         norms = linalg.norm( s[0], inf )


      prered  = -value
   
      if ( msgTR[0:5] == 'Error' ):
         msg = [ 'Error (number 0', str( msgTR[5] ), ') in local solver' ];
         if ( level == 'toplevel' ):
            disp( msg )
         
         return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model,gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval
   
      if ( interpol_TR == 1 ):     # save for computing backtracking interpolation of TR

         gTs  = dot(gx,s.T)[0][0]

      ############################################################
      ####### Evaluate the function at the trial point ###########
      ############################################################

      xplus  = x + s.T

      #  Set index of new point
         
      m = m + 1;
      
      #  Include point in X and evaluate f
      #  (xstatus(m) is set to 0 but is updated later on)

      if ( logical_or(useMethod > 0, whichmodel >= 7) ):
         compute_g = 1;
      else:
         compute_g = 0;
   
      [X, fX, neval, xstatus, sstatus, dstatus, gmat, msg] = bcdfo_augmX_evalf(
         f, xplus, m, X, fX, nfix, xfix, indfix, indfree,
         fxmax, neval, maxeval, verbose, show_errg, xstatus, 0, sstatus, dstatus,
         scaleX, scalefacX, compute_g, gmat );
      fxplus = fX[m]

      if ( msg[0:5] == 'Error' ):
         if ( level == 'toplevel' ):
            disp( msg )
         
         return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval
   
      rho     = ( fx - fxplus + eps_rho ) / ( prered + eps_rho );

      if ( prered <= 0 or abs( fx - fxplus ) < eps_rho ):
         rho = -1   

      if ( rho >= eta1 ):
         succ = 1
      else:
         succ = 0  
   
      if ( verbose >= 10 ):
         print ' rho = '+str( rho )
   
   
      ##################################################################
      ###### Include the new point in the interpolation set Y ##########
      ##################################################################
   
      oldDelta = Delta;                # save old radius for printing
      oldY     = Y;                    # save in case of bad condition
      oldRZ    = RZ;                   # save in case of bad condition
      oldQZ    = QZ;                   # save in case of bad condition
      oldscale = scale;                # save in case of bad condition
      i_xold   = i_xbest;              # save to compute min-frob-norm model
   
      pos = -1

      if ( useMethod > 0 ):

         itype = 'BFGS';

         if ( rho >= eta1 ):
         
            #  Successful BFGS-iteration

            if ( verbose >=2 ):
               print ' Successful iteration.'
         
            # Update the current iterate

            x        = copy(xplus)
            fx       = copy(fxplus)
            ind_Y[0] = copy(m)
            i_xbest  = copy(m)
            yk       = gmat[indfree,m] - gx;
            gx       = gmat[indfree,m]
            normgx,pgx = bcdfo_projgrad(n,x,gx,xl[0,indfree],xu[0,indfree]);
            pos      = 0

            # Update trust-region radius

            Delta = min( max( gamma3*norms, Delta ), Deltamax );

         else:

            yk       = gmat[indfree,m] - gx;

      elif ( useMethod == 0 ):

         #  Try to replace a dummy point if model degree is quadratic

         dummiesY  = array([],int)
         dummy_set = where( dstatus == c.dummy() )[0]
         
         if ( verbose == 2 ):
            print ' nbr of dummy points in X: ' + str(len(dummy_set))

         if ( dummy_set and (cur_degree >= pquad) ):

            for j in range(0,cur_degree):
               if ( len(nonzero( dummy_set == ind_Y[j])[0]) != 0 ):
                  dummiesY = append( dummiesY, j )
      
            if ( verbose == 2 ):
               print ' found the following dummies in Y: ' + str(dummiesY)
         
            #  Replace a dummy interpolation point using Lagrange maximization
            
            Lambda_DP    = 1e-10;
            criterion_DP = 'standard';
            
            # to do !!
            
            [ QZ, RZ, Y, pos, x, scale ] = bcdfo_include_in_Y( xplus, QZ, RZ, Y,  
               dummiesY, Lambda_DP, criterion_DP, x, whichmodel, succ, scale, 
               shift_Y, Delta, normgx, kappa_ill );
      
         ndummyY = len(dummiesY)
   
         #  If dummy point found to replace

         if ( pos > -1 ):
   
            itype  = 'repD';

            if ( verbose >= 2 ):
               print ' replacing interpolation point '+str( pos )+' (dummy)'
               condRZ_afterinclude = linalg.cond(RZ)
         

            #  Update status and position of the new point

            xstatus[ind_Y[pos]] = c.unused()
            xstatus[m]          = c.inY()
            ind_Y[pos]          = copy(m)
            fY[pos]             = copy(fxplus)
            ndummyY             = ndummyY - 1

            #  Swap points if included a successful point

            if ( rho >= eta1 ):
               [ QZ, RZ, Y, ind_Y, fY, x, scale ] = bcdfo_swap_in_Y(
                  0, pos, QZ, RZ, Y, ind_Y, fY, x, whichmodel,
                  scale, shift_Y, Delta, normgx, kappa_ill );
               fx           = copy(fxplus)
               i_xbest      = copy(m)
               
               if ( not shift_Y ):
                  x = Y[:,0]
            
               poised_model = 0;
               if ( verbose >= 2 ):
                  print ' swapped point to position 1'
            
         
               itype = 'repDs';
         
               #  Update the trust-region radius.

               Delta = min( max( gamma3*norms, Delta ), Deltamax );

            else:
         
               #  Shrink trust region in unsuccessful iteration

               if ( shrink_Delta == 1 and Delta > epsilon ):
                 Delta = gamma2 * Delta;    

            #  Compute the associated polynomial interpolation model.

            model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, model, ind_Y, 
               i_xold, m, gx, scale, shift_Y, Delta, indfree, gmat, 
               hessian, epsilon, noisy );
            gx     = bcdfo_gradP( model, x, x, scale, shift_Y );
            normgx,pgx = bcdfo_projgrad(n,x,gx,xl[0,indfree],xu[0,indfree]);
     
         #  No dummy point to replace

         else:   
      
            ########################################
            ####### Successful iteration ###########
            ########################################

            if ( rho >= eta1 ):

               if ( verbose >=2 ):
                  print ' Successful iteration.'

               #  Augment interpolation set if not fully quadratic yet

               if ( cur_degree < pquad or ( whichmodel == 3 and cur_degree < pquad+pquad) ):

                  [ cur_degree, QZ, RZ, Y, xbase, scale ] = bcdfo_augment_Y( xplus,  
                       Y, whichmodel, shift_Y, Delta, normgx, kappa_ill);
                  pos = cur_degree-1
                  
                  if ( pos > -1 ):
                     ind_Y   = append(ind_Y,m)
                     fY      = append(fY,fxplus)

               #  Include xplus in the interpolation set, by replacing another point if
               #  the model is already fully quadratic.

               else:

                  [ QZ, RZ, Y, pos, x, scale ] = bcdfo_include_in_Y( xplus, QZ, RZ, Y,
                     range(0,cur_degree), Lambda_XN, criterion_S, x, whichmodel, succ,
                     scale, shift_Y, Delta, normgx, kappa_ill );

                  if ( pos > -1 ):
                     xstatus[ind_Y[pos]] = c.unused()
                     ind_Y[pos]   = copy(m)
                     fY[pos]      = copy(fxplus)
                    

               # If xplus could/should be included in the interpolation set

               if ( pos > -1 ):

                  itype  = 'succ';

                  if ( verbose >= 2 ):
                     print ' replacing/including interpolation point '+str(pos)+' (successful)'
               
                  xstatus[m]   = c.inY()

                  #  Move it in the first position, redefining the base point.
                  
                  [ QZ, RZ, Y, ind_Y, fY, x, scale ] = bcdfo_swap_in_Y( 0, pos, QZ,
                     RZ, Y, ind_Y, fY, x, whichmodel, scale, shift_Y, Delta, 
                     normgx, kappa_ill );

                  fx           = copy(fxplus)
                  i_xbest      = copy(m)
                  
                  if ( shift_Y == 0 ):
                     x         = copy(Y[:,0])                     
                                    
                  poised_model = 0;

                  #  Compute the associated polynomial interpolation model.
                  
                  model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, model, 
                     ind_Y, i_xold, m, gx, scale, shift_Y, Delta, indfree, 
                     gmat, hessian, epsilon, noisy );

                  gx     = bcdfo_gradP( model, x, x, scale, shift_Y );
                  normgx,pgx = bcdfo_projgrad(n,x,gx,xl[0,indfree],xu[0,indfree]);

                  #  Update the trust-region radius.

                  Delta = min( max( gamma3*norms, Delta ), Deltamax );          
   
            ##########################################
            ####### Unsuccessful iteration ###########
            ##########################################
   
            #  Enter if iteration unsuccessful or the point could not be included in Y yet

            if ( rho < eta1 or pos == -1 ):

               if ( verbose >= 2 and rho < eta1 ):
                  print ' Unsuccessful iteration. '           
         
               #  The model is not fully quadratic yet: add (if possible)
               #  the new point to the interpolation set and recompute the model.

               if ( ( ( cur_degree < pquad ) or ( whichmodel == 3 and cur_degree < pquad+pquad) ) and ( rho < eta1 ) ):
               
                  [ cur_degree, QZ, RZ, Y, xbase, scale ] = bcdfo_augment_Y(
                     xplus, Y[:,0:cur_degree], whichmodel, shift_Y,
                     Delta, normgx, kappa_ill );
            
                  if ( verbose >= 2 ):
                     print ' including interpolation point '+str( cur_degree )+' (augm)'
               

                  # Update status and position of the new point
               
                  xstatus[m]        = c.inY()
                  ind_Y = append(ind_Y,m)
                  fY    = append(fY,fxplus)
                  poised_model      = 0

                  model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, model, 
                        ind_Y, i_xold, m, gx, scale, shift_Y, Delta, indfree, 
                        gmat, hessian, epsilon, noisy );
                  gx     = bcdfo_gradP( model, x, x, scale, shift_Y );
                  normgx,pgx = bcdfo_projgrad(n,x,gx,xl[0,indfree],xu[0,indfree]);

                  itype  = 'augm'
                  pos    = copy(m)

                  #  Shrink trust region in unsuccessful iteration

                  #if ( shrink_Delta == 1 and Delta > epsilon ):
                  #   Delta     = gamma2 * Delta;            
      
               #  Enter if the model is already fully quadratic *or*
               #  xplus could not yet be included in the set.
               #  The decision to include xplus here depends on possibly eliminating
               #  another point.

               if ( ( cur_degree >= pquad ) or ( pos == -1 ) ):

                  if ( ( pos == -1 ) and ( poised_model == 0 or Delta <= eps_current ) ):

                     #  Compute the distance of the interpolation points to the current
                     #  iterate. (Distinguish between the successful badcond and the
                     #  unsuccessful case!)

                     d = zeros((1,cur_degree))[0]
               
                     if ( rho >= eta1 ):
                        for j in range(0,cur_degree):
                           if ( lSolver == 1 ):
                              d[j] = linalg.norm( Y[:,j] - xplus[:,0] );
                           else:
                              d[j] = linalg.norm( Y[:,j] - xplus[:,0],inf );
                                                            
                     else:
                        
                        #  do not replace the current iterate
                        
                        for j in range(1,cur_degree):
                           if ( lSolver == 1 ):
                              d[j] = linalg.norm( Y[:,j] - x[:,0] );
                           else:
                              d[j] = linalg.norm( Y[:,j] - x[:,0],inf );
                  
                     #  Compute the basic distance used to define far/close points.

                     FPlength = factor_FPU * ( 1 + eps_TR ) * Delta;

                     #  Replace a far interpolation point.
                  
                     if ( rho >= eta1 ):
                        criterion_FPn = 'weighted';      # use weighted measure, not furthest point
                     else:
                        criterion_FPn = criterion_FP;

                     if ( where( d > FPlength )[0].any() ):
                        [ QZ, RZ, Y, pos, x, scale ] = bcdfo_include_in_Y( xplus, QZ, RZ, 
                           Y, where( d > FPlength )[0], Lambda_FP, criterion_FPn, x,
                           whichmodel, succ, scale, shift_Y, Delta, normgx, kappa_ill );

                     if ( pos > -1 ):
               
                        itype  = 'repF';
                  
                        if ( verbose >= 2 ):
                           print ' replacing interpolation point '+str( pos )+' (far)'                     

                        #  Update status and position of the new point

                        xstatus[ind_Y[pos]] = c.unused()
                        xstatus[m]          = c.inY()
                        ind_Y[pos]          = copy(m);
                        fY[ pos ]           = copy(fxplus);
                   
                        #  Swap points if included a successful point
         
                        if ( rho >= eta1 ):
                           [ QZ, RZ, Y, ind_Y, fY, x, scale ] = bcdfo_swap_in_Y( 0, 
                              pos, QZ, RZ, Y, ind_Y, fY, x, whichmodel, scale, 
                              shift_Y, Delta, normgx, kappa_ill );
                           fx           = copy(fxplus);
                           i_xbest      = copy(m);
                           if ( not shift_Y ):
                              x = copy(Y[:,0])
                        
                           poised_model = 0;
                           if ( verbose >= 2 ):
                              print ' swapped point to position 1'                        
                     
                           itype = 'repFs';
                     
                           #  Update the trust-region radius.

                           Delta = min( max( gamma3*norms, Delta ), Deltamax );
 
                        else:
   
                           #  Shrink trust region in unsuccessful iteration

                           if ( shrink_Delta == 1 and Delta > epsilon ):
                              Delta = gamma2 * Delta;
                                       
                        #  Compute the associated polynomial interpolation model.

                        model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, model, 
                           ind_Y, i_xold, m, gx, scale, shift_Y, Delta, 
                           indfree, gmat, hessian, epsilon, noisy );
                        gx     = bcdfo_gradP( model, x, x, scale, shift_Y );
                        normgx,pgx = bcdfo_projgrad(n,x,gx,xl[0,indfree],xu[0,indfree]);
                   
                  

                     # Replace a close interpolation point.

                     if ( pos == -1 ):

                        if ( rho >= eta1 ):
                           criterion_CPn = 'standard';    # find best improvement
                        else:
                           criterion_CPn = criterion_CP;
                     
                        if ( rho >= eta1 ):
                           Lambda_CPn = 1e-15;   # try hard to include a successful badcond point
                        else:
                           Lambda_CPn = Lambda_CP;
                           d[0] = 2 * FPlength;           # excludes the current iterate
                     
                  
                        [ QZ, RZ, Y, pos, x, scale ] = bcdfo_include_in_Y( xplus, QZ, 
                           RZ, Y, where( d <= FPlength )[0], Lambda_CPn, criterion_CPn, x,
                           whichmodel, succ, scale, shift_Y, Delta, normgx, kappa_ill );

                        if ( pos > -1 ):
                  
                           itype  = 'repC';

                           #  Safeguard frobenius model type 4 when replacing point 1

                           if ( pos == 0 ):
                              i_xold = copy(ind_Y[0])
                     
                           if ( verbose >= 2 ):
                              print ' replacing interpolation point '+str( pos )+' (close)'

                           #  Update status and position of the new point

                           xstatus[ind_Y[pos]] = c.unused()
                           xstatus[m]          = c.inY()
                           ind_Y[pos]          = copy(m)
                           fY[pos]             = copy(fxplus)
                     
                           #  Swap points if included a successful point
         
                           if ( rho >= eta1 ):
                              [ QZ, RZ, Y, ind_Y, fY, x, scale ] = bcdfo_swap_in_Y( 
                                 0, pos, QZ, RZ, Y, ind_Y, fY, x, whichmodel,
                                 scale, shift_Y, Delta, normgx, kappa_ill );
                              fx           = copy(fxplus)
                              i_xbest      = copy(m)
                              if ( not shift_Y ):
                                 x = copy(Y[:,0])
                           
                              poised_model = 0;
                        
                              if ( verbose >= 2 ):
                                 print ' swapped point to position 0'
                           
                        
                              itype = 'repCs';
                        
                              #  Update the trust-region radius.

                              Delta = min( max( gamma3*norms, Delta ), Deltamax );

                           else:

                              #  Shrink trust region in unsuccessful iteration

                              if ( shrink_Delta == 1 and Delta > epsilon ):
                                 Delta = gamma2 * Delta;                           
                        
                           #  Compute the associated polynomial interpolation model.

                           model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, model, 
                              ind_Y, i_xold, m, gx, scale, shift_Y, Delta, 
                              indfree, gmat, hessian, epsilon, noisy );
                           gx     = bcdfo_gradP( model, x, x, scale, shift_Y );
                           normgx,pgx = bcdfo_projgrad(n,x,gx,xl[0,indfree],xu[0,indfree]);
      
          # end currently used method

               # Decrease the radius.

               if ( pos == -1 ):

                  if ( verbose >= 2 ):
                     print ' decreasing the TR radius'

                  #  Set status of the new point

                  xstatus[m] = c.unused()

                  #  Compute new trust-region radius

                  if ( interpol_TR == 1 ):
                     curvature = - prered - gTs;
                     gam_inter = ( eta2 - 1 ) * gTs / ( fxplus - fx - gTs - eta2 * curvature );
                     Delta     = max( gamma1, min( gam_inter, gamma2 ) ) * min( Delta, norms );
                  else:
                     Delta     = gamma2 * Delta;
               
                  itype        = 'redD';

                  #  Check that the trust-region radius has not become so small that a step
                  #  of this size will not be significant.

                  if ( Delta < stallfact * linalg.norm( x ) and level == 'toplevel' ):

                     if ( verbose >= 1 ):
                        if ( show_errg ):
                           print '************************************* Trust-region radius small *********************************'
                        else:
                           print '******************************** Trust-region radius small ****************************'
                     
                        if ( verbose == 3 ):
                           bcdfo_print_summary_vector( 'x ', x )
                           bcdfo_print_summary_vector( 'gx', gx )
                        elif ( verbose > 3 ):
                           bcdfo_print_vector( 'x ', x )
                           bcdfo_print_vector( 'gx', gx )
                     
                  
                     msg = 'Algorithm stopped after '+str( neval )+' evaluations of the objective function because Delta small.'
                     return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval

      ########################################
      ######## Iteration printout ############
      ########################################

      if ( verbose >= 1 ):

         if ( show_errg ):

            #  Compute poisedness and gradient error, if unavailable so far.

            if ( not poisedness_known ):
               [ poised, Y_radius ] = bcdfo_poisedness_Y( QZ, RZ, Y, eps_L, x, 
               lSolver, whichmodel, hardcons, xl, xu, indfree, stratLam, scale, shift_Y );
               poisedness_known     = 1;
         
            errg = poised * Y_radius / factor_CV;

            #  Iteration summary

            print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % (nit,
            neval, fx, normgx, errg, norms, oldDelta, rho, n, itype )

         else:

            print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % (nit,
            neval, fx, normgx, norms, oldDelta, rho, n, itype )

         # printout to file convhist.m

#      if ( verbose ):
#         fid = fopen('convhist.m','a');
#         fprintf(fid,'%6d  %+.14e %.2e' % neval, fx, normgx)
#         fclose(fid);
      
         if ( verbose == 3 ):
            bcdfo_print_summary_vector( 'x ', x )
            bcdfo_print_summary_vector( 'gx', gx )
         elif ( verbose > 3 ):
            bcdfo_print_vector( 'x ', x )
            bcdfo_print_vector( 'gx', gx )
            bcdfo_print_vector( 's', s )
      
         poisedness_known = 0;           # reset for the next iteration

         if ( verbose >= 2 ):
            # print ' cond(Z) = '+str( linalg.cond(RZ) )

            if ( verbose >= 10 ):
               bcdfo_print_vector( 'model', model );

      #  Test for non-suitable model
      
      if ( useMethod == 0 and ( isnan(model).any() or not( isreal(model).all()) or isinf(model).any() )):

         msg = 'Error: model contains NaN or Inf or nonreal components!!';
         if ( verbose and level == 'toplevel' ):
            disp(msg);
      
         return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval
   

      #  Test for the maximum number of function evaluations.

      if ( neval >= maxeval ):
         if ( verbose >= 1 ):
            print ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAX. EVALUATIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      
         msg  = 'Error: Maximum number of '+str( maxeval )+' function evaluations reached. (bcdfo_main)'
         
         return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval
   

   ###############################################################################
   ################  END OF MAIN ITERATION #######################################
   ###############################################################################

   if ( verbose >= 1 ):
      if ( show_errg ):
         print ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAX. ITERATIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      else:
         print ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAX. ITERATIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

   msg  = 'Error: Maximum number of '+str( maxit )+' iterations reached.'

   return nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval
   
   
