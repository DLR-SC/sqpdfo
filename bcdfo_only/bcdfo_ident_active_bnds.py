#! /usr/bin/env python
from numpy import *

def bcdfo_ident_active_bnds(f, i_xbest, xl, xu, m, X, fX, ind_Y, QZ, RZ, Delta, 
   cur_degree, neval, maxeval, nit, maxit, model, gx, normgx, show_errg, pquad,
   pdiag, plin, stallfact, eps_rho, Deltamax, rep_degree, epsilon, 
   verbose_main, eta1, eta2, gamma1, gamma2, gamma3, interpol_TR, factor_CV, 
   Lambda_XN, Lambda_CP, factor_FPU, factor_FPR, Lambda_FP, criterion_S, 
   criterion_FP, criterion_CP, mu, theta, eps_TR, eps_L, lSolver, stratLam, 
   eps_current, vstatus, xstatus, sstatus, dstatus, sspace_save, xspace_save,
   xfix, fxmax, poised_model, kappa_ill, kappa_th, eps_bnd, poised, Y_radius,
   c, level, whichmodel, hardcons, noisy, scaleX, scalefacX, CNTsin, 
   shrink_Delta, gmat, hessian, hessian_bfgs, g_bfgs, useMethod, stage, Deltamin,
   scale, shift_Y):
###############################################################################
#
#  function bcdfo_ident_active_bnds
#
#  Identify active and nearly-active bounds, build a model in the subspace
#  of the remaining free variables and call BC-DFO recursively.
#
#  INPUT:
#
#  see inputs of function bcdfo() and:
#
#  i_xbest        : index of current best iterate
#  m              : cardinality of X
#  X              : set of all points
#  fX             : function values of X
#  ind_Y          : indices of the interpolation set Y out of the set X
#  neval          : number of function evaluations
#  nit            : number of iterations
#  model          : coefficients of the polynomial model
#  gx             : gradient at the current iterate
#  normgx         : infinity norm of the projected gradient
#  stallfact      : termination (stalled) when ||s|| <= stallfact * norm( x )
#  eps_rho        : stabilization constant for rho
#  vstatus        : status of the variables (free/fixed)
#  xstatus        : points in Y or not
#  sstatus        : points in subspace or not
#  dstatus        : points are dummy points or not
#  ndummyY        : nbr of dummy points in the set Y
#  sspace_save    : all subspaces already explored
#  xspace_save    : the entering and exiting points corresponding to each subspace
#  xfix           : n-dim. vector that stores the values of fixed variables
#  poised_model   : it is currently known that model is well-poised
#  poised         : poisedness of the set
#  Y_radius       : current radius of the set Y
#  c              : structure which contains a bunch of constants
#  level          : "toplevel" or "sublevel"
#  scale          : scaling factor of the matrix Z
#  CNTsin         : variable to produce a quasi-random vector in bcdfo_choose_lin
#  gmat           : a matrix containing the computed gradients
#
#  OUTPUT:
#
#  updated inputs and
#  msg            : (error) message
#
#  DEPENDENCIES: bcdfo_choose_lin, bcdfo_main, bcdfo_gradP, bcdfo_projgrad,
#                bcdfo_augmX_evalf, bcdfo_evalP, bcdfo_find_smallf,
#                bcdfo_build_QR_of_Y
#
#  PROGRAMMING: A. Troeltzsch, S. Gratton, December 2009.
#               (This version 3 XI 2011)
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
###############################################################################

   #  Initialization

   verbose      = 0
   msg          = 'empty message'              # return message
   full_n,m     = shape( X )                   # full dimension / nbr of points in X
   indfree      = where( vstatus == c.free() )[0]    # currently free variables
   indfix       = where( vstatus >= c.fixed() )[0]   # currently fixed variabled
   n            = len( indfree )            # dimension of current (sub)space
   nfix         = full_n - n                   # number of fixed variables
   Y            = X[indfree,:]
   Y            = Y[:,ind_Y]                # interpolation set
   fY           = fX[ind_Y]                  # function values of the points in Y
   x            = X[indfree,i_xbest]        # current best point x
   x            = x[:,newaxis]
   xl_free      = xl[0,indfree]                # lower bounds in the current subspace
   xu_free      = xu[0,indfree]                # upper bounds in the current subspace
   pgnorm,pgx   = bcdfo_projgrad( n, x, gx, xl_free, xu_free )
   ndummyY      = 0
   
   if ( verbose_main >= 2 ):
      verbose = 1
   
   
   if ( verbose ):
       print '****** enter ident_active_bounds: nfree = '+str(n)+' ****'
   
   
   #  Initiate new subspace
   
   xnew      = copy(x)
   xbd       = zeros((n,1))
   i_free    = array([])
   I_L       = array([])
   I_U       = array([])
   nbr_ss    = size(sspace_save,2)
   
   if (strcmp(level,'toplevel')):
      sspace = zeros((full_n))
   else:
      sspace = sspace_save[:,nbr_ss-1]
   
   #  Identify active and nearly-active bounds.
   #  If a bound is only nearly-active, move current point to the bound
   #  and save it as xnew. Save value of the bound in xbd.
   
   for i in range(0,n):
      ii = indfree[i]
      if ( logical_and( x[i,0] - gx[0,i] < xl_free[i], x[i,0] - xl_free[i] <= min(eps_bnd,pgx[0,i]) ) ):
         I_L         = append(I_L, i)
         xnew[i,0]     = xl_free[i]
         xbd[i,0]      = xl_free[i]
         vstatus[ii] = c.fixed()
         sspace[ii]  = -1;
      elif ( logical_and( x[i,0] - gx[0,i] > xu_free[i], xu_free[i] - x[i,0] <= min(eps_bnd,pgx[0,i]) ) ):
         I_U         = append(I_U, i)
         xnew[i,0] 	   = xu_free[i]
         xbd[i,0]      = xu_free[i]
         vstatus[ii] = c.fixed()
         sspace[ii]  = 1;
      else:
         i_free      = append(i_free, i)
      
   
   I_LandU = sort( append( I_L, I_U ) )
   
   if ( verbose and I_LandU ):
      print 'lower active bounds: '+str(I_L)
      print 'upper active bounds: '+str(I_U)
   

   #  Return if no active bounds

   if ( len( I_LandU ) == 0 ):
      msg = 'no active bounds'
      if ( verbose ):
         print msg
         print '****** return from ident_active_bounds **********'
   
      return nit, i_xbest, m, X, fX, QZ, RZ, Delta, eps_current, cur_degree, sstatus, xstatus, dstatus, ndummyY, sspace_save, xspace_save, ind_Y, Y, model, gx, normgx, neval, scale, CNTsin, gmat, hessian, hessian_bfgs, stage, msg


   #  Return if new subspace was already tried before from this position

   xspace_save[2,:] = 0
   for i in range(0,size(sspace_save,2)):
      if ( logical_and( linalg.norm(sspace - sspace_save[:,i]) == 0, logical_or(i_xbest == xspace_save[0,i], i_xbest == xspace_save[1,i]) ) ):
         msg = 'subspace already explored'
         xspace_save[2,i] = 1
         
         if ( verbose ):
            print msg
            print '****** return from ident_active_bounds **********'
      
         return nit, i_xbest, m, X, fX, QZ, RZ, Delta, eps_current, cur_degree, sstatus, xstatus, dstatus, ndummyY, sspace_save, xspace_save, ind_Y, Y, model, gx, normgx, neval, scale, CNTsin, gmat, hessian, hessian_bfgs, stage, msg
   
   

   #  Save actual subspace and index of current iterate
   
   sspace_save             = append( sspace_save, sspace )
   xspace_save[0,nbr_ss] = i_xbest;
   xspace_save[1,nbr_ss] = 0;          # place holder for index of exiting point
     
   #  Check if x was moved to the bound
   #  and evaluate f at xnew if necessary
   
   if ( linalg.norm(x - xnew) > 1e-10 ):
   
      #  Check if new point coincides with another point
   
      coincides = 0;
      nodummyX = where( dstatus == c.nodummy() )[0]
   
      for i in range(0,len(nodummyX)):
         ii = nodummyX[i]
         if ( linalg.norm(xnew - X[indfree,ii]) < 1e-10 ):
            coincides = nodummyX[i]
            break         
   
      if ( coincides ):
         msg = 'current iterate moved to the bounds but this point already exists';
         if ( verbose ):
           disp(msg)
      
         return nit, i_xbest, m, X, fX, QZ, RZ, Delta, eps_current, cur_degree, sstatus, xstatus, dstatus, ndummyY, sspace_save, xspace_save, ind_Y, Y, model, gx, normgx, neval, scale, CNTsin, gmat, hessian, hessian_bfgs, stage, msg
   

      #  Set index of new point

      m = m + 1;

      #  Update X and evaluate f

      if ( logical_or( whichmodel >= 7, useMethod > 0 ) ):
         compute_g = 1
      else:
         compute_g = 0
  
      [X, fX, neval, xstatus, sstatus, dstatus, gmat, msg] = bcdfo_augmX_evalf(f, xnew,
         m, X, fX, nfix, xfix, indfix, indfree, fxmax, neval, maxeval, verbose, 
         show_errg, xstatus, 0, sstatus, dstatus, scaleX, scalefacX, compute_g, gmat)
      if ( msg[0:5] == 'Error'  ):
         if ( level == 'toplevel' ):
            print msg
      
         return nit, i_xbest, m, X, fX, QZ, RZ, Delta, eps_current, cur_degree, sstatus, xstatus, dstatus, ndummyY, sspace_save, xspace_save, ind_Y, Y, model, gx, normgx, neval, scale, CNTsin, gmat, hessian, hessian_bfgs, stage, msg
   

      #  Check if new function value is smaller
      #  and adjust indices if necessary
   
      if ( fX[m] < fX[i_xbest] ):    
         if ( verbose > 0 ):
            print ' i_xbest: '+str(i_xbest)
            print 'after moving best point to the bounds: '
            print ' i_xbest: '+str(m)
      
         xstatus[m] = c.inY()
         i_xbest    = copy(m)
         ind_Y[0]   = copy(m)
         fY[0]      = fX[m]
         x          = xnew
         fx         = fX[m]
      else:
         msg = 'current iterate moved to the bound but new fvalue not smaller'
         if ( verbose ):
            print msg
      
         return nit, i_xbest, m, X, fX, QZ, RZ, Delta, eps_current, cur_degree, sstatus, xstatus, dstatus, ndummyY, sspace_save, xspace_save, ind_Y, Y, model, gx, normgx, neval, scale, CNTsin, gmat, hessian, hessian_bfgs, stage, msg
         
   #  Return if all bounds active
   
   if ( len( I_LandU ) == n ):
      msg = 'all bounds active'
      
      if ( verbose ):
         print msg
   
      return nit, i_xbest, m, X, fX, QZ, RZ, Delta, eps_current, cur_degree, sstatus, xstatus, dstatus, ndummyY, sspace_save, xspace_save, ind_Y, Y, model, gx, normgx, neval, scale, CNTsin, gmat, hessian, hessian_bfgs, stage, msg
   
   #  There is at least one (nearly-) active bound,
   #  thus prepare for the subspace
   
   else:

      if ( logical_or( useMethod == 0, logical_and( useMethod == 2, stage > 1 ) ) ):  #  BCDFO or HYBRID

         #  Identify points which lie in the current subspace
         #  and identify points which lie very close to the
         #  current subspace (these are going to be dummy points)

         astatus = ones((m))
         for j in range(0,m):

            # exclude far points from the procedure

            if ( linalg.norm( x[:,0] - X[indfree,j] ) > Delta ):
               astatus[j] = 0
               sstatus[j] = c.out()

            # exclude dummy points from the procedure

            elif ( dstatus[j] == c.dummy() ):
               astatus[j] = 0
               sstatus[j] = c.out()

            else:

               for i in range(0,n):
                  ii = indfree[i]

                  if ( where(I_L == i)[0] >= 1 ):

                     # break if variable is NOT close to active bound
                     if ( abs(X[ii,j] - xl[0,ii]) > min(eps_bnd,pgx[0,i]) ):
                        astatus[j] = 0
                        sstatus[j] = c.out()
                        break

                     # not in subspace if variable is NOT exactly at bound
                     elif ( abs(X[ii,j] - xl[0,ii]) > 1e-10 ):
                        sstatus[j] = c.out()
                  

                  elif ( where(I_U == i)[0] >= 1 ):

                     # break if variable is NOT close to active bound
                     if ( abs(xu[0,ii] - X[ii,j]) > min(eps_bnd,pgx[0,i]) ):
                        astatus[j] = 0
                        sstatus[j] = c.out()
                        break

                     # not in subspace if variable is NOT exactly at bound
                     elif ( abs(xu[0,ii] - X[ii,j]) > 1e-10 ):
                        sstatus[j] = c.out()
                           
         #  Collect indices of points in the subspace and those close to a subspace

         ind_inSubspc = where( sstatus >= c.insubs() )[0] 
         ind_ctSubspc = where( logical_and( astatus == 1, sstatus == c.out() ) )[0]

         if ( verbose ):
           print 'nbr of points exact at a bound: '+str(len(ind_inSubspc))+' with indices: '+str(ind_inSubspc)
           print 'nbr of points close to a bound: '+str(len(ind_ctSubspc))+' with indices: '+str(ind_ctSubspc)
      
         #  Project nearby points onto their active bound(s) and get model-values
         #  of the projected points --> dummy points

         i_xbest_dummy = copy(i_xbest)
         for j in range(0,len(ind_ctSubspc)):

            k = ind_ctSubspc[j]
            y = X[indfree,k]
            
            for  i in range(0,n):
               ii = indfree[i]
               if ( vstatus[ii] == 1 ):
                  if ( abs( y[i] - xu_free[i] ) <= min(eps_bnd,pgx[0,i]) ):
                     y[i] = xu_free[i]
                  elif ( abs( y[i] - xl_free[i] ) <= min(eps_bnd,pgx[0,i]) ):
                     y[i] = xl_free[i]
                                    
            # check if dummy point coincides with a real point

            coincides = 0
            for i in range(0,len(ind_inSubspc)):
               ii = ind_inSubspc[i]
               if ( linalg.norm(y - X[indfree,ii]) < 1e-10 ):
                  coincides = ind_inSubspc[i]
                  if ( verbose ):
                     print 'dummy of point '+str(k)+' coincides with another point and is not taken into account'
                  break
                     
            if ( coincides == 0 ):

               mvalue = bcdfo_evalP( model, y, x, scale, shift_Y )

               # add dummy points to X

               m     = m + 1
               I     = eye(full_n)
               yfull = dot(I[:,indfix], xfix[indfix,0]) + dot(I[:,indfree], y)
               X     = append( X, yfull )
               fX    = append( fX, mvalue )

               # update status

               sstatus[m] = c.insubs()
               dstatus[m] = c.dummy()

               if ( verbose ):
                  print 'dummy of point '+str(k)+' is added as point '+str(m)+' to the set'            

               ind_inSubspc = append( ind_inSubspc, m )
               
      #else:  #  BFGS
      
         # no dummies in BFGS method

      # useMethod
    
      #  Merge all fixed variables and define new subspace
      
      I       = eye(full_n)
      xfix    = dot(I[:,indfix], xfix[indfix,0]) + dot(I[:,indfree], xbd[:,0])

      indfree = where( vstatus == c.free() )[0]
      indfix  = where( vstatus >= c.fixed() )[0]
      new_n   = length(indfree);

      #  Adjust the constants related to the dimension of the problem

      rep_degree_before = rep_degree
    
      if (rep_degree == plin):
         rep_degree = new_n + 1; 
      elif (rep_degree == pdiag):
         rep_degree = 2 * new_n + 1; 
      elif (rep_degree == pquad):
         rep_degree = ( (new_n + 1) * (new_n + 2) ) / 2; 
   
      plin  = new_n + 1;
      pdiag = 2 * new_n + 1;
      pquad = ( ( new_n + 1 ) * ( new_n + 2 ) ) / 2;

      #  Keep Hessian information for the return from the subspace

      if (logical_or( useMethod > 0, whichmodel == 7 )):

         hessian_before_ss = hessian
         hessian = hessian[ix_(i_free,i_free)]

         hessian_bfgs_before_ss = hessian_bfgs
         hessian_bfgs = hessian_bfgs[ix_(i_free,i_free)]
    
      #  Construct interpolation set in the subspace (n+1 points)
   
      if ( logical_or( useMethod == 0, logical_and( useMethod == 2, stage > 1 ) )):  #  BCDFO  or HYBRID

         if ( verbose ):
            print '****** choose new set of points in subspace nfree = '+str(new_n)+' ********'
      

         [ i_xbest, m, X, fX, cur_degree, ind_Y, Y, neval, xstatus, sstatus, dstatus,
            ndummyY, CNTsin, scale, gmat, msg ] = bcdfo_choose_lin( 
            f, i_xbest_dummy, xl, xu, X, fX, fxmax, Delta, kappa_th, 
            eps_L, plin, lSolver, vstatus, sstatus, dstatus, xfix, neval, maxeval, scale,
            shift_Y, c, 'subspace', plin, pdiag, pquad, whichmodel, stratLam, hardcons, 
            show_errg, normgx, kappa_ill, CNTsin, scaleX, scalefacX, gmat );

         if ( msg[0:5] == 'Error' ):
            return nit, i_xbest, m, X, fX, QZ, RZ, Delta, eps_current, cur_degree, sstatus, xstatus, dstatus, ndummyY, sspace_save, xspace_save, ind_Y, Y, model, gx, normgx, neval, scale, CNTsin, gmat, hessian, hessian_bfgs, stage, msg
      

         fY = fX[ind_Y]
         fx = copy(fX[i_xbest])
         x  = copy(X[indfree,i_xbest])

         #  Compute subspace-model and gradient

         [ QZ, RZ, x, scale ] = bcdfo_build_QR_of_Y( Y, whichmodel, shift_Y, 
            Delta, normgx, kappa_ill );

         subsmodel = zeros((1,pquad))

         model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, subsmodel, ind_Y, 
               0, 0, gx, scale, shift_Y, Delta, indfree, gmat, hessian_bfgs, 
               epsilon, noisy, stage );
         gx     = bcdfo_gradP( model, x, x, scale, shift_Y );
         normgx,pgp = bcdfo_projgrad( new_n, x, gx, xl[0,indfree], xu[0,indfree] );

         if ( logical_or( useMethod > 0, whichmodel == 7 ) ):
            g_bfgs = gmat[indfree,i_xbest]
         else:
            g_bfgs = []
      

         if ( verbose ):
            print '****** computed new model ********'
            print 'ind_Y = '
            print ind_Y      
      
      else:  #  BFGS

         cur_degree = 1
         x  = X[indfree,i_xbest]
         fx = fX[i_xbest]
         g_bfgs = gmat[indfree,i_xbest]
         gx = g_bfgs
         normgx,pgp = bcdfo_projgrad( new_n, x, gx, xl[0,indfree], xu[0,indfree] )
         ind_Y = copy(i_xbest)
         Y = X[indfree,i_xbest]

         if ( logical_and(useMethod == 2, stage == 1) ):
            [ QZ, RZ, x, scale ] = bcdfo_build_QR_of_Y( Y, whichmodel, shift_Y, 
            Delta, normgx, kappa_ill );

            subsmodel = zeros((1,pquad));
            model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, subsmodel, ind_Y, 
               0, 0, gx, scale, shift_Y, Delta, indfree, gmat, hessian_bfgs, 
               epsilon, noisy, stage );
            
      # useMethod

      #  Iteration printout after computing the subspace model

      if ( logical_and(logical_or( useMethod == 0, useMethod == 2 ), verbose_main >= 1) ):
         itype = 'subs'

         if ( show_errg ):

            #  Compute poisedness and gradient error, if unavailable so far.

            [ poised, Y_radius ] = bcdfo_poisedness_Y( QZ, RZ, Y, eps_L, x, 
            lSolver, whichmodel, hardcons, xl, xu, indfree, stratLam, scale, shift_Y );
            poisedness_known     = 1;

            errg = poised * Y_radius / factor_CV;

            #  Iteration summary

            print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % ( nit, neval, fx, normgx, errg, 0.0, Delta, 0.0, new_n, itype )

         else:

            print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % (nit, neval, fx, normgx, 0.0, Delta, 0.0, new_n, itype )
     
         # printout to file convhist.m

#         if ( verbose_main )
#            fid = fopen('convhist.m','a');
#            fprintf(fid,'#6d  #+.14e #.2e \n', neval, fx, normgx);
#            fclose(fid);
      
         if ( verbose_main == 3 ):
            bcdfo_print_summary_vector( 'x ', x )
            bcdfo_print_summary_vector( 'gx', gx )
         elif ( verbose_main > 3 ):
            bcdfo_print_vector( 'x ', x )
            bcdfo_print_vector( 'gx', gx )
            bcdfo_print_vector( 's', s )
      
         poisedness_known = 0;           # reset for the next iteration
   
         if ( verbose_main >= 2 ):
            #print ' cond(Z) = '+str( cond(RZ) )

            if ( verbose_main >= 10 ):
               bcdfo_print_vector( 'model', model );
                  
      if ( verbose ):
         print '****** solve subspace problem (nfree = '+str(new_n)+') ********'
         print '##################################################'
   
      #  Return if problem already converged (after evaluating f at possible dummy points)
    
      if ( normgx <= epsilon ) :
         if ( verbose ):
            print 'problem already converged'
            print ''       
      
         # if dummies in Y, recompute their real function values

         if ( ndummyY > 0 ):
            dummiesY = where( logical_and(dstatus == c.dummy(), xstatus == c.inY()) )[0]
         
            if (verbose):
               print 'dummies in Y - recompute fvalues. dummies: '+str(dummiesY)         
         
            for i in range(0,len( dummiesY )):

               k = dummiesY[i]
               j = where( k == ind_Y )[0]

               #  Update X and evaluate function

               if ( whichmodel >= 7 ):
                  compute_g = 1
               else:
                  compute_g = 0            
            
               [X, fX, neval, xstatus, sstatus, dstatus, gmat, msg] = bcdfo_augmX_evalf(
                  f, Y[:,j], k, X, fX, nfix, xfix, indfix, indfree,
                  fxmax, neval, maxeval, verbose, show_errg, xstatus, c.inY, sstatus,
                  dstatus, scaleX, scalefacX, compute_g, gmat)
               fY[j] = fX[k]
               
               if ( msg[0:5] == 'Error' ):
                  if ( level == 'toplevel' ):
                     print msg
               
                  return nit, i_xbest, m, X, fX, QZ, RZ, Delta, eps_current, cur_degree, sstatus, xstatus, dstatus, ndummyY, sspace_save, xspace_save, ind_Y, Y, model, gx, normgx, neval, scale, CNTsin, gmat, hessian, hessian_bfgs, stage, msg
                     
            ndummyY  = 0;

            #  Move to the best point found if different from x.

            [x, fx, QZ, RZ, Y, fY, ind_Y, i_xbest, scale] = bcdfo_find_smallf(
               c, QZ, RZ, Y, fY, ind_Y, i_xbest, cur_degree, 
               indfree, x, xl, xu, fx, dstatus, whichmodel, scale, shift_Y, 
               Delta, normgx, kappa_ill);

            #  Compute the associated polynomial interpolation model.
            #  (always n+1 points in Y)

            model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, model, 
                  ind_Y, 0, 0, gx, scale, shift_Y, Delta, indfree, gmat, 
                  hessian_bfgs, epsilon, noisy, stage );
            gx     = bcdfo_gradP( model, x, x, scale, shift_Y );
            normgx,pgp = bcdfo_projgrad(new_n,x,gx,xl[0,indfree],xu[0,indfree])

            if ( whichmodel >= 7 ):
               g_bfgs = gmat[indfree,i_xbest]
            else:
               g_bfgs = []
                            
      #  Solve problem in the subspace
  
      if ( normgx > epsilon ):

         [nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree,  
            model, gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, 
            sspace_save, xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, 
            stage, neval] = bcdfo_main(
         f, nit-1, nit-1, i_xbest, xl, xu, m, X, fX, ind_Y, QZ, RZ, Delta,
         cur_degree, neval, maxeval, maxit, model, gx, normgx, show_errg, pquad, 
         pdiag, plin, stallfact, eps_rho, Deltamax, rep_degree, epsilon, 
         verbose_main, eta1, eta2, gamma1, gamma2, gamma3, interpol_TR, factor_CV, 
         Lambda_XN, Lambda_CP, factor_FPU, factor_FPR, Lambda_FP, criterion_S, 
         criterion_FP, criterion_CP, mu, theta, eps_TR, eps_L, lSolver, stratLam, 
         eps_current, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save, 
         xspace_save, xfix, fxmax, poised_model, kappa_ill, kappa_th, eps_bnd, 
         poised, Y_radius, c, 'sublevel', whichmodel, hardcons, noisy, scaleX, 
         scalefacX, CNTsin, shrink_Delta, gmat, hessian, hessian_bfgs, g_bfgs, useMethod, 
         stage, Deltamin, scale, shift_Y);
      
      if ( verbose ):
         print '****** exit subspace problem (nfree = '+str(new_n)+') ********'
         print '##################################################'  

      if ( msg[0:5] == 'Error' ):
         #  in case of an error, re-assemble gradient
         i_fix = setdiff(range(0,n),i_free)
         I     = eye(n)
         nfix  = n - len( where( indfree > 0 )[0])
         gx    = dot( I[:,i_fix], zeros((nfix,1))) + dot(I[:,i_free], gx)
         return nit, i_xbest, m, X, fX, QZ, RZ, Delta, eps_current, cur_degree, sstatus, xstatus, dstatus, ndummyY, sspace_save, xspace_save, ind_Y, Y, model, gx, normgx, neval, scale, CNTsin, gmat, hessian, hessian_bfgs, stage, msg
       
      #  Store index of best point when exiting all subspaces and going to the toplevel
      #  (to avoid reentering the same subspace later on)
	
      nbr_ss_last = size(sspace_save,2)

      for j in range(nbr_ss, nbr_ss_last):
         xspace_save[1,j] = i_xbest  

   #  Compute full-dimensional model and gradient if in the toplevel,
   #  otherwise return to next upper level

   if ( level == 'toplevel' ):

      #  fix only alwaysfixed variables, others set to free
    
      indfree          = where( vstatus < c.alwaysfixed() )[0]
      vstatus[indfree] = c.free()
      sstatus          = ones((m))
      n                = len(indfree)

      #  Reset some constants

      rep_degree = rep_degree_before

      plin  = n + 1
      pdiag = 2 * n + 1
      pquad = ( ( n + 1 ) * ( n + 2 ) ) / 2

      #  Re-assemble Hessian from full space and subspace

      if ( logical_or( useMethod > 0, whichmodel == 7 ) ):

         hesse   = hessian_before_ss
         hesse[ix_(i_free,i_free)] = hessian
         hessian = hesse

         h_bfgs  = hessian_bfgs_before_ss
         h_bfgs[ix_(i_free,i_free)] = hessian_bfgs
         hessian_bfgs = h_bfgs
      
      if ( logical_or( useMethod == 0, logical_and( useMethod == 2, stage > 1 ) ) ):  #  BCDFO or HYBRID

         if ( verbose ):
            print '****** choose Y in the full space ********'
     
         Delta = epsilon
         [ i_xbest, m, X, fX, cur_degree, ind_Y, Y, neval, xstatus, sstatus, dstatus,
             ndummyY, CNTsin, scale, gmat, msg ] = bcdfo_choose_lin( 
             f, i_xbest, xl, xu, X, fX, fxmax, Delta, kappa_th, eps_L, 
            rep_degree, lSolver, vstatus, sstatus, dstatus, xfix, neval, maxeval, scale, 
            shift_Y, c, 'fullspace', plin, pdiag, pquad, whichmodel, stratLam, hardcons,
            show_errg, normgx, kappa_ill, CNTsin, scaleX, scalefacX, gmat )

         if ( msg[0:5] == 'Error' ):
            return nit, i_xbest, m, X, fX, QZ, RZ, Delta, eps_current, cur_degree, sstatus, xstatus, dstatus, ndummyY, sspace_save, xspace_save, ind_Y, Y, model, gx, normgx, neval, scale, CNTsin, gmat, hessian, hessian_bfgs, stage, msg
      
         fY = fX[ind_Y]
         fx = fX[i_xbest]
         x  = X[indfree,i_xbest]

         #  Compute QR factorization of new coefficient matrix

         [ QZ, RZ, x, scale ] = bcdfo_build_QR_of_Y( Y, whichmodel, shift_Y, 
            Delta, normgx, kappa_ill )

         #  Compute model and gradient

         if ( verbose ):
            disp('****** compute full-dimensional model ********')
      
         restartmodel = zeros((1,pquad))
         model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, restartmodel, 
               ind_Y, 0, 0, gx, scale, shift_Y, Delta, indfree, gmat, hessian_bfgs, 
               epsilon, noisy, stage )
         gx     = bcdfo_gradP( model, x, x, scale, shift_Y )
         normgx,pgp = bcdfo_projgrad( n, x, gx, xl[0,indfree], xu[0,indfree] )

         if ( logical_or( useMethod > 0, whichmodel == 7 ) ):
            g_bfgs = gmat[indfree,i_xbest]
         else:
            g_bfgs = []
      
      
      else:  #  BFGS
   
         cur_degree = 1
         fY      = fX[ind_Y]
         fx      = fX[i_xbest]
         x       = X[indfree,i_xbest]
         g_bfgs  = gmat[indfree,i_xbest]
         gx      = g_bfgs
         normgx,pgp  = bcdfo_projgrad( n, x, gx, xl[0,indfree], xu[0,indfree] )
         ind_Y   = copy(i_xbest)
         Y       = X[indfree,i_xbest]

         if ( logical_and( useMethod == 2, stage == 1 ) ):
            [ QZ, RZ, x, scale ] = bcdfo_build_QR_of_Y( Y, whichmodel, shift_Y, 
            Delta, normgx, kappa_ill );

            restartmodel = zeros((1,pquad))
            model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, restartmodel, 
               ind_Y, 0, 0, gx, scale, shift_Y, Delta, indfree, gmat, hessian_bfgs, 
               epsilon, noisy, stage )
            
       # end useMethod

      #  Iteration output after exiting a subspace

      if ( logical_and( logical_or( useMethod == 0, useMethod == 2 ), verbose_main >= 1 ) ):
         itype = 'back'
   
         if ( show_errg ):

            #  Compute poisedness and gradient error, if unavailable so far.

            [ poised, Y_radius ] = bcdfo_poisedness_Y( QZ, RZ, Y, eps_L, x, 
               lSolver, whichmodel, hardcons, xl, xu, indfree, stratLam, scale, shift_Y )
            poisedness_known     = 1
            errg = poised * Y_radius / factor_CV

            #  Iteration summary

            print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % (nit, neval, fx, normgx, errg, 0.0, Delta, 0.0, n, itype)

         else:

            print '%5d  %5d  %+.14e  %.2e  %.2e  %.2e  %+.2e  %3d  %4s' % (nit, neval, fx, normgx, 0.0, Delta, 0.0, n, itype )
      
         # printout to file convhist.m

#         if ( verbose_main )
#            fid = fopen('convhist.m','a');
#            fprintf(fid,'#6d  #+.14e #.2e \n', neval, fx, normgx);
#            fclose(fid);
      
         if ( verbose_main == 3 ):
            bcdfo_print_summary_vector( 'x ', x )
            bcdfo_print_summary_vector( 'gx', gx )
         elif ( verbose_main > 3 ):
            bcdfo_print_vector( 'x ', x )
            bcdfo_print_vector( 'gx', gx )
            bcdfo_print_vector( 's', s )
      
         poisedness_known = 0;           # reset for the next iteration

         if ( verbose_main >= 2 ):
            #print ' cond(Z) = '+str( linalg.cond(RZ) )

            if ( verbose_main >= 10 ):
               bcdfo_print_vector( 'model', model );
                    
      #  Adjust Delta if not converged in the full-space

      if ( logical_and( normgx > epsilon, strcmp(level,'toplevel') ) ):
         Delta = max( 0.1 * normgx, epsilon )   

      if ( verbose ):
         normgx
         cur_degree
         Delta
   
   else:

      if ( logical_or(useMethod > 0, whichmodel == 7) ):  # BFGS and HYBRID

         #  Re-assemble the Hessian coming from each subspace to the next higher one

         hesse   = hessian_before_ss;
         hesse[ix_(i_free,i_free)] = hessian;
         hessian = hesse;

         h_bfgs  = hessian_bfgs_before_ss;
         h_bfgs[ix_(i_free,i_free)] = hessian_bfgs;
         hessian_bfgs = h_bfgs;

   if ( verbose ):
      print '****** return from ident_active_bounds **********'

