#!/usr/local/bin/python
from numpy import *
from bcdfo_solve_TR_MS import *
from helper import generalDecorator
from runtime import matlabarray, char

class bcdfo_solve_TR_MS_bc_decorator(generalDecorator):
	def __init__(self, f):
		self.f = f

	def __call__(self, gx, H, lb, ub, Delta, eps_D, stratLam, nargout=None):
		
		gx = self.convert(gx)
		#print "gx", gx
		#print "type(gx)", type(gx)
		H = self.convert(H)
		#print "H", H
		#print "type(H)", type(H)
		lb = self.convert(lb)
		#print "lb", lb
		#print "type(lb)", type(lb)
		ub = self.convert(ub)
		#print "ub", ub
		#print "type(ub)", type(ub)
		Delta = self.convert(Delta)
		#print "Delta", Delta
		#print "type(Delta)", type(Delta)
		eps_D = self.convert(eps_D)
		#print "eps_D", eps_D
		#print "type(eps_D)", type(eps_D)
		stratLam = self.convert(stratLam)
		#print "stratLam", stratLam
		#print "type(stratLam)", type(stratLam)
		
		
		s, lamb, norms, value, gplus, nfact, neigd, msg = self.f(gx, H, lb, ub, Delta, eps_D, stratLam)
		#self.printTypes([s, lamb, norms, value, gplus, nfact, neigd, msg])

		
		s = matlabarray(s).T
		lamb = matlabarray(lamb)
		norms = matlabarray(norms)
		value = matlabarray(value)
		gplus = matlabarray(gplus)
		nfact = matlabarray(nfact)
		neigd = matlabarray(neigd)
		msg = char(msg)
				
		return s, lamb, norms, value, gplus, nfact, neigd, msg
		
@bcdfo_solve_TR_MS_bc_decorator
def bcdfo_solve_TR_MS_bc( gx, H, lb, ub, Delta, eps_D, stratLam, ):
###############################################################################
#
#  An implementation of exact trust-region minimization based on the
#  More-Sorensen algorithm subject to bound constraints.
#
#  INPUT: 
#
#  g        : the model's gradient
#  H        : the model's Hessian
#  lb       : lower bounds on the step
#  ub       : upper bounds on the step
#  Delta    : the trust-region's radius
#  eps_D    : the accuracy required on the equation ||s|| = Delta for a
#             boundary solution
#  stratLam : the strategy to adjust lamb to find an active bound
#             (1 - Newtons method, 0 - bisection method)
#
#  OUTPUT:
#
#  s        : the trust-region step
#  lamb   : the Lagrange multiplier corresponding to the trust-region constraint
#  norms    : the norm of s
#  value    : the value of the model at the optimal solution
#  gplus    : the value of the model's gradient at the optimal solution
#  nfact    : the number of Cholesky factorization required
#  neigd    : the number of eifenvalue decompositions required
#  msg      : an information message
#
#  DEPENDENCIES: bcdfo_solve_TR_MS
#
#  PROGRAMMING: A. Troeltzsch, S. Gratton, July 2009. 
#              ( This version 14 I 2010 )
#
#  TEST:
#
#  bcdfo_solve_TR_MS_bc( [2; 3], [4 6; 6 5], [-10; -10], [10; 10], 1.0, 0.001, 1 )
#  should give
#    0.5153
#   -0.8575
#
#  bcdfo_solve_TR_MS_bc( [2; 3], [4 6; 6 5], [-0.1; -0.1], [10; 10], 1.0, 0.001, 1 )
#  should give
#   -0.1
#   -0.1
#
#  bcdfo_solve_TR_MS_bc( [2; 3], [4 6; 6 5], [-10; -10], [0; 0], 1.0, 0.001, 1 )
#  should give
#    0
#   -0.6
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
###############################################################################

   verbose    = 0;
   theta      = 1.0e-13;          # accuracy of the interval on lamb
   eps_bound  = 1.0e-5;           # the max distance | bound - x | < eps_bound for 
                                  # a boundary solution
                                    
   # initialization 
   
   msg        = '';
   lamb       = 0;                # initial lamb
   lamb_save  = []
   norms_b_save = []
   value      = 0;                # initial model value
   nfact      = 0;                # initialize in case of error
   neigd      = 0;                # initialize in case of error
   Delta0     = Delta;
   g          = copy(gx)
   g0         = copy(gx)                # copy initial gradient
   gplus      = copy(gx)                # initialize gplus
   s          = zeros_like(g)     # initial zero step
   norms      = 0;
   n          = len(g[0]);        # space dimension
   I          = eye(n);           # identity matrix used for projection
   ind_active = array([],int)     # indices of active bounds
   ind_free   = arange(0,n)        # indices of inactive bounds
   nfree      = n;                # nbr of inactive indeces

   if ( isnan(H).any() ):
      if (verbose ):
         disp( 'Error in bcdfo_solve_TR_MS_bc: H contains NaNs!' )
      
      msg = 'Error1';
      return s, lamb, norms, value, gplus, nfact, neigd, msg
   
   if ( not isreal(H).all() ):
      if (verbose ):
         disp( 'Error in bcdfo_solve_TR_MS_bc: H contains imaginary parts!' )
      
      msg = 'Error2';
      return s, lamb, norms, value, gplus, nfact, neigd, msg
   
   if ( isinf(H).any() ):
      if (verbose ):
         disp( 'Error in bcdfo_solve_TR_MS_bc: H contains infinite elements!' )
      
      msg = 'Error3';
      return s, lamb, norms, value, gplus, nfact, neigd, msg
   
   #  Fix active components
   
   ind_g_crit = where( logical_or( logical_and( g[0] > 0, abs(lb[0]) <= 1e-10 ), logical_and( g[0] < 0, ub[0] <= 1e-10 )))[0]
   
   if ( len(ind_g_crit) != 0):
       ind_active = ind_free[ind_g_crit]
       ind_free   = list(set(ind_free)-set(ind_active) )
       nfree      = len(ind_free)

   #  Loop until no free variables anymore
   
   j=0;
   while nfree > 0:
   
      #  Loop until all active variables are detected and fixed
       
      new_call_to_MS = 1
      while ( new_call_to_MS == 1 ):
         j=j+1
            
         #  Minimize system in the (possibly) reduced subspace
   
         if ( verbose >= 1 ):
            print '('+str(j)+') ---- minimizing in the (sub)space of '+str(len(ind_free))+' variable(s)'

         g_reduced = g[:,ind_free]
         
         H_reduced = H[ix_(ind_free,ind_free)]

         #  Call unconstrained MS in (possibly) reduced subspace
         
         [ s_deltaMS, lamb, norms_deltaMS, value_red, gplus_red, nfact_r, neigd_r, msg, hardcase ] = bcdfo_solve_TR_MS(
            g_reduced, H_reduced, Delta, eps_D )
         

         nfact = nfact + nfact_r
         neigd = neigd + neigd_r
        
         gplus[0,ind_free] = gplus[0,ind_free] + gplus_red[0]

         
         s_after_reduced_ms = s + dot( s_deltaMS, I[ind_free,:] )
         
         #  Compute critical components which became active during the last MS iteration
           
         ind_u_crit = where( logical_and(ub[0,ind_free]-s_after_reduced_ms[0,ind_free] <= eps_bound, ub[0,ind_free] <= 1e-10) )[0]
         ind_l_crit = where( logical_and(s_after_reduced_ms[0,ind_free]-lb[0,ind_free] <= eps_bound, lb[0,ind_free] >= -1e-10) )[0]
           
         #  Fix these new active components
           
         if ( len(ind_u_crit) + len(ind_l_crit) != 0 ):
              
            ind_active = append( ind_active, ind_free[ind_u_crit] )
            ind_active = append( ind_active, ind_free[ind_l_crit] )
            ind_free   = list( set(range(0,n))-set(ind_active) )
            nfree      = len(ind_free)
            
            if ( verbose ):
               disp('fixed one or more variables')            
               
            #  If no inactive variables anymore --> exit
               
            if ( nfree == 0 ):
               norms = linalg.norm( s )
               value = 0.5 * dot( dot( s, H ), s.transpose() ) + dot( s, g0.transpose() )
               if ( verbose ):
                  disp('no inactive variables anymore - return')
               
               return s, lamb, norms, value, gplus, nfact, neigd, msg
                           
         else:            
            new_call_to_MS = 0;
                      
      #  Check if step is outside the bounds
       
      if ( verbose == 2 ):
         disp('check if step inside bounds')      
   
      out_of_ubound = where( (ub[0,ind_free]-s_after_reduced_ms[0,ind_free]) < 0.0 )[0]
      out_of_lbound = where( (s_after_reduced_ms[0,ind_free]-lb[0,ind_free]) < 0.0 )[0]
      out_of_ubound_init = out_of_ubound
      out_of_lbound_init = out_of_lbound
   
      if ( len(out_of_ubound) + len(out_of_lbound) != 0 ):
   
         back_inside = 0;
         lamb0 = lamb;
         if ( verbose == 2 ):
            disp( 'step outside bounds!' )
            print 'lambda_0='+str(lamb0)         
   
         #  Set lower bound on lamb.
   
         lower = lamb
           
         #  New lamb for bisection method
          
         if ( stratLam == 0 ):
            lamb = max( 2.0, 2 * lamb )
         
           
         #  Compute upper bound on lamb (using the closest bound out of the hit bounds) 
           
         gnorm = linalg.norm( g )
         if ( len(out_of_ubound) > 0 ):
            delta_b = min( abs( ub[0,ind_free[out_of_ubound]]-s[0,ind_free[out_of_ubound]] ) )
         
         if ( len(out_of_lbound) > 0 ):
            delta_b = min( abs( lb[0,ind_free[out_of_lbound]]-s[0,ind_free[out_of_lbound]] ) )
            if ( len(out_of_ubound) > 0 ):
               delta_b  = min( min( abs( ub[0,ind_free[out_of_ubound]]-s[0,ind_free[out_of_ubound]] ) ), delta_b )
         
         goverD   = gnorm / delta_b
         Hnorminf = linalg.norm( H, inf )
         if ( Hnorminf > 0 ):        # because Octave generates a NaN for null matrices.
            HnormF   = linalg.norm( H, 'fro' )
         else:
            HnormF   = 0
         
         upper    = max( 0, goverD + min( Hnorminf, HnormF ) );
   
         #  Compute active components
   
         ind_u_active = where( abs( ub[0,ind_free]-s_after_reduced_ms[0,ind_free] ) <= eps_bound )[0]
         ind_l_active = where( abs( s_after_reduced_ms[0,ind_free]-lb[0,ind_free] ) <= eps_bound )[0]
   
         #  Loop on successive trial values for lamb.
   
         i = 0
         while ((( len(ind_u_active) + len(ind_l_active)) == 0) or (len(out_of_lbound) + len(out_of_ubound) != 0 )):
            i = i + 1
   
            #  Print lambda value
   
            old_lamb = lamb
            new_lamb = -1
            if ( verbose ):
               print 'solve_TR_MS_bc ('+str(i)+'): lower = '+str(lower)+' lamb = '+str(lamb)+' upper = '+str(upper)            
   
            #  Factorize H + lamb * I.
            
            try:
               R = linalg.cholesky( H[ix_(ind_free,ind_free)] + lamb * eye( nfree ) ).T
               p = 0
            except:
               R = array([[]])
               p = 1   
  
            if ( isnan( R.any() ) ):
               disp( 'Error in bcdfo_solve_TR_MS_bc: NaNs in Cholesky factorization' )
               msg = 'Error4';
               return s, lamb, norms, value, gplus, nfact, neigd, msg
            
            nfact    = nfact + 1;
   
            #  Successful factorization 
   
            if ( p == 0 and hardcase == 0 ):
               
               s_deltaH = - linalg.solve( R, linalg.solve( R.transpose(), g[:,ind_free].transpose() ) ).T
               
               s_duringH = s + dot( s_deltaH,I[:,ind_free] )
   
               #  Find components which are at its bound and became active
   
               ind_u_crit  = where( logical_and(ub[0,ind_free]-s_duringH[0,ind_free] <= eps_bound, ub[0,ind_free] <= 1e-10 ) )[0]
               ind_l_crit  = where( logical_and(s_duringH[0,ind_free]-lb[0,ind_free] <= eps_bound, lb[0,ind_free] >= -1e-10 ) )[0]
   
               #  Set these new active components to zero for one iteration
   
               if ( len(ind_u_crit) != 0 ):
                  s_deltaH[0,ind_u_crit]  = 0.0
                  s_duringH[0,ind_free[ind_u_crit]] = 0.0
               
               if ( len(ind_l_crit) != 0):
                  s_deltaH[0,ind_l_crit]  = 0.0
                  s_duringH[0,ind_free[ind_l_crit]] = 0.0
               
   
               out_of_ubound = nonzero( (ub[0,ind_free]-s_duringH[0,ind_free]) < 0.0 )[0]
               out_of_lbound = nonzero( (s_duringH[0,ind_free]-lb[0,ind_free]) < 0.0 )[0]
               
               # find an appropriate bound for the next homotopy-step when using Newton's method
   
               if ( stratLam == 1 or verbose > 0 ):
   
                  if ( logical_or( len(out_of_ubound) != 0, len(out_of_lbound) != 0 )):
   
                     #  OUTSIDE the bounds: find the furthest step component 
                     
                     outside = 1    
                                                                 
                     if ( len(out_of_ubound) != 0 ):
                        diff_b_u = max( abs( ub[0,ind_free[out_of_ubound]] - s_duringH[0,ind_free[out_of_ubound]] ) )
                        ind_b_u  = abs( ub[0,ind_free[out_of_ubound]] - s_duringH[0,ind_free[out_of_ubound]] ).argmax()
                        norms_b = abs( s_deltaH[0,out_of_ubound[ind_b_u]] )
                        delta_b = abs( ub[0,ind_free[out_of_ubound[ind_b_u]]] - s[0,ind_free[out_of_ubound[ind_b_u]]] )
                        ind_b   = out_of_ubound[ind_b_u]
                        sign_b  = sign( ub[0,ind_free[out_of_ubound[ind_b_u]]] - s[0,ind_free[out_of_ubound[ind_b_u]]] )
                        out_of_ubound_init = append(out_of_ubound_init, out_of_ubound )
                     
                     if ( len(out_of_lbound) != 0 ):
                        diff_b_l = max( abs( s_duringH[0,ind_free[out_of_lbound]] - lb[0,ind_free[out_of_lbound]] ) )
                        ind_b_l  = abs( s_duringH[0,ind_free[out_of_lbound]] - lb[0,ind_free[out_of_lbound]] ).argmax()
                        norms_b = abs( s_deltaH[0,out_of_lbound[ind_b_l]] )
                        delta_b = abs( lb[0,ind_free[out_of_lbound[ind_b_l]]] - s[0,ind_free[out_of_lbound[ind_b_l]]] )
                        ind_b   = out_of_lbound[ind_b_l]
                        sign_b  = sign( lb[0,ind_free[out_of_lbound[ind_b_l]]] - s[0,ind_free[out_of_lbound[ind_b_l]]] )
                        out_of_lbound_init = append(out_of_lbound, out_of_lbound)
                        
                     if ( len(out_of_ubound) != 0 and len(out_of_lbound) != 0 ):
                        if (diff_b_u > diff_b_l):
                           norms_b = abs( s_deltaH[0,out_of_ubound[ind_b_u]] )
                           delta_b = abs( ub[0,ind_free[out_of_ubound[ind_b_u]]] - s[0,ind_free[out_of_ubound[ind_b_u]]] )
                           ind_b   = out_of_ubound[ind_b_u]
                           sign_b  = sign( ub[0,ind_free[out_of_ubound[ind_b_u]]] - s[0,ind_free[out_of_ubound[ind_b_u]]] )
                        else:
                           norms_b = abs( s_deltaH[0,out_of_lbound[ind_b_l]] )
                           delta_b = abs( lb[0,ind_free[out_of_lbound[ind_b_l]]] - s[0,ind_free[out_of_lbound[ind_b_l]]] )
                           ind_b   = out_of_lbound[ind_b_l]
                           sign_b  = sign( lb[0,ind_free[out_of_lbound[ind_b_l]]] - s[0,ind_free[out_of_lbound[ind_b_l]]] )
   
                  else:
                  
                     #  INSIDE the bounds but no step component active:
                     #  find the closest components to its bound from the
                     #  set of components which were initially outside
   
                     outside = 0
                      
                     if ( len(out_of_ubound_init) != 0 ):
                        diff_b_u = min( abs( ub[0,ind_free[out_of_ubound_init]] - s_duringH[0,ind_free[out_of_ubound_init]] ) );
                        ind_b_u = abs( ub[0,ind_free[out_of_ubound_init]] - s_duringH[0,ind_free[out_of_ubound_init]] ).argmin()
                        norms_b = abs( s_deltaH[0,out_of_ubound_init[ind_b_u]] );
                        delta_b = abs( ub[0,ind_free[out_of_ubound_init[ind_b_u]]] - s[0,ind_free[out_of_ubound_init[ind_b_u]]])
                        ind_b   = out_of_ubound_init[ind_b_u]
                        sign_b  = sign( ub[0,ind_free[out_of_ubound_init[ind_b_u]]] - s[0,ind_free[out_of_ubound_init[ind_b_u]]] )
                        
                     if ( len(out_of_lbound_init) != 0 ):
                        diff_b_l = min( abs( s_duringH[0,ind_free[out_of_lbound_init]] - lb[0,ind_free[out_of_lbound_init]] ) )
                        ind_b_l = abs( s_duringH[0,ind_free[out_of_lbound_init]] - lb[0,ind_free[out_of_lbound_init]] ).argmin()
                        norms_b = abs( s_deltaH[0,out_of_lbound_init[ind_b_l]] )
                        delta_b = abs( lb[0,ind_free[out_of_lbound_init[ind_b_l]]] - s[0,ind_free[out_of_lbound_init[ind_b_l]]] )
                        ind_b   = out_of_lbound_init[ind_b_l]
                        sign_b  = sign( lb[0,ind_free[out_of_lbound_init[ind_b_l]]] - s[0,ind_free[out_of_lbound_init[ind_b_l]]] )
                          
                     if ( len(out_of_ubound_init) != 0 and len(out_of_lbound_init) != 0 ):
                        if ( diff_b_u < diff_b_l ):
                           norms_b = abs( s_deltaH[0,out_of_ubound_init[ind_b_u]] )
                           delta_b = abs( ub[0,ind_free[out_of_ubound_init[ind_b_u]]] - s[0,ind_free[out_of_ubound_init[ind_b_u]]] )
                           ind_b   = out_of_ubound_init[ind_b_u]
                           sign_b  = sign( ub[0,ind_free[out_of_ubound_init[ind_b_u]]] - s[0,ind_free[out_of_ubound_init[ind_b_u]]] )
                           
                        else:
                           norms_b = abs( s_deltaH[0,out_of_lbound_init[ind_b_l]] )
                           delta_b = abs( lb[0,ind_free[out_of_lbound_init[ind_b_l]]] - s[0,ind_free[out_of_lbound_init[ind_b_l]]] )
                           ind_b   = out_of_lbound_init[ind_b_l]
                           sign_b  = sign( lb[0,ind_free[out_of_lbound_init[ind_b_l]]] - s[0,ind_free[out_of_lbound_init[ind_b_l]]] )     
                                                      
                  
               #  Iteration printout
   
               if ( verbose ):
                  lamb_save.append(lamb)
                  norms_b_save.append(norms_b)
                  
                
                  if ( outside == 0 ):
                     print '%s%d%s %12.8e %s %12.8e %s' % ('solve_TR_MS_bc (', i,'): |s_i| = ', 
                     norms_b, ' |bound_i| = ', delta_b, ' s < bounds')
                  else:
                     print '%s%d%s %12.8e %s %12.8e' % ('solve_TR_MS_bc (', i,'): |s_i| = ', 
                     norms_b, ' |bound_i| = ', delta_b)
                  
               
   
               #  Test if step inside bounds +/- eps_bound
                    
               out_of_uEpsbound = where( (ub[0,ind_free] - s_duringH[0,ind_free]) < -eps_bound )[0]
               out_of_lEpsbound = where( (s_duringH[0,ind_free] - lb[0,ind_free]) < -eps_bound )[0]
   
               if ( len(out_of_uEpsbound) + len(out_of_lEpsbound) == 0 ):
                       
                  if ( verbose >= 2 ):
                     disp('all components inside the bounds + eps_bound')
                                          
                  back_inside = 1
   
                  # check if at least one component active
   
                  ind_u_active = where( abs( ub[0,ind_free] - s_duringH[0,ind_free] ) <= eps_bound )[0]
                  ind_l_active = where( abs( s_duringH[0,ind_free] - lb[0,ind_free] ) <= eps_bound )[0]
   
                  if ( (len(ind_u_active) + len(ind_l_active)) != 0 ):
                     if ( verbose >= 2 ):
                        print 'all components inside the bounds + eps_bound, '+str(len(ind_u_active)+len(ind_l_active))+' component/s close to one of its bounds'
                     
   
                     #  Compute the current step after the homotopy-step
   
                     s_afterH = s + dot( s_deltaH, I[:,ind_free] )
   
                     #  Move active components to their bounds
                           
                     if ( len(ind_u_active) > 0 ):
                        s_afterH[0,ind_free[ind_u_active]] = ub[0,ind_free[ind_u_active]]
                     
                     if ( len(ind_l_active) > 0 ):
                        s_afterH[0,ind_free[ind_l_active]] = lb[0,ind_free[ind_l_active]]
                     
   
                     #  Define information message.
                           
                     msg = 'boundary solution';
   
                     break;
                  
               
   
               #  Compute new lamb
                    
               if ( stratLam == 0 ):		#  Using bisection method
     
                  if ( back_inside == 0 ):
                     lamb = 2 * lamb
                     if ( upper < lamb ):
                        upper = 2*lamb; 
                     
                  else:
                     if ( not out_of_ubound and not out_of_lbound ):
                        upper = lamb;
                     else:
                        lower = lamb;
                     
                     new_lamb  = (lamb + old_lamb) / 2;
                     theta_range = theta * ( upper - lower );
                     if ( logical_and( new_lamb > lower + theta_range, new_lamb < upper - theta_range) ):
                        lamb = new_lamb;
                     else:
                        lamb = max( sqrt( lower * upper ), lower + theta_range );
                     
                  
                        
               else:   #  Using  Newton's iteration on the modified secular equation
   
                  #  Reset bounds on lamb
   
                  if ( len(out_of_ubound ) + len(out_of_lbound) == 0 ):
                     upper = lamb;
                     
                  else:
                     lower = lamb;
                  
                  
                       
                  #  check the sign of the chosen step component
                       
                  unitvec        = zeros((1,nfree));
                  unitvec[0,ind_b] = 1;
                  es             = dot( unitvec, s_deltaH.transpose() )[0]
                  
                  if ( sign(es) != sign_b ):
                           
                     #  if step component has the wrong sign
                     #  (other than the active bound): one bisection step
                           
                     new_lamb  = (lower + upper) / 2;
                     
                  else:
                  
   	               #  Newton step
   
                     w1 = linalg.solve( R.transpose(), unitvec.transpose())
                     w2 = linalg.solve( R.transpose(), s_deltaH.transpose())
                    
                     new_lamb = lamb + ( ( norms_b - delta_b ) / delta_b ) * ( norms_b**2 / ( es * dot(w1.transpose(), w2)[0] ) )[0]
                     if ( back_inside == 0 and upper <= new_lamb ):
                        upper = 2 * new_lamb;
                     
                  
   
                  #  Check new value of lamb wrt its bounds
   
                  theta_range = theta * ( upper - lower );
                  if ( new_lamb > lower + theta_range and new_lamb <= upper - theta_range ):
                     lamb = new_lamb;
                  else:
                     lamb = max( sqrt( lower * upper ), lower + theta_range );
                  
   
               
   
            else: #  Unsuccessful factorization: compute new lamb 
                   
               if ( verbose ):
                  disp('unsuccessful factorization')
               
               hardcase = 0;
               lower  = lamb;
               t      = 0.5;
               lamb = ( 1 - t ) * lower + t * upper;
            
   
            #  Return with error message after 100 iterations
   
            if ( i >= 10 ):
               s      = zeros((1,n))
               norms  = 0;
               msg    = 'Error6 in bcdfo_solve_TR_MS_bc: iteration limit in bc-MS exceeded!'
               print msg
               return s, lamb, norms, value, gplus, nfact, neigd, msg
            
   
          # end of while-loop
               
      else:
           
         # MS was successful inside the bounds
           
         if ( verbose >= 2 ):
            disp( 'step inside bounds!' ) 
         
           
         #  Define information message.
   
         msg = '(partly) interior solution';
                  
         #  update the step and model value

         s     = s_after_reduced_ms;
         norms = linalg.norm( s );         
         value = 0.5 * dot(dot(s, H), s.transpose()) + dot(s, g0.transpose())

         return s, lamb, norms, value, gplus, nfact, neigd, msg
          
      #  update the step and model value
   
      s     = s_afterH;
      norms = linalg.norm( s );   
      value = 0.5 * dot(dot(s, H), s.transpose()) + dot(s, g0.transpose())
       
      #  update trust-region radius
   
      Delta = Delta0 - norms
       
      if ( Delta < -eps_bound ):
         disp('Error in bcdfo_solve_TR_MS_bc: delta smaller than zero !!!!!!')
         msg = 'Error7'; 
         return s, lamb, norms, value, gplus, nfact, neigd, msg
      elif ( Delta < 0 ):
         return s, lamb, norms, value, gplus, nfact, neigd, msg
      
       
      #  update gradient 
      
      g     = g0 + dot( s,H )
       
      #  update active bounds
       
      ind_active = where( logical_or( (ub[0]-s[0]) <= eps_bound, (s[0]-lb[0]) <= eps_bound ) )[0]
      ind_free   = list( set(range(0,n)) - set(ind_active) )
      nfree      = len(ind_free)
       
      if ( nfree > 0 ):
          
         #  check first-order criticality
          
         ng_reduced = linalg.norm( g[0,ind_free], inf )
         if ( ng_reduced <= 1e-5 ):
            if ( verbose >= 2 ):
               disp('point first order critical - return')
               ng_reduced
            
            return s, lamb, norms, value, gplus, nfact, neigd, msg
         
       
         #  Fix components which became active
       
         ind_g_crit = where( logical_or( logical_and( abs( lb[0,ind_free] ) <= 1e-10, g[0,ind_free] > 0 ), logical_and( ub[0,ind_free] <= 1e-10, g[0,ind_free] < 0 )))[0]
   
         if ( len(ind_g_crit) != 0):
            ind_active = append( ind_active, ind_free[ind_g_crit] )
            ind_free   = list( set(ind_free) - set(ind_active) )
            nfree      = len( ind_free )
         
      
   
       
             
