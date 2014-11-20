#! /usr/bin/env python
from numpy import *

def bcdfo_projected_tcg( h, gg, H, lb, ub, cgitmax, epsi, epsi_max, verbose ):
################################################################################
# A very simple (bound-)projected truncated conjugate-gradient algorithm,
# without preconditioning
#
# Ph. Toint, 9 December 2005.
#
# verbose: output amount: 0 (silent), 1, 2, 3 (debug)
################################################################################

   # Initializations

   verbose   = 0
   value     = 0
   g         = copy(gg)
   n         = len( gg[0] )
   cgitstart = 1
   s_save    = array([[]])
   g_save    = array([[]])
   v_save    = array([])
   eps       = finfo(float).eps

   #  Define the active bounds.

   xfree = ones(( 1, n ))
   nfree = n
   for i in range(0,n):
      if ( logical_or(logical_and( abs( lb[0,i] ) <= eps, g[0,i] > 0 ), logical_and( abs(ub[0,i]) <= eps, g[0,i] < 0 ) ) ):
         xfree[0,i] = 0
         nfree      = nfree - 1

   ng = linalg.norm( g * xfree )

   #  Define TCG adapt

   eps_term = max ( min( epsi, sqrt( ng ) ) * ng, 0.95 * epsi_max )

   # Initialize the iteration counter, the interior indicator and the step.

   cgits    = 0
   interior = 1
   s        = zeros_like( g )

   # Loop on successive sets of CG iterations (in different faces).

   for iface in range(1,n+1):

      if ( verbose > 0 ):
         disp(' ')
         disp( ' ==== (re)starting CG on face '+str(iface)+' [ nfree = '+str(nfree)+' ]' )
         disp( '      It       value       ||Pg||    face  nfree ' )
         print '      %2d  %.7e  %.2e  %5d  %5d' % (cgits, value, ng, iface, nfree )

      # Reinitialize the first direction to the free negative gradient
      # in the current face.
  
      p = - g * xfree
#      print 'p'
#      print p

      # Nothing more to do in this face if the projected gradient
      # (in p) is small enough.

      if ( iface > 1 ):
         ng = linalg.norm( p )


      if ( ng  < eps_term ):
         if ( interior ):
            msg = 'pr_tcg1 - interior solution'
         else:
            msg = 'pr_tcg1 - boundary solution'

         if ( verbose > 1 ):
            disp(msg)
            
         gplus = g
         return s, msg, cgits, value, gplus, interior

      # CG loop

      for cgits in range(cgitstart,cgitmax+1):

         Hp  = dot(p,H)
         pHp = dot(p, Hp.T)[0][0]
         gp  = dot(g, p.T)[0][0]

         # Piecewise linesearch loop along the current
         # (projected) CG direction.  The loop is on
         # successive breakpoints.
      
         newfix = 0
      
         for j in range(0,n):

            # Find the step to the first bound.

            minsig = 1.0e99
            imin   = -1
            for  i in range(0,n):
               if ( xfree[0,i] ):
                  if ( p[0,i] > 0 ):
                     sigbound = ( ub[0,i] - s[0,i] ) / p[0,i]
                  elif ( p[0,i] < 0 ):
                     sigbound = ( lb[0,i] - s[0,i] ) / p[0,i]
                  else:
                     sigbound = 1.0e100

                  if ( sigbound < minsig ):
                     minsig = sigbound
                     imin   = i

            # Check if negative curvature encountered.

            if ( pHp <= 0 ):
               sigma = minsig
            else:
               alpha = - gp / pHp

               if ( abs(alpha) < minsig ):
                  sigma = alpha
                  imin  = -1
               else:
                  sigma = minsig

            # Compute the trial point, the associated gradient and model value.

            s      = s + sigma * p
            g      = g + sigma * Hp
            value  = value + sigma * gp + 0.5 * sigma**2 * pHp
         
            s_save = append(s_save, s)
            g_save = append(g_save, g)
            v_save = append(v_save, value)

            # The trial point is on the boundary of the current face.

            if ( imin > -1 ):

               # Move to the first bound exactly and update
               # the activity indicators.

               nfb = nfree
               for  i in range(0,n):
                  if ( xfree[0,i] ):
                     if ( logical_and(p[0,i] > 0, abs( s[0,i]-ub[0,i] ) <= eps ) ):
                        s[0,i]     = ub[0,i]
                        xfree[0,i] = 0
                        nfree      = nfree - 1
                     elif ( logical_and(p[0,i] < 0, abs( s[0,i]-lb[0,i] ) <= eps ) ):
                        s[0,i]     = lb[0,i]
                        xfree[0,i] = 0
                        nfree      = nfree - 1

               # Update other quantities. If only one bound
               # has been hit, then the update is cheaper.

               nfix   = nfb - nfree
               newfix = newfix + nfix
               if ( nfix == 1 ):
                  Hp        = Hp  - p[0,imin] * H[:,imin].T
                  p[0,imin] = 0
               else:
                  p   = p * xfree
                  Hp  = dot( p, H )

               gp  = dot( g, p.T)[0][0]
               pHp = dot( p, Hp.T)[0][0]

               
               # Terminate the search loop if the slope is positive or zero
               # (this may happen if no free variable is left).

               if ( gp >= 0 ):
                  if ( nfree == 0 ):
                     if ( verbose > 1 ):
                        disp('      break if no free variables anymore: gp='+str(gp))

                     break
                  else:
                     if ( verbose > 1 ):
                        disp('      !!! slope is positive or zero: gp='+str(gp)+' but there are still free variables')

                     break

            # The new point is inside the current face.

            else:
               if ( verbose > 1 ):
                  disp('      break if new point is inside the current face')

               break;

            if ( verbose > 0 ):
               disp(' ')
               disp( ' ---- new face [ nfree = '+str(nfree)+' ] in piecewise search' )

            # No need to pursue the piecewise search if there are no
            # more free variables.

            if ( nfree == 0 ):
               if ( verbose > 1 ):
                  disp('      break if nfree=0')

               break
      
         # End of the piecewise linesearch: exit CG if new variables
         # have been activated.

         if ( newfix > 0 ):
            interior = 0
            if ( verbose > 1 ):
               disp('      end of piecewise linesearch: exit if new variables have been activated')

            break

         # CG body (1): compute the new free gradient and its norm.

         gfree  = g * xfree
         ngfree = linalg.norm( gfree )

         #  Print out

         if ( verbose > 0 ):
            disp( '      It       value       ||Pg||    face  nfree ' )
            print '      %2d  %.7e  %.2e  %5d  %5d' % (cgits, value, ngfree, iface, nfree )

         # CG termination ?
         # 1) accuracy obtained
            
         if  ( ngfree < eps_term ):
            if ( interior ):
               msg = 'pr_tcg2 - interior solution '
            else:
               msg = 'pr_tcg2 - boundary solution '

            if ( verbose > 1 ):
               disp(msg)
               
            gplus = g
            return s, msg, cgits, value, gplus, interior

         # 2) too many CG iterations

         if ( cgits >= cgitmax ):
            msg   = 'pr_tcg - iteration limit reached'
            if ( verbose > 1 ):
               disp(msg)
               
            gplus = g
            return s, msg, cgits, value, gplus, interior

         # CG body (2): compute the new search direction.
         
         beta   = ( ngfree / ng )**2
         p      = - gfree + beta * p
         ng     = ngfree

      # Update the iteration counter for starting the next CG loop.

      cgitstart = cgits

   if ( nfree == 0 ):
      msg   = 'pr_tcg3 - boundary solution'
   else:
      msg   = 'pr_tcg3 - no solution in max iterations'

   if ( verbose > 1 ):
      disp(msg)
      
   gplus = g
   return s, msg, cgits, value, gplus, interior
