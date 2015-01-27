#!/usr/local/bin/python
from numpy import *
import helper

@helper.convertingDecorator
def bcdfo_solve_TR_MS_( gg, HH, Delta, eps_D ):
	return bcdfo_solve_TR_MS( gg, HH, Delta, eps_D )

def bcdfo_solve_TR_MS( gg, HH, Delta, eps_D ):

###############################################################################

#  A simple implementation of exact trust-region minimization based on the
#  More-Sorensen algorithm.

#  INPUT: 

#  g        : the model's gradient
#  H        : the model's Hessian
#  Delta    : the trust-region's radius
#  eps_D    : the accuracy required on the equation ||s|| = Delta for a
#             boundary solution

#  OUTPUT:

#  s        : the trust-region step
#  lamb   : the Lagrange multiplier corresponding to the trust-region constraint
#  norms    : the norm of s
#  value    : the value of the model at the optimal solution
#  gplus    : the value of the model's gradient at the optimal solution
#  nfact    : the number of Cholesky factorization required
#  neigd    : the number of eigenvalue decompositions required
#  msg      : an information message

#  DEPENDENCIES: -

#  PROGRAMMING: Ph. Toint and S. Gratton, April 2009. (This version 14 I 2010)

#  TEST:

#  bcdfo_solve_TR_MS( [ 2 ; 3 ], [ 4 6; 6 5 ], 1.0, 0.001 )
#  should give
#    0.5153
#   -0.8575
#  bcdfo_solve_TR_MS( [ 2 ; 0 ], [ 4 0; 0 -15 ], 1.0, 0.001 )
#  should give
#   -0.1053
#    0.9944

#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.

###############################################################################

   verbose  = 0
   theta    = 1.0e-13;         # accuracy of the interval on lambda
   eps_g    = 1.0e-15;        # theshold under which the gradient is considered as zero
   nitmax   = 300;            # the maximum number of MS iterations
   n,dum    = shape( HH );    # space dimension
   g        = copy(gg)
   H        = copy(HH)
   s        = zeros_like(gg)     # initialize zero step
   norms    = 0;              # initial stepnorm
   lamb     = 0;              # initial lambda
   value    = 0;              # initial model value
   gplus    = copy(gg)              # initialize gplus
   nfact    = 0;              # factorization counter
   neigd    = 0;              # eigen decomposition counter
   hardcase = 0;              # hard case

   if ( verbose ):
      disp( ' bcdfo_solve_TR_MS : ============ enter' )

   #  Compute initial bounds on lambda.

   gnorm    = linalg.norm( g );
   goverD   = gnorm / Delta;
   Hnorminf = linalg.norm( H, inf );
   if ( Hnorminf > 0 ):        # because Octave generates a NaN for null matrices.
      HnormF   = linalg.norm( H, 'fro' );
   else:
      HnormF   = 0;

   lower    = max( 0, goverD - min( Hnorminf, HnormF ) );
   upper    = max( 0, goverD + min( Hnorminf, HnormF ) );

   #  Compute the interval of acceptable step norms.

   Dlower   = ( 1 - eps_D ) * Delta;
   Dupper   = ( 1 + eps_D ) * Delta;

   #  Zero gradient

   if ( gnorm < eps_g ):
      if ( verbose ):
         disp( ' bcdfo_solve_TR_MS : zero gradient:' )
         
      D, V    = linalg.eig( H )
      neigd   = neigd + 1;
      mu      = min( D )
      imu     = argmin( D )
  
      if ( mu < 0 ):
         s = Delta * V[:,imu].T
      else:
         if ( gnorm == 0 ):
            s = zeros_like(g)
         else:
            s = - Delta * ( g / gnorm );
   
      sfound = 1;
      norms  = linalg.norm( s );
      lamb = -mu;

   #  Nonzero gradient

   else:

      if ( verbose ):
         disp( ' bcdfo_solve_TR_MS : nonzero gradient:' )

      #  Compute initial lambda.

      if ( lower == 0 ):
         lamb = 0;
      else:
         lamb = max( sqrt( lower * upper ), lower + theta * ( upper - lower ) );

      #  Loop on successive trial values for lambda.

      for i in range(0,nitmax):
         new_lamb = -1;
         sfound     = 0;
         if ( verbose ):
            print ' bcdfo_solve_TR_MS ('+str(i)+'): lower = '+  str( lower )+' lambda = '+ str( lamb )+' upper = '+  str( upper )

         #  Factorize H + lambda I.

         try:
             R = linalg.cholesky( H + lamb * eye( n ) ).T
             p = 0
         except:
             R = array([[]])
             p = 1

         if ( isnan(R.any()) ):
            disp( 'NaN in Cholesky factorization' )
            msg = 'error4';
            if ( verbose ):
               disp( ' bcdfo_solve_TR_MS : ============ error exit' )
        
            return s, lamb, norms, value, gplus, nfact, neigd, msg, hardcase
     
         nfact    = nfact + 1;

         #  Successful factorization

         if ( p == 0 ):
            gT     = g.T
            s      = - linalg.solve(R, linalg.solve( R.transpose(), gT )).T
            sfound = 1;
            norms  = linalg.norm( s );
            if ( verbose ):
               print ' bcdfo_solve_TR_MS ('+str(i)+'): ||s|| = '+ str( norms )+ ' Delta  = '+ str( Delta )
         

            #  Test for successful termination.

            if ( logical_or( logical_and( lamb <= eps_g, norms <= Dupper ), logical_and( norms >= Dlower, norms <= Dupper ) ) ):

               #  Compute the optimal model value and its gradient.

               w     = dot(H, s.T)
               sT    = s.T
               #print "g:\n", g, "\nsT:\n", sT, "\ns:\n", s, "\nw:\n", w															
               value = dot(g, sT)  + 0.5 * dot(s, w)
               gplus = g + w.T;
               norms = linalg.norm( s );
   
               #  Define information message.

               if ( norms < ( 1 - eps_D ) * Delta ):
                  msg = 'interior solution';
               else:
                  msg = 'boundary solution';
           
               if ( verbose ):
                  disp( ' bcdfo_solve_TR_MS : ============ successful exit' )
            
               return s, lamb, norms, value, gplus, nfact, neigd, msg, hardcase
         

            #  Newton's iteration on the secular equation

         
            w = linalg.solve(R.transpose(), s.T)
            normw2     = dot(w.transpose(),w)[0][0]
            new_lamb = lamb + ( ( norms - Delta ) / Delta ) * ( norms**2 / normw2 )

            #  Check feasibility wrt lambda interval.

            if ( norms > Dupper ):
               lower = lamb
            else:
               upper = lamb
         
            theta_range = theta * ( upper - lower );
            if ( new_lamb > lower + theta_range and new_lamb < upper - theta_range ):
               lamb = new_lamb
            else:
               lamb = max( sqrt( lower * upper ), lower + theta_range );

         #  Unsuccessful factorization: take new lambda as the middle 
         #  of the allowed interval

         else:
            lower  = lamb;
            t      = 0.5;
            lamb = ( 1 - t ) * lower + t * upper;

         #  Terminate the loop if the lambda interval has shrunk to meaningless.

         if ( upper - lower < theta * max( 1, upper ) ):
            break

   #  The pseudo-hard case

   #  Find eigen decomposition and the minimal eigenvalue.

   D, V    = linalg.eig( H )
   neigd   = neigd + 1;
   mu      = min( D )
   imu     = argmin( D )

   #  create a matrix with the eigenvalues on the diagonal
   D = diag(D)

   if ( verbose ):
      gamma   = abs(dot(V[:,imu],g.T)[0])
      print ' bcdfo_solve_TR_MS : pseudo hard case: gamma = '+ str(gamma)+ ' ||g|| = '+ str(linalg.norm(g))

   #  Compute the critical step and its orthogonal complement along the
   #  eigenvectors corresponding to the most negative eigenvalue
  
   D        = D - mu * eye(n)
   maxdiag  = max( diag( D ) );
   ii       = where( abs( diag( D ) ) < 1e-10 * maxdiag )
   if ( len( ii[0] ) < n and len( ii[0] ) > 0 ):
      D[ii,ii]    = 0.5 * maxdiag * eye(len(ii[0]))
      Dinv        = linalg.inv(D)
      Dinv[ii,ii] = 0
      scri        = - dot(dot(dot(V,Dinv),V.transpose()),g.transpose())
      nscri       = linalg.norm( scri );
   else:
      scri  = zeros(( n, 1 ))
      nscri = 0

   if ( nscri <= Delta ):
      p2 = poly1d([linalg.norm(V[:,imu] )**2, 2*dot(V[:,imu].transpose(),scri)[0], nscri**2 - Delta**2 ],r=0)
      root = max(p2.r)
      s     = scri.T + root*V[:,imu]
   else:
      s     = Delta * scri.T / nscri
 
   lamb   = -mu;
   if ( verbose ):
      print ' bcdfo_solve_TR_MS : ||scri|| = '+ str( linalg.norm(scri) )+ ' lambda = '+ str( lamb)

   hardcase = 1;

   #  Compute the model value and its gradient.
   sT    = s.T
   w     = dot( H, sT)
   value = dot( g, sT ) + 0.5 * dot( s, w)
   gplus = g + w.T
   norms = linalg.norm( s )

   #  Define information message.

   if ( norms < ( 1 - eps_D ) * Delta ):
      msg = [ 'interior solution ( '+ str( nfact )+' factorizations,  lambda = '+ str( lamb )+ ')' ]
   else:
      msg = [ 'boundary solution ( '+ str( nfact )+ ' factorizations, '+str( neigd )+ ' eigen decomposition, lambda = '+ str( lamb )+ ' )' ]

   if ( verbose ):
      disp( ' bcdfo_solve_TR_MS : ============ hard case exit' )

   return s, lamb, norms, value, gplus, nfact, neigd, msg, hardcase
   
