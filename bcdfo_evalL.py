#! /usr/bin/env python
from numpy import *
from bcdfo_evalP import *
import helper

@helper.convertingDecorator
def bcdfo_evalL_( QZ, RZ, Y, choice_set, x, xbase, whichmodel, scale, shift_Y ):
	print "choice_set", choice_set
	print "choice_set[0]", choice_set[0]
	print "type choiceset0", type(choice_set[0])
	choice_set =  [int(gnlph) for gnlph in choice_set[0]]
	return bcdfo_evalL( QZ, RZ, Y, choice_set, x, xbase, whichmodel, scale, shift_Y )

def bcdfo_evalL( QZ, RZ, Y, choice_set, x, xbase, whichmodel, scale, shift_Y ):
###############################################################################
#
#  Computes the values of the Lagrange polynomials at x.
#
#  INPUTS:
#
#  QZ, RZ      : the QR factors of the (possibly shifted) matrix Z containing
#                the polynomial expansion of the interpolation points
#  Y           : current interpolation set Y
#  choice_set  : the indices (of Y's columns) for which to calculate the 
#                Lagrange polynomial values
#  x           : the point at which the model must be evaluated
#  xbase       : the current base point
#  whichmodel  : kind of model/Lagrange polynomial to compute
#  scale       : the current model scaling
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#
#  OUTPUT:
#
#  values      : the values of the Lagrange polynomials at x.
#
#  PROGRAMMING: A. Troeltzsch, September 2010. 
#
#  DEPENDENCIES: bcdfo_evalP
#
#  TEST:
#  Y = [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; 
#  [ QZ, RZ, xbase, scale] = bcdfo_build_QR_of_Y( Y, 0, 1 );
#  values = bcdfo_evalL( QZ, RZ, Y, [2:6], [-1;1], xbase, 0, scale, 1 )
#  should give 
#     values =
#         0   97.0000    2.9900    1.0000 -100.0000   -0.4950 
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
###############################################################################

   n,p1 = shape( Y )
   I      = eye( p1 )
   lc     = len(choice_set)
   q      = ( ( n + 1 ) * ( n + 2 ) ) / 2
   values = zeros((p1,1))
   
   if ( whichmodel == 3 and p1 < q ):
   
      #  for underdetermined regression model use min l2-norm model
      whichmodel = 2;
   
   ###############################################################################
   
   if ( whichmodel == 0 ):
   
      # evaluate (sub-basis) Lagrange polynomial (p1 = q)
      print "choice_set", choice_set  
      #print "I[choice_set,:]", I[choice_set,:]
      #print "linalg.solve( RZ, QZ.transpose())", linalg.solve( RZ, QZ.transpose())						
      print "bcdfo_evalP",bcdfo_evalP( dot(I[choice_set,:], linalg.solve( RZ, QZ.transpose())), x, xbase, scale, shift_Y );
      print " values[choice_set] ",  values[choice_set] 
      print "values", values						
      values[choice_set] = bcdfo_evalP( dot(I[choice_set,:], linalg.solve( RZ, QZ.transpose())), x, xbase, scale, shift_Y );
   
   elif ( whichmodel == 1 ):
   
      if ( p1 == n+1 or p1 == q  ):
   
         # evaluate Minimum L2 norm Lagrange polynomial (p1 == n+1 or p1 == q)
   
         values[choice_set] = bcdfo_evalP( dot(I[choice_set,:], linalg.solve( RZ, QZ.transpose()).transpose()), x, xbase, scale, shift_Y );
   
      else:
   
         # evaluate Minimum Frobenius norm Lagrange polynomial (p1 <= q)
   
         #  If shifting is active, the matrix of the interpoaltion points is first
         #  shifted and scaled (to compute M correctly).
         
         Y2 = copy( Y )
         x2 = copy ( x )
   
         if ( shift_Y ):
            xbase  = Y2[:,0]
            scaleY = 0
            for i in range(0,p1):
               Y2[:,i] = Y2[:,i] - xbase;
               scaleY = max( scaleY, linalg.norm( Y2[:,i] ) )
    
            Y2     = Y2 / scaleY
            
            #xbase = reshape(xbase,
   
         M  = bcdfo_evalZ( Y2, q ).transpose()
         MQ = M[:,n+1:q]
   
         if ( shift_Y ):
            phi = bcdfo_evalZ( ( x2 - xbase ) * scale[0][2], q );
         else:
            phi = bcdfo_evalZ( x2, q );
   
         values[choice_set] =  dot( hstack([I[choice_set,:], zeros((lc,n+1))]), dot( QZ, linalg.solve(RZ.transpose(), vstack([dot( MQ, phi[n+1:q]), phi[0:n+1]]) ) ))
   
   elif ( whichmodel == 2 ):
   
      # evaluate Minimum L2 norm Lagrange polynomial (p1 <= q)
   
      values[choice_set] = bcdfo_evalP( dot(I[choice_set,:], linalg.solve( RZ, QZ.transpose())), x, xbase, scale, shift_Y );
   
   elif ( whichmodel == 3 ):
   
      # evaluate Regression Lagrange polynomial (p1 >= q)
   
      values[choice_set] = bcdfo_evalP( dot( I[choice_set,:], linalg.solve( RZ, QZ.transpose()).transpose()), x, xbase, scale, shift_Y );
   
   return values
   
