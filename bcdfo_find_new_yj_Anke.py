#!/usr/local/bin/python
from numpy import *
#import fpconst
from bcdfo_computeLj import *
from bcdfo_gradP import *
from bcdfo_hessP import *
from bcdfo_solve_TR_MS import *
import helper

@helper.convertingDecorator
def bcdfo_find_new_yj_( QZ, RZ, Y, j, Delta, eps_L, xbase, lSolver, whichmodel,
   scale, shift_Y ):
# j is an index, so we have to decrease it.				
	return bcdfo_find_new_yj( QZ, RZ, Y, j-1, Delta, eps_L, xbase, lSolver, whichmodel,
   scale, shift_Y )

def bcdfo_find_new_yj( QZ, RZ, Y, j, Delta, eps_L, xbase, lSolver, whichmodel,
   scale, shift_Y ):
###############################################################################
#
#  Computes a point which is best to replace yj, the j-th (j>1) column of Y (the
#  base for the j-th polynomial) in a ball of radius Delta centered at the
#  first column of Y.  This is achieved by maximizing the absolute value of
#  the j-th Lagrange polynomial in that ball.
#
#  For conventions on how polynomals are represented, see the documentation of
#  evalZ.
#
#  INPUT:
#
#  QZ          : the Q matrix of the QR decomposition of Z(Y)
#  RZ          : the R matrix of the QR decomposition of Z(Y)
#  Y           : a matrix whose columns contain the current interpolation points
#  j           : the index of the interpolation point one wishes to replace (j > 1)
#  Delta       : the radius of the ball centered at Y(:,1) in which the
#                replacement point must be found
#  eps_L       : the relative accuracy on the trust-region constraint for
#                maximization of the Lagrange polynomials
#  xbase       : the current base point
#  lSolver     : linear solver used for the minimization of the model
#  whichmodel  : kind of model/Lagrange polynomial to compute
#  scale       : the current interpolation set scaling
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#
#  OUTPUT:
#
#  ynew        : the best replacement for Y(:,j)
#  improvement : the improvement in poisedness obtained by the update, which
#                is equal to |L_j(new y)|. If this value is smaller than the
#                threshold input parameter, L and X are unchanged by the
#                procedure.
#
#  PROGRAMMING: Ph. Toint, February 2009. (This version 22 VI 2009)
#
#  USES: bcdfo_gradP, bcdfo_hessP, bcdfo_solve_TR_MS
#
#  TEST:
#  Y = [ 3 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; whichmodel = 0;
#  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y , whichmodel, 0 );
#  [ ynew, improvement ] = bcdfo_find_new_yj( QZ, RZ, Y, 5, 1.0, 0.001, xbase, 1, ...
#       whichmodel, scale, 0 )
#  should give
#  ynew =
#
#    3.2808
#   -0.9598
#
#  improvement =
#
#  314.8825
#
#  the same must be obtained by the shifted and scaled version:
#  Y = [ 3 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; whichmodel = 0;
#  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y , whichmodel, 1 );
#  [ ynew, improvement ] = bcdfo_find_new_yj( QZ, RZ, Y, 5, 1.0, 0.001, xbase, 1, ...
#       whichmodel, scale, 1 )
#
################################################################################

   verbose     = 0             # 1 for debug
   n,P1        = shape(Y)
   improvement = 0
   Ycp         = copy(Y)

   if ( verbose ):
      disp('--------- enter find_new_yj ')

   if ( j < 1 ):           # never attempt to replace the current iterate.
      msgTR = 'The index to replace is the current iterate! Return.'
      if ( verbose ):
         disp( msgTR )
      return zeros((n,1),float), improvement

   #  Get the j-th Lagrange polynomial

   Lj = bcdfo_computeLj( QZ, RZ, j, Ycp, whichmodel, scale, shift_Y )
   
   if ((not (isreal(Lj)).any()) or isnan(Lj).any() or isinf(Lj).any()):
      msgTR = 'Error0: Lagrange polynomial contains NaN or Inf or nonreal components!!';
      if ( verbose ):
         disp( msgTR )

      return zeros((n,1),float), improvement

   #  Maximize Lj in a larger 2-norm TR if using infty-norm in the local solver (CG)

   if ( lSolver == 2 ):
      Delta = sqrt( n ) * Delta

   #  Get the polynomial's gradient and Hessian at the current iterate.

   if ( shift_Y ):

      #  When shifted, the origin in the scaled variables corresponds
      #  to Y(:,0) in the unscaled space

      g  = bcdfo_gradP( Lj, zeros((n,1),float), xbase, scale, 0 )
      H  = bcdfo_hessP( Lj, zeros((n,1),float), xbase, scale, 0 )
   
      #  Minimize this polynomial and its opposite.

      [ pstep, lamb, norms, pvalue, gplus, nfact, neigd, msgTR, hardcase ] = bcdfo_solve_TR_MS(
         g, H, Delta*scale[2][0], eps_L )
      pstep = pstep / scale[2][0]
      [ mstep, lamb, norms, mvalue, gplus, nfact, neigd, msgTR, hardcase ] = bcdfo_solve_TR_MS(
         -g, -H, Delta*scale[2][0], eps_L )
      mstep = mstep / scale[2][0]
      
   else:

      #  When no shift occurs, the current iterate is Y(:,0)
      #print "Ycp\n", Ycp
      #print "Ycp[:,0]\n", Ycp[:,0]

      g  = bcdfo_gradP( Lj, Ycp[:,0].reshape((n,1)), xbase, scale, 0 );
      H  = bcdfo_hessP( Lj, Ycp[:,0].reshape((n,1)), xbase, scale, 0 );

      #  Minimize this polynomial and its opposite.

      [ pstep, lamb, norms, pvalue, gplus, nfact, neigd, msgTR, hardcase ] = bcdfo_solve_TR_MS(
         g, H, Delta, eps_L )
      [ mstep, lamb, norms, mvalue, gplus, nfact, neigd, msgTR, hardcase ] = bcdfo_solve_TR_MS(
         -g, -H, Delta, eps_L )


   if ( verbose ):
      print ' === find_new_yj: j = '+ str(j)+ ' positive value = '+ str(pvalue)+' step:'
      print pstep
      print ' === find_new_yj: j = '+ str(j)+ ' negative value = '+ str(mvalue)+' step:' 
      print mstep

   #  Select the maximum in absolute value.  
   
   if ( mvalue < pvalue ):
      improvement = abs( mvalue );
      ynew        = Ycp[:,0] + mstep;
   else:
      improvement = abs( pvalue );
      ynew        = Ycp[:,0] + pstep;
      
   ynew = ynew.T

   if ( verbose ):
   
      print 'ynew'
      print ynew
      disp('--------- exit find_new_yj ')

   return ynew, improvement
   
