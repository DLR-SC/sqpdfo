#!/usr/local/bin/python
from numpy import *
from bcdfo_poisedness_Y import *
from bcdfo_replace_in_Y import *

def bcdfo_repair_Y( QZ, RZ, Y, Delta, farfact, farthr, closethr, 
           eps_L, xbase, lSolver, whichmodel, hardcons, xl, xu, indfree, 
           stratLam, scale, shift_Y, normgx, kappa_ill ):

###############################################################################
#
#  Repairs the interpolation set Y by first replacing the interpolation points
#  at a distance larger that farfact*Delta from Y(:,0), and then replaces 
#  interpolation points with the best possible point obtained by taking the best
#  optimal replacement, until this improvement no longer causes a relative simplex
#  volume increase by at least threshold.  For conventions on how polynomials
#  are represented, see the documentation of evalZ.
#
#  INPUT:
#
#  QZ          : the Q matrix of the QR decomposition of Z(Y)
#  RZ          : the R matrix of the QR decomposition of Z(Y)
#  Y           : a matrix whose columns contain the current interpolation points
#  Delta       : the radius of the ball centered at Y(:,1) in which the
#                replacement point must be found
#  farfact     : multiple of TR radius defining far points
#  farthr      : the minimal acceptable improvement in volume simplex
#                (poisedness) for the new point to be acceptable in
#                replacement of a far away point.
#  closethr    : the minimal acceptable improvement in volume simplex
#                (poisedness) for the new point to be acceptable in
#                replacement of a close point.
#  eps_L       : the relative accuracy on the trust-region constraint for
#                maximization of the Lagrange polynomials
#  xbase       : the current base point
#  lSolver     : linear solver used for the minimization of the model
#  whichmodel  : kind of model/Lagrange polynomial to compute
#  scale       : the current interpolation set scaling
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#  normgx      : infinity norm of the projected gradient
#  kappa_ill   : threshold to declare a system matrix as ill-conditioned
#
#  OUTPUT:
#
#  QZ          : the Q matrix of the QR decomposition of Z(Y) for updated Y
#  RZ          : the R matrix of the QR decomposition of Z(Y) for updated Y
#  Y           : the updated set of interpolation points (column-wise)
#  replaced    : a vector containing the indices of the points of Y which have
#                been replaced
#  maximprove  : the best improvement still possible by replacing a single
#                interpolation point (gives an estimate of poisedness)
#  Y_radius    : the maximum distance from any interpolation point to the 
#                base point
#  xbase       : the updated base point
#  scale       : the updated interpolation set scaling
#
#  PROGRAMMING: Ph. Toint, S. Gratton, April 2009. (This version 22 VI 2009)
#
#  DEPENDENCIES: bcdfo_find_new_yj, bcdfo_gradP, bcdfo_hessP, bcdfo_replace_in_Y
#                bcdfo_poisedness_Y
#
#  TEST:
#  Y = [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ];
#  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y(  Y, 0, 0, 1, 1, 1e15 );
#  [ QZplus, RZplus, Yplus, replaced, maximprove, Y_radius, xbase, scale ] = ...
#            bcdfo_repair_Y( QZ, RZ, Y, 0.7, 10, 1.0e-10, 1.1, 0.001, xbase,...
#            1, 0, 0, [-10,-10], [10,10], [1,2], 1, scale, 0, 1, 1e15 )
#  QZplus =
#
#   -1.0000         0         0         0         0         0
#         0    0.6596    0.5831    0.2590   -0.3837   -0.1032
#         0    0.6401   -0.7030    0.1160    0.1723   -0.2299
#         0   -0.1657    0.0581    0.9007    0.3837    0.1032
#         0   -0.1560    0.2178   -0.0539    0.2743   -0.9220
#         0   -0.3216   -0.3391    0.3244   -0.7750   -0.2752
#
#  RZplus =
#
#   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000
#         0   -0.7616   -0.1207    0.9878   -0.1412    0.9682
#         0         0    0.7473    1.2823   -0.5045   -0.9705
#         0         0         0    2.3194   -0.0241    0.1242
#         0         0         0         0    0.5459    0.8931
#         0         0         0         0         0   -2.3039
#
#  Yplus =
#
#         0   -0.5023    0.3561    2.0000   -0.6030         0
#         0   -0.4875   -0.6026         0    0.3555    2.0000
#
#  replaced =
#
#     2     3     5
#
#  maximprove =
#
#    1.0002
# 
#  Y_radius =
#
#     2
#
#  xbase =
#
#     0
#     0
#
#  scale =
#
#     1
#     1
#     1
#     1
#     1
#     1
#
#  The scaled version: 
#  Y = [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ];
#  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y(  Y, 0, 1, 1, 1, 1e15 );
#  [ QZplus, RZplus, Yplus, replaced, maximprove, Y_radius, xbase, scale ] = ...
#           bcdfo_repair_Y( QZ, RZ, Y, 0.7, 10, 1.0e-10, 1.1, 0.001, xbase,...
#           1, 0, 0, [-10,-10], [10,10], [1,2], 1, scale, 1, 1, 1e15 )
#  QZplus =
#
#   -1.0000         0         0         0         0         0
#         0    0.7017    0.6633    0.1499   -0.2060   -0.0529
#         0    0.6810   -0.7144    0.0589    0.0925   -0.1178
#         0   -0.0881    0.0262    0.9003    0.4120    0.1058
#         0   -0.0830    0.1119   -0.0352    0.2945   -0.9448
#         0   -0.1710   -0.1910    0.4028   -0.8322   -0.2820
#
#
#  RZplus =
#
#   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000
#         0   -0.3580   -0.0762    0.6576   -0.0867    0.6395
#         0         0    0.3491    0.6764   -0.3138   -0.6584
#         0         0         0    0.6001   -0.0159    0.0413
#         0         0         0         0    0.1465    0.2397
#         0         0         0         0         0   -0.5902
#
#  Yplus =
#
#         0   -0.5023    0.3561    2.0000   -0.6030         0
#         0   -0.4875   -0.6026         0    0.3555    2.0000
#
#  replaced =
#
#     2     3     5
#
#  maximprove =
#
#    1.0002
#
#  Y_radius =
#
#     2
#
#
#  xbase =
#
#     0
#     0
#
#  scale =
#
#    1.0000
#    0.5000
#    0.5000
#    0.2500
#    0.2500
#    0.2500
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
################################################################################

   verbose           = 0         # 1 for debug
   Ycopy = copy(Y)
   max_improve_loops = 20        # the max number of times every close point is 
                                  # considered for improvement
   
   if ( verbose ):
     print '--- enter bcdfo_repair_Y for Delta = '+str(Delta)+' ---'
   
   n, p1 = shape( Ycopy )
   replaced  = array([])
   
   #  Compute the distance of the current interpolation points to the base point.
   
   d = zeros((p1))
   for j in range(1,p1):
      d[j] = linalg.norm( Ycopy[:,j]- Ycopy[:,0] )

   dsorted = sort( d )[ : :-1]
   jsorted = argsort( d )[ : :-1]
   
   if ( verbose ):
      print 'd='
      print d
      print 'dsorted='
      print dsorted
      print 'jsorted='
      print jsorted

   #  First replace the distant points.

   for j in range(0,p1):
      if ( dsorted[j] > farfact * (1 + eps_L) * Delta ):
         jmax = jsorted[j]
              
         if ( hardcons == 1 ):
            y, improvement = bcdfo_find_new_yj_bc( QZ, RZ, Ycopy, jmax, Delta, 
               eps_L, xbase, lSolver, whichmodel, xl, xu, indfree, 
               stratLam, scale, shift_Y )
         else:
            y, improvement = bcdfo_find_new_yj( QZ, RZ, Ycopy, jmax, Delta, eps_L,
               xbase, lSolver, whichmodel, scale, shift_Y)
      
         if ( verbose ):
            print 'lambda = '+str(improvement)+' at j='+str(jmax)
      
         if ( improvement > farthr ):        
            QZ, RZ, Ycopy, xbase, scale = bcdfo_replace_in_Y( QZ, RZ, y, Ycopy, jmax, xbase,
               whichmodel, scale, shift_Y, Delta, normgx, kappa_ill )
            replaced  = append( replaced, jmax )
            d[jmax] = linalg.norm( y - Ycopy[:,0] ) 

      else:
         Y_radius = dsorted[j]
         break

   if ( verbose ):
      print 'replaced='
      print replaced
      poised, Y_radius = bcdfo_poisedness_Y( QZ, RZ, Ycopy, eps_L, xbase, lSolver,
             whichmodel, hardcons, xl, xu, indfree, stratLam, scale, shift_Y )
      print 'after distant replace: poisedness(Y) = '+str(poised)+' Y_radius  = '+str(Y_radius)

   #  Perform a loop over possible optimal improvements.

   for k in range(0,max_improve_loops):

      maximprove = 0

      #  Loop on all the possible replacements to find the best one.

      for j in range(1,p1):

         if ( hardcons == 1 ):
            y, improvement = bcdfo_find_new_yj_bc( QZ, RZ, Ycopy, j, Delta, 
               eps_L, xbase, lSolver, whichmodel, xl, xu, indfree, 
               stratLam, scale, shift_Y );
         else:
            y, improvement = bcdfo_find_new_yj( QZ, RZ, Ycopy, j, Delta, eps_L, 
               xbase, lSolver, whichmodel, scale, shift_Y )

         #  Remember the current polynomial value, index and replacement point is
         #  this is the best so far

         if ( verbose > 1 ):
            print ' ==> j = '+str(j)+' improve = '+str(improvement)
            print 'y='
            print y

         if ( improvement > maximprove ):
            maximprove = improvement
            jmax       = j
            ymax       = y

      #  If no significant improvement was found, return after updating the 
      #  interpolation radius.

      if ( maximprove  < closethr or jmax == 0 ):
         Y_radius = max( d )           # recompute Y_radius after the exchanges
         if ( verbose ):
            print 'replaced='
            print replaced
            print 'maximprove(small)= '+str(maximprove)+', jmax= '+str(jmax)
            poised, Y_radius = bcdfo_poisedness_Y( QZ, RZ, Ycopy, eps_L, xbase, lSolver,
               whichmodel, hardcons, xl, xu, indfree, stratLam, scale, shift_Y )
            print 'after everything: poisedness(Y) = '+str(poised)+' Y_radius  = '+str(Y_radius)
            print '--- exit 1 bcdfo_repair_Y in round k='+str(k)+' ---'

         if ( not replaced.any() ):
            maximprove = 0
             
         return QZ, RZ, Ycopy, replaced, maximprove, Y_radius, xbase, scale

      #  Perform the best replacement.

      if ( verbose ):
         print 'maximprove= '+str(maximprove)+', jmax= '+str(jmax)

      QZ, RZ, Ycopy, xbase, scale = bcdfo_replace_in_Y( QZ, RZ, ymax, Ycopy, jmax, xbase, 
         whichmodel, scale, shift_Y, Delta, normgx, kappa_ill )

      d[jmax] = linalg.norm( ymax - Ycopy[:,0] )

      if ( len( where( replaced == jmax )[0] ) == 0 ):
         replaced = append( replaced, jmax )
         
   # recompute Y_radius after the exchanges
   
   Y_radius = max( d )

   if ( verbose ):
      print 'replaced='
      print replaced
      poised, Y_radius = bcdfo_poisedness_Y( QZ, RZ, Ycopy, eps_L, xbase, lSolver, 
            whichmodel, hardcons, xl, xu, indfree, stratLam, scale, shift_Y );
      print 'after everything: poisedness(Y) = '+str(poised)+' Y_radius  = '+str(Y_radius)
      print '--- exit 2 bcdfo_repair_Y after round k='+str(k)+'  ---'

   return QZ, RZ, Ycopy, replaced, maximprove, Y_radius, xbase, scale

