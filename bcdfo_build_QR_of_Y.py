#! /usr/bin/env python
from numpy import *
from bcdfo_evalZ import *
from bcdfo_checkZ import *
import helper
#from runtime import matlabarray

def bcdfo_build_QR_of_Y_( Y, whichmodel, shift_Y, Delta, normgx, kappa_ill, nargout=None ):
		Y = helper.convert(Y)
		whichmodel = helper.convert(whichmodel)
		shift_Y = helper.convert(shift_Y)
		Delta = helper.convert(Delta)
		normgx = helper.convert(normgx)
		kappa_ill = helper.convert(kappa_ill)
		
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y(Y, whichmodel, shift_Y, Delta, normgx, kappa_ill)
				
		QZ = matlabarray(QZ)
		RZ = matlabarray(RZ)
		xbase = matlabarray(xbase)
		scale = matlabarray(scale)
		
		return QZ, RZ, xbase, scale
	

def bcdfo_build_QR_of_Y( Y, whichmodel, shift_Y, Delta, normgx, kappa_ill ):
###############################################################################
#
#  Computes the QR factorization of the (possibly shifted) matrix containing
#  the polynomial expansion of the interpolation points. If shifting is
#  required, (i.e. if shift_Y is true) the matrix being factorized has columns
#  containing (Y(:,j)-Y(:,1))/ scale, where scale is the max( norm(Y(:,j)-Y(:,1)).
#
#  INPUT:
#
#  Y           : a matrix whose columns contain the current interpolation points
#  whichmodel  : kind of model to build
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#  Delta       : trust-region radius
#  normgx      : infinity norm of the projected gradient
#  kappa_ill   : threshold to declare a system matrix as ill-conditioned
#
#  OUTPUT:
#
#  QZ, RZ      : the QR factors of the (possibly shifted) matrix containing
#                the polynomial expansion of the interpolation points,
#  xbase       : the base point,
#  scale       : the model diagonal scaling.
#
#  PROGRAMMING: A. Troeltzsch, Ph. Toint, S. Gratton, 2009-2011.
#               (This version 14 I 2011)
#
#  DEPENDENCIES: bcdfo_evalZ
#
#  TEST:
#  Y = [ 1 2 1 3 3 1 ; 1 2 2 1 1.01 3 ];
#  [QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y, 0, 1, 1, 1, 1e15 );
#  model = ( QZ * ( RZ' \ [1 2 3 4 5 6 ]' ) )';
#  model * bcdfo_evalZ( ([1;3]-xbase)*scale(2),6)
#  should give 6.0
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
################################################################################

   n,p1 = shape( Y )

   #  Compute and check the size of the tolerance.

   if ( normgx == 0.0 ):
      delt = Delta**2 / 1e-12;
   else:
      delt = Delta**2 / normgx;

   if ( delt > 1.0e-1 ):
      delt = 1.0e-1;

   if ( delt < 1.0e-10 ):
      delt = 1.0e-10;

   #  Define the size of the polynomial

   if ( whichmodel == 0 ):
      q = p1;
   else:
      q = ( ( n + 1 ) * ( n + 2 ) ) / 2;

   #  for underdetermined regression model use min l2-norm model

   if ( whichmodel == 3 and p1 < q ):
      whichmodel = 2;

   #  If shifting is active, the matrix of the interpoaltion points is first
   #  shifted and scaled, and the base point and scaling factors defined.

   Y2 = copy(Y)

   if ( shift_Y and ( p1 > 1 ) ):
      xbase  = copy(Y2[:,0])
      scaleY = 0;
      scale = array([[1]])
      lst = range(0,p1)
      for i in lst:
         Y2[:,i] = Y2[:,i] - xbase;
         if ( linalg.norm( Y2[:,i] ) > scaleY ):
            scaleY = linalg.norm( Y2[:,i] )

      scale = append(scale, scaleY**(-1) * ones((1,min(n,q-1)),float))
      scale = append(scale, scaleY**(-2) * ones((1,q-n-1), float))
      scale = scale.reshape(q,1)
      
      Y2     = Y2 / scaleY; 
      xbase = xbase.reshape(n,1)
      
      #  Otherwise, the base point and scaling factors are trivial.

   else:
      scale = ones(( q, 1 ),float);
      xbase = zeros(( n, 1 ), float);

   ###############################################################################

   #  Perform factorization of the Z matrix corresponding to the (possibly
   #  shifted and scaled) interpolation matrix.

   if ( whichmodel == 0 ):

      # QR of (sub-basis) interpolation matrix (p1 = q) or
      # (factorize matrix Z = M')

      Z = bcdfo_evalZ( Y2, q );

      #  Check condition of Z and cure if ill-conditioned

      badcond = bcdfo_checkZ( Z, kappa_ill )

      if ( badcond ):
         [U,Sdiag,Vh]        = linalg.svd(Z)
         V = Vh.T
         Sdiag[Sdiag < delt] = delt
         S              = diag(Sdiag)
         M              = (dot(dot(V, S), U.T)).T
         [ QZ, RZ ]     = qr( M );
      else:
         [ QZ, RZ ]     = qr( Z );
   
   elif ( whichmodel == 1 ):

      #  Mixed model: minimum Frobenius-norm model (when underdetermined) and
      #  minimum l2-norm model (at linear and quadratic degree)

      if ( p1 == n+1 or p1 == q ):

         F = bcdfo_evalZ( Y2, p1 ).transpose();

      else:

         # QR of Minimum Frobenius norm interpolation matrix (p1 <= q)
         # (factorize matrix F = [MQMQ' ML; ML' 0])

         M  = bcdfo_evalZ( Y2, q ).transpose();
         ML = M[ :, 0:n+1 ]
         MQ = M[ :, n:q ]
         F  = MQ * MQ.transpose()
         F  = append( F, ML, axis=1 )
         F1 = ML.transpose()
         F1 = append( F1, zeros((n+1,n+1),float), axis=1 )
         F  = append( F, F1, axis=0 )

      #  Check condition of Z and cure if ill-conditioned

      badcond = bcdfo_checkZ( F, kappa_ill )
      
      if ( badcond ):
         [U,S,V]        = svd(F);
         Sdiag          = diag(S);
         indices        = where(Sdiag < delt);
         Sdiag[indices] = delt;
         S              = diag(Sdiag);
         M              = (V * S * U.transpose()).transpose();
         [ QZ, RZ ]     = qr( M );
      else:
         [ QZ, RZ ]     = qr( F );

   elif ( whichmodel == 2 ):

      # QR of Minimum L2 norm interpolation matrix (p1 <= q)
      # (factorize matrix Z = M')

      Z = bcdfo_evalZ( Y2, q );

      #  Check condition of Z and cure if ill-conditioned

      badcond = bcdfo_checkZ( Z, kappa_ill )
      
      if ( badcond ):
         [U,S,V]        = svd(Z);
         Sdiag          = diag(S);
         indices        = where(Sdiag < delt);
         Sdiag[indices] = delt;
         S              = diag(Sdiag);
         M              = (V * S * U.transpose()).transpose();
         [ QZ, RZ ]     = qr( M );
      else:
         [ QZ, RZ ]     = qr( Z );   

   elif ( whichmodel == 3 ):

      # QR of Regression interpolation matrix (p1 >= q)
      # (factorize matrix Z = M)

      Z =  bcdfo_evalZ( Y2, q ).transpose() ;

      #  Check condition of Z and cure if ill-conditioned

      badcond = bcdfo_checkZ( Z, kappa_ill )
      
      if ( badcond ):
         [U,S,V]        = svd(Z);
         Sdiag          = diag(S);
         indices        = where(Sdiag < delt);
         Sdiag[indices] = delt;
         S              = diag(Sdiag);
         M              = (V * S * U.transpose()).transpose();
         [ QZ, RZ ]     = qr( M );
      else:
         [ QZ, RZ ]     = qr( Z );
      

   elif ( whichmodel == 8 ):

      # QR of gmodel interpolation matrix
      # (factorize matrix Z = M')

      Z = bcdfo_evalZ( Y2, q );

      #  Check condition of Z and cure if ill-conditioned

      badcond = bcdfo_checkZ( Z, kappa_ill )

      if ( badcond ):
         [U,S,V]        = svd(Z);
         Sdiag          = diag(S);
         indices        = where(Sdiag < delt);
         Sdiag[indices] = delt;
         S              = diag(Sdiag);
         M              = (V * S * U.transpose()).transpose();
         [ QZ, RZ ]     = qr( M );
      else:
         [ QZ, RZ ]     = qr( Z );
        

   elif ( whichmodel == 9 ):

      # QR of gmodel interpolation matrix
      # (factorize matrix Z = M')

      Z = bcdfo_evalZ( Y2, q );

      #  Check condition of Z and cure if ill-conditioned

      badcond = bcdfo_checkZ( Z, kappa_ill )
      
      if ( badcond ):
         [U,S,V]        = svd(Z);
         Sdiag          = diag(S);
         indices        = where(Sdiag < delt);
         Sdiag[indices] = delt;
         S              = diag(Sdiag);
         M              = (V * S * U.transpose()).transpose();
         [ QZ, RZ ]     = qr( M );
      else:
         [ QZ, RZ ]     = qr( Z );

   return QZ, RZ, xbase, scale
