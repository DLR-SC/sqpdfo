#! /usr/bin/env python
from numpy import *

def bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, P_old, ind_Y,
         i_xold, i_xplus, g, scale, shift_Y, Delta0, indfree, gmat, hessian, 
         epsilon, noisy ):
###############################################################################
#
#  Computes the polynomial P, where P is then represented by a row vector
#  containing its coefficients for the successive monomials.
#  More specifically, these values are:
#  P(1)        : constant coefficient,
#  P(2:n+1)    : coefficients for the linear terms in x(1)... x(n),
#  P(n+2:2n+2) : coefficients for the squared terms in x(1)^2 ... x(n)^2
#  P(2n+3,3n+2): coefficients for the quadratic terms of the first subdiagonal:
#                in x(1)*x(2) ... x(n-1)*x(n)
#  P(3n+3,4n+1): coefficients for the quadratic terms of the second subdiagonal:
#                in x(1)*x(3) ... x(n-2)*x(n)
#  etc.
#
#  INPUTS:
#
#  QZ, RZ      : the QR factors of the (possibly shifted) matrix Z containing
#                the polynomial expansion of the interpolation points
#  Y           : current interpolation set
#  fY          : function values of the interpolation points
#  whichmodel  : the kind of model to build
#  P_old       : polynomial coefficients of the former polynomial
#  scale       : scaling of the interpolation matrix
#  shift_Y     : shift of the matrix of interpolation points
#
#  OUTPUT:
#
#  P           : a row vector containing the coefficients of the polynomial
#
#  DEPENDENCIES: -
#
#  PROGRAMMING: A. Troeltzsch, September 2010.
#
#  TEST:
#  Y = [ 0 1 0 2 0 ; 0 0 1 0 2 ]; fY = [ 1 2 3 4 5 ];
#  whichmodel = 0;
#  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y(  Y, whichmodel, 1 );
#  P = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, [0, 0, 0, 0, 0], 0, ...
#      1, 1, [0, 0], scale, 1, 1, [1 2], [], [1 0,0 1], 1e-5 )
#  should give:
#  P =
#   1.0000   1.0000   4.0000   4.0000       0
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
###############################################################################
#   print QZ 
#   print RZ

   Y2 = copy(Y)
   [n,p1] = shape(Y)
   badcond = 0
   q = ( ( n + 1 ) * ( n + 2 ) ) / 2
   
   #  as long as regression model is underdetermined: use min l2-norm model

   if ( whichmodel == 3 and p1 < q ):
      whichmodel = 2

      # don't use gradient if Delta too small in a noisy function

   if ( whichmodel == 9 and noisy == 1 and Delta0 <= epsilon ):
      whichmodel = 2

###############################################################################

   if ( whichmodel == 0 ):

      # build (sub-basis) model (p1 = q)
      # (QZ and RZ are the factors of Z = M')

      P = reshape(dot(QZ, linalg.solve( RZ.transpose(),fY)),(1,p1))

#      P = ( QZ * ( RZ.transpose() \ fY.transpose() ) ).transpose()

   elif ( whichmodel == 1 ):

      #  build mixed model: Minimum Frobenius norm model (p1 <= q)
      #  and l2-norm model (p1 == n+1 or p1 == q)

      if ( p1 == n+1 or p1 == q ):

         P = linalg.solve(RZ,dot(QZ.transpose(), fY.transpose()))
#         P = ( RZ \ ( QZ.transpose() * fY.transpose() ) ).transpose()

         if ( p1 == n+1 ):
            P[n+1:q] = 0

         P = reshape(P,(1,q))

      else:

         #  build Minimum Frobenius norm model (p1 <= q)
         #  minimizing the norm of the Hessian
         #  (QZ and RZ are the factors of F = [MQMQ' ML; ML' 0])

         #  If shifting is active, the matrix of the interpoaltion points is first
         #  shifted and scaled (to compute M correctly).

         if ( shift_Y ):
            xbase  = Y[:,0];
            scaleY = 0;
            for i in range(0,p1):
               Y[:,i] = Y[:,i] - xbase;
               scaleY = max( scaleY, linalg.norm( Y[:,i] ) );
            
            Y = Y / scaleY         

         # compute right-hand side with function values and zero

         rhs = concatenate( fY, zeros((1,n+1)) )


         mualpha = linalg.solve( RZ,dot(QZ.transpose(), rhs.transpose()))
#         mualpha = ( RZ \ ( QZ.transpose() * rhs.transpose() ) ).transpose();

         # constant and linear part of P

         P[1:n+1] = mualpha[p1+1:p1+n+1].transpose();

         # quadratic part of P

         M        = bcdfo_evalZ( Y, q ).transpose();
         P[n+1:q] = dot(M[:,n+1:q].transpose(), mualpha[0:p1].transpose());

   elif ( whichmodel == 2 ):

      # Builds minimum L2 norm model (p1 <= q)
      # (QZ and RZ are the factors of Z = M')
      # Take pseudo-inverse for underdetermined system because the result
      # is different from using backslash-operator

      P =  dot( QZ, linalg.lstsq(RZ.transpose(),fY))
#      P = ( QZ * ( RZ.transpose() \ fY.transpose() ) ).transpose();

   elif ( whichmodel == 3 ):

      #  Build Regression model (p1 >= q)
      #  (QZ and RZ are the factors of Z = M)
      #  Take pseudo-inverse for solving the system because the result
      #  is different from using backslash-operator (except for p1==q)


      #   P        = ( pinv(RZ) *  QZ' * fY' )';
      P = linalg.lstsq(RZ,QZ.transpose() * fY)
      #P = ( RZ \ ( QZ.transpose() * fY.transpose() ) ).transpose();

   elif ( whichmodel == 8 ):

      if ( shift_Y ):
         xbase  = Y[:,0];
         scaleY = 0;
         for i in range(0,p1):
            Y[:,i] = Y[:,i] - xbase;
            scaleY = max( scaleY, linalg.norm( Y[:,i] ) );

         Y     = Y / scaleY;

      Z = bcdfo_evalZ( Y, q );

      Pinter = linalg.lstsq(Z.transpose(),fY)
#      Pinter = ( Z.transpose() \ fY.transpose() ).transpose();
      P[1:n+1] = Pinter[1:n+1];

      Zg = [];
      g  = [];

      for i in range( 0, min( p1, 10*n )):
#      for i = 1

         Zg0 = diag( ones(n,1) );        # part of g

         Zg1 = diag( Y[:,i] );           # diagonal

         Zg2 = [];                       # the (i+1)-th subdiagonal
         jj=1;
         for k in range( 0,n-1):
            nsd = n-k;
            for j in range( 0,nsd):
               Zg2[j,jj] = Y[j,i];
               Zg2[k+j,jj] = Y[k+j,i];
               jj=jj+1;

         Zg = append( Zg, concatenate(Zg1, Zg2), axis = 0 );
#        Zg = [ Zg; Zg0 Zg1 Zg2 ];
#        Zg = [ Zg; Zg0 zeros(n,q-n-1) ];
         g  = append( g, gmat[indfree,ind_Y[i]] / scale[2])

         #  Adjusting the impact of the gradient on the resulting model

#        Zg = Zg .* impact_g;
#       g  = g .* impact_g;



#       Zg = [ Z(n+2:q,:)'; Zg ];
#        g  = [ fY'; g ];
      P[n+2:q] = linalg.lstsq(Zg,g)
#      P(n+2:q) = Zg \ g;

   elif ( whichmodel == 9 ):

      #  Build gmodel with gradient information into account
      # (QZ and RZ are the factors of Z = M')

      #  If shifting is active, the matrix of the interpoaltion points is first
      #  shifted and scaled (to compute Z correctly).

      if ( shift_Y ):
         xbase  = Y[:,1];
         scaleY = 0;
         for i in range(0,p1):
            Y[:,i] = Y[:,i] - xbase;
            scaleY = max( scaleY, linalg.norm( Y[:,i] ) );

         Y     = Y / scaleY;


      # compute standard system matrix

      Z = bcdfo_evalZ( Y, q );

      # build linear part of model

#     warning off
#     P(1:n+1) = ( Z(1:n+1,:)' \ fY' )';

      # build quadratic part of model

      Zg = [];
      g  = [];

      # compute impact of the gradient on the system

      if ( Delta0 > 1 ):
         impact_g = 3.0;
      else:
         impact_g = expm1(Delta0)**2;

      # compute error of the gradient and modify the impact on the system

      if ( n > 1 and p1 > n+1 ):
         K = min( p1 - 1, q - 1 );
         nneigh = p1 - 1;
         x0 = Y[:,0];
         f0 = fY[0];
         distance = sum( ( Y.transpose() - ones( p1, 1 ) * x0.transpose() )**2, 2 );
         [ dd, ind ] = sort( distance );
         ind = ind[2:K+1];
         d = max( abs( Y[:,1:p1-1].transpose() - ones( nneigh, 1 ) * x0.transpose() ), [], 1 );
         S = Y[:,ind].transpose();

         for i in range(0,K):
            Snorm[i]=linalg.norm(S[i,:]);

         near1=range(1,K)
         S1=S
         [ Sclosest, c ] = min( Snorm );
         while Snorm(c) < 1.0e-2 and size( near1, 2 ) > 1:
            near1[c] = [];
            S1[c,:] = [];
            Snorm[c] = [];
            [ Sclosest, c ] = min( Snorm );

         findiff = ( f0 - fY(near1(c)) ) / Sclosest;
         grad = ( gmat(indfree,ind_Y(1)).transpose() * S1[c,:].transpose() ) / Sclosest;
         ratio = abs( findiff - grad );
         qbar = S1[c,:] * diag( ones( size( gmat(indfree,ind_Y(1)), 2 ), 1 ) * ratio ) * S1[c,:].transpose();
         impact_g = impact_g / qbar;

#disp(['Delta: ',num2str(Delta0), ', ratio: ', num2str(ratio), ', qbar: ', num2str(qbar) , ', impact_g: ', num2str(impact_g)]);
#disp(' ');

#      else:

#disp(['Delta: ',num2str(Delta0), ', ratio: none, qbar: none, impact_g: ', num2str(impact_g)]);
#disp(' ');


      # loop over the gradients of all points / the gradient of the current iterate

#      for i in range(0,p1)
      for i in range(0):

         Zg0 = diag( ones((n,1),float) );        # part of g

         Zg1 = diag( Y[:,i] );           # diagonal

         Zg2 = [];                       # the (i+1)-th subdiagonal
         jj=1;
         for k in range(0,n-1):
            nsd = n-k;
            for j in range(0,nsd):
               Zg2[j,jj] = Y[j,i];
               Zg2[k+j,jj] = Y[k+j,i];
               jj=jj+1;


#         Zg = [ Zg; Zg1 Zg2 ];
         Zg = append( Zg, concatenate(Zg0, Zg1, Zg2),axis=0)
#         Zg = [ Zg; Zg0 zeros(n,q-n-1) ];
         g  = append( g, gmat[indfree,ind_Y(i)] / scale(2) )

         #  Adjusting the impact of the gradient on the resulting model

         Zg = Zg * impact_g;
         g  = g * impact_g;


      #  Enhance quadratic part by gradient information

#      Zg = [ Z(n+2:q,:)'; Zg ];
#      g  = [ fY'; g ];

#      P(n+2:q) = Zg \ g;
#      P(n+2:q) = pinv(Zg) * g;
#      P(n+2:q) = zeros(1,q-n-1);

       #  Compute linear and quadratic part at the same time.

#      P = ( Z' \ fY' )';

       #  Compute linear and quadratic part at the same time.
       #  Enhance quadratic part by exact gradient information

#      Zg = [ Z'; zeros(n*p1,1) Zg ];
      Z = append( Z.transpose(), concatenate(zeros((n,1)), Zg),axis=0)
      g  = g.insert(0, fY.transpose())

      P = linalg.lstsq(Z,g)
#      P = ( Zg \ g ).transpose();
#      P = (pinv(Zg) * g).transpose()

   return P
