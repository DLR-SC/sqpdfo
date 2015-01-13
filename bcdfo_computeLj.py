#!/usr/local/bin/python
from numpy import *
import helper

@helper.convertingDecorator
def bcdfo_computeLj_( QZ, RZ, j, Y, whichmodel, scale, shift_Y ):
	#j is an index, so we have to decrease it.
	return bcdfo_computeLj( QZ, RZ, j-1, Y, whichmodel, scale, shift_Y )

def bcdfo_computeLj( QZ, RZ, j, Y, whichmodel, scale, shift_Y ):
###############################################################################
#
#  Computes the coefficients of the j-th Lagrange polynomial.
#
#  INPUTS:
#
#  QZ, RZ      : the QR factors of the (possibly shifted) matrix Z containing
#                the polynomial expansion of the interpolation points
#  j           : the index of the interpolation point one wishes to compute
#                the associated Lagrange polynomial
#  Y           : current interpolation set Y
#  whichmodel  : the method to compute the polynomial
#  scale       : scaling factor of the interpolation matrix
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#
#  OUTPUT:
#
#  Lj          : a row vector containing the coefficients of the polynomial
#
#  DEPENDENCIES: -
#
#  PROGRAMMING: A. Troeltzsch, September 2010.
#
#  TEST:
#  Y = [ 0 1 0 2 0 ; 0 0 1 0 2 ];  whichmodel = 0;
#  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y(  Y, whichmodel, 1 );
#  Lj = bcdfo_computeLj( QZ, RZ, 1, Y, whichmodel, scale, 1 )
#  should give:
#  Lj =
#   1.0000  -3.0000  -3.0000   4.0000   4.0000
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
###############################################################################

    n, p1 = shape( Y )
    q = ( ( n + 1 ) * ( n + 2 ) ) / 2
    
    vzeros = zeros((j,1))
    vzeros = append(vzeros,1)
    vzeros = append(vzeros, zeros((p1-j-1,1)))

    #  as long as regression model is underdetermined: use min l2-norm model

    if ( whichmodel == 3 and p1 < q ):

        whichmodel = 2

    ###############################################################################

    if ( whichmodel == 0  ):

        # build (sub-basis) Lagrange polynomials (p1 = q)
        # (QZ and RZ are the factors of Z = M')
        
        Lj = reshape(dot(QZ, linalg.solve( RZ.transpose(), vzeros)).transpose(),(1,p1))
                      
#        P = reshape(dot(QZ, linalg.solve( RZ.transpose(),reshape(fY,(p1,1)))),(1,p1))

    elif ( whichmodel == 1 ):

        # build mixed Lagrange polynomials:
        # Minimum Frobenius norm Lagrange polynomials (p1 <= q) or
        # L2-norm Lagrange polynomials (p1 == n+1 or p1 == q )

        if ( p1 == n+1 or p1 == q ):

            # build linear/quadratic Lagrange polynomials

            Lj = linalg.solve( RZ,  dot(QZ.transpose(), vzeros) ).transpose()

            if ( p1 == n+1 ):
                Lj[n+1:q+1] = 0

            Lj = reshape(Lj,(1,q))

        else:

            # build Minimum Frobenius-norm Lagrange polynomials (p1 <= q)

            #  If shifting is active, the matrix of the interpoaltion points is first
            #  shifted and scaled (to compute M correctly).

            Y2 = copy( Y )
            if ( shift_Y ):
                xbase  = Y2[:,1]
                scaleY = 0
                for i in range(0,p1):
                    Y2[:,i] = Y2[:,i] - xbase
                    scaleY = max( scaleY, linalg.norm( Y2[:,i] ) )

            Y2     = Y2 / scaleY;

            # compute the right-hand side of the system

            rhs = zeros((j-1,1))
            rhs.append( 1 )
            rhs.append( zeros((p1+n+1-j,1)) )

            mualpha    = linalg.solve( RZ, dot( QZ.transpose(), rhs ) ).transpose()

            # constant and linear part of P

            Lj[0:n+1] = mualpha[p1:p1+n+1].transpose()

            # quadratic part of P

            M        = bcdfo_evalZ( Y2, q ).transpose();
            Lj[n+1:q] = dot(M[:,n+1:q].transpose(), mualpha[0:p1].transpose())     

    elif ( whichmodel == 2 ):

        # build Minimum L2 norm Lagrange polynomials (p1 <= q)
        # (QZ and RZ are the factors of Z = M')
        
        if ( p1 < q ):
            Lj = dot( QZ, linalg.lstsq( RZ.transpose(), vzeros) ).transpose()
        else:
            Lj = dot( QZ, linalg.solve( RZ.transpose(), vzeros) ).transpose()

    elif ( whichmodel == 3 ):

        # build Regression Lagrange polynomials (p1 >= q)
        # (QZ and RZ are the factors of Z = M)

        if ( p1 > q ):
            Lj = linalg.lstsq( RZ,  dot(QZ.transpose(), vzeros) ).transpose()
        else:
            Lj = linalg.solve( RZ,  dot(QZ.transpose(), vzeros) ).transpose()

    elif ( whichmodel == 8 ):

        if ( p1 < q ):
            Lj = dot( QZ, linalg.lstsq( RZ.transpose(), vzeros) ).transpose()
        else:
            Lj = dot( QZ, linalg.solve( RZ.transpose(), vzeros) ).transpose()

    elif ( whichmodel == 9 ):

        # build gmodel Lagrange polynomials
        # (for the interpolation matrix without the gradient)

        if ( p1 < q ):
            Lj = dot( QZ, linalg.lstsq( RZ.transpose(), vzeros) ).transpose()
        else:
            Lj = dot( QZ, linalg.solve( RZ.transpose(), vzeros) ).transpose()

    return Lj
