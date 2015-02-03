# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 15:58:53 2014

@author: jaco_da
"""
import numpy as np
from bcdfo_evalZ import *

import helper
from runtime import matlabarray

#@helper.convertingDecorator
#def bcdfo_computeP_(QZ, RZ, Y, fY, whichmodel, P_old, ind_Y,
#         i_xold, i_xplus, g, scale, shift_Y, Delta0=None, indfree=None, gmat=None, hessian=None, 
#         epsilon=None, noisy=None):
#										
#    return bcdfo_computeP(QZ, RZ, Y, fY, whichmodel, P_old, ind_Y,
#         i_xold, i_xplus, g, scale, shift_Y, Delta0=None)#, indfree=None, gmat=None, hessian=None, 
#         #epsilon=None, noisy=None)

def bcdfo_computeP_(QZ, RZ, Y, fY, whichmodel, P_old, ind_Y,
         i_xold, i_xplus, g, scale, shift_Y, Delta0=None, indfree=None, gmat=None, hessian=None, 
         epsilon=None, noisy=None):
		
		QZ = helper.convert(QZ)
		RZ = helper.convert(RZ)
		Y = helper.convert(Y)
		#print "fY", fY
		#print "fY[1,:]", fY[1,:]
		#print "type fY1,:", type(fY[1,:])
		#for lst in fY:
	#		break
		#lst = [gnlph for gnlph in fY[1, :]]
	#	print "gnlph lst:", lst
		fY = helper.convertArray(fY)#lst)
		whichmodel = helper.convert(whichmodel)
		P_old = helper.convert(P_old)
		ind_Y = helper.convert(ind_Y)
		i_xold = helper.convert(i_xold)
		i_xplus = helper.convert(i_xplus)
		g = helper.convert(g)
		scale = helper.convert(scale)
		shift_Y = helper.convert(shift_Y)
		Delta0 = helper.convert(Delta0)
		indfree = helper.convert(indfree)
		gmat = helper.convert(gmat)
		hessian = helper.convert(hessian) 
		epsilon = helper.convert(epsilon)
		noisy = helper.convert(noisy)
		
		P = bcdfo_computeP(QZ, RZ, Y, fY, whichmodel, P_old, ind_Y,
			i_xold, i_xplus, g, scale, shift_Y, Delta0)#), indfree, gmat, hessian, 
		      #epsilon, noisy)
				
		P = matlabarray(P)
						
		return P

#return P


def bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, P_old, ind_Y, i_xold, i_xplus, g, scale, shift_Y, Delta0 ):
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%
#%  Computes the polynomial P, where P is then represented by a row vector
#%  containing its coefficients for the successive monomials.
#%  More specifically, these values are:
#%  P(1)        : constant coefficient,
#%  P(2:n+1)    : coefficients for the linear terms in x(1)... x(n),
#%  P(n+2:2n+2) : coefficients for the squared terms in x(1)^2 ... x(n)^2
#%  P(2n+3,3n+2): coefficients for the quadratic terms of the first subdiagonal:
#%                in x(1)*x(2) ... x(n-1)*x(n)
#%  P(3n+3,4n+1): coefficients for the quadratic terms of the second subdiagonal:
#%                in x(1)*x(3) ... x(n-2)*x(n)
#%  etc.
#%
#%  INPUTS:
#%
#%  QZ, RZ      : the QR factors of the (possibly shifted) matrix Z containing
#%                the polynomial expansion of the interpolation points
#%  Y           : current interpolation set
#%  fY          : function values of the interpolation points
#%  whichmodel  : the kind of model to build
#%  P_old       : polynomial coefficients of the former polynomial
#%  scale       : scaling of the interpolation matrix
#%  shift_Y     : shift of the matrix of interpolation points
#%
#%  OUTPUT:
#%
#%  P           : a row vector containing the coefficients of the polynomial
#%
#%  DEPENDENCIES: -
#%
#%  PROGRAMMING: A. Troeltzsch, September 2010.
#%
#%  TEST:
#%  Y = [ 0 1 0 2 0 ; 0 0 1 0 2 ]; fY = [ 1 2 3 4 5 ];
#%  whichmodel = 0;
#%  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y(  Y, whichmodel, 1 );
#%  P = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, [0, 0, 0, 0, 0], 0, ...
#%      1, 1, [0, 0], scale, 1, 1 )
#%  should give:
#%  P =
#%   1.0000   1.0000   4.0000   4.0000       0
#%
#%  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[n,p1] = Y.shape#size(Y)
	#badcond = 0
	q = ( ( n + 1 ) * ( n + 2 ) ) / 2

	if ( whichmodel == 3 and p1 < q ):
		#  %  for underdetermined regression model use min l2-norm model
		whichmodel = 2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if ( whichmodel == 0 ):

		   #% build (sub-basis) model (p1 = q) 
		   #% (QZ and RZ are the factors of Z = M')

		   #warning off

		P        = np.dot( QZ , np.linalg.solve( RZ.T , fY.T ) ).T

	elif ( whichmodel == 1 ):

		#%  build mixed model: Minimum Frobenius norm model (p1 <= q)
		#%  and l2-norm model (p1 == n+1 or p1 == q)
   
		n_rhs = fY.shape[0]#size(fY,1)

		if ( p1 == n+1 or p1 == q ):
      
			P[0:n_rhs-1,0:p1-1] = np.linalg.solve( RZ , np.dot( QZ.T , fY.T ) ).T

			if ( p1 == n+1 ):
				P[0:n_rhs-1,n+1:q-1] = 0.0
		else:

			#%  build Minimum Frobenius norm model (p1 <= q)
			#%  minimizing the norm of the Hessian
			#%  (QZ and RZ are the factors of F = [MQMQ' ML; ML' 0])
	
			#%  If shifting is active, the matrix of the interpoaltion points is first
			#%  shifted and scaled (to compute M correctly).

			if ( shift_Y ):
				xbase  = Y[:,0]
				scaleY = 0
				for i in range(p1):
					Y[:,i] = Y[:,i] - xbase
					scaleY = max( scaleY, np.linalg.norm( Y[:,i] ) )

				Y     = Y / scaleY
      

	      #% compute right-hand side with function values and zero
			P = []
			for i in range(n_rhs):
				rhs = fY[i,:].append( np.zeros((1,n+1)))

				mualpha    = np.linalg.solve( RZ , np.dot( QZ.T , rhs.T ) ).T

				#% constant and linear part of P

				P_i[0:n] = mualpha[p1:p1+n].T

				#% quadratic part of P

				M        = bcdfo_evalZ( Y, q ).T
				P_i[n+1:q-1] = np.dot(M[ :, n+2 : q ].T , mualpha[ 1: p1 ].T)
         
				P = P.append(P_i, axis=1)
         
	elif ( whichmodel == 2 ):
		

		#% Minimum L2 norm model (p1 <= q)
		#% (QZ and RZ are the factors of Z = M')
		#% Take pseudo-inverse for underdetermined system because the result
		#% is different from using backslash-operator
   
#%   if ( p1 < q )

#%      warning off
#%      P        = ( QZ * ( pinv(RZ') * fY' ) )';

#%   else

#      warning off
		P        = np.dot( QZ , np.linalg.solve( RZ.T , fY.T ) ).T

#%   end
   
	elif ( whichmodel == 3 ):

		#% build Regression model (p1 >= q)
		#% (QZ and RZ are the factors of Z = M)
		#% Take pseudo-inverse for solving the system because the result
		#% is different from using backslash-operator (except for p1==q)

#   warning off
		P        = np.dot(np.dot( np.linalg.pinv(RZ) ,  QZ.T) , fY.T ).T


	return P