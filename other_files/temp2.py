# -*- coding: utf-8 -*-
"""
Created on Thu Jan 08 16:23:35 2015

@author: jaco_da
"""

import numpy as np
from runtime import matlabarray
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from bcdfo_checkZ import bcdfo_checkZ
from bcdfo_evalZ import bcdfo_evalZ
from bcdfo_computeLj import bcdfo_computeLj_

B = np.array([[1.0,1.0,1.0],[1.0, 77.0, 122.0],[1.0, 122.0, 194.0]])#[[14,32,50], [32, 77, 122], [50, 122, 194]])

#print B
Q, R = np.linalg.qr(B)

print "B\n", B
print "Q\n", Q
print "R\n", R

Y = matlabarray([[1.0, 77.0, 122.0], [1.0, 122.0, 194.0]])#, [50.0, 122.0, 194.0]])#matlabarray([[ 0, 1, 0, 2, 0], [0, 0, 1, 0, 2 ]])
whichmodel = 0
QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(  Y, whichmodel, 0, 1, 1, 1e15 )
print "QZ\n", QZ
print "RZ\n", RZ
Lj = bcdfo_computeLj_( QZ, RZ, 1, Y, whichmodel, scale, 1 )

#print "type(scale)", type(scale)
print "Lj", Lj
#correctLj = matlabarray([0.833333333333333,  -4.455569580721621,  -0.146502319909943])#matlabarray([   1.0000,  -3.0000,  -3.0000,   4.0000,   4.0000])		
#print "Warning: Different Results due to not unique QR decomposition (?)"
#self.assertEqual(Lj, correctLj)

#kappa_ill = 1e15

#Y = matlabarray([[1.0, 77.0, 122.0], [1.0, 122.0, 194.0]])
#npY = np.array([[1.0, 77.0, 122.0], [1.0, 122.0, 194.0]])
#QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 0, 1, 1, kappa_ill )

#Z = bcdfo_evalZ(npY, 3)
#print "Z", Z
#print "badcond =", bcdfo_checkZ( B, kappa_ill )

#print "QZ\n", QZ
#print "RZ\n", RZ


#Y = np.array([[ 1.0, 2.0, 1.0, 3.0, 3.0, 1.0], [1.0, 2.0, 2.0, 1.0, 1.01, 3.0 ]])
#Z = bcdfo_evalZ(Y, 3)