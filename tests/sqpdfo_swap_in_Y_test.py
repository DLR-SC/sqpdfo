# -*- coding: utf-8 -*-
"""

Created on Fri Nov 14 18:53:38 2014
% TESTS :

%  Y = [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; ind_Y = [1 2 3 4 5 6];
%  fY = [2 5 1 3 2 6];
%  [QZ,RZ,xbase,scale] = bcdfo_build_QR_of_Y( Y, 0, 0, 1, 1, 1e15 );
%  Z = QZ*RZ;
%  [ QZ, RZ, Y, xbase, scale ] = sqpdfo_swap_in_Y( 1, 3, QZ, RZ, Y, ind_Y, ...
%     fY, [],[],xbase, 0, scale, 0, 1, 1, 1e15 )
%  [ QZ, RZ, Y, xbase, scale ] = sqpdfo_swap_in_Y( 1, 3, QZ, RZ, Y, ind_Y, ...
%     fY, [],[],xbase, 0, scale, 0, 1, 1, 1e15 )
%  norm( Z - QZ*RZ)
%  should give something very small.
%
%  The same holds for the scaled version:
%  Y = [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; ind_Y = [1 2 3 4 5 6];
%  fY = [2 5 1 3 2 6];
%  [QZ,RZ,xbase,scale] = bcdfo_build_QR_of_Y( Y, 0, 1, 1, 1, 1e15 );
%  Z = QZ*RZ;
%  [ QZ, RZ, Y, xbase, scale ] = sqpdfo_swap_in_Y( 1, 3, QZ, RZ, Y, ind_Y, ...
%     fY, [], [], xbase, 0, scale, 1, 1, 1, 1e15 )
%  [ QZ, RZ, Y, xbase, scale ] = sqpdfo_swap_in_Y( 1, 3, QZ, RZ, Y, ind_Y, ...
%     fY, [],[],xbase, 0, scale, 1, 1, 1, 1e15 )
%  norm( Z - QZ*RZ)
@author: jaco_da
"""

import unittest
from sqpdfo.sqpdfo_swap_in_Y import *
from sqpdfo.bcdfo_build_QR_of_Y import *
from numpy import array


class Test_sqpdfo_swap_in_Y(unittest.TestCase):
    """
      Reminder :
      This class is a test for sqpdfo_swap_in_Y which swaps the position of interpolation points i and j in Y, and updates the
     factorization of Z(Y) accordingly.
    """ 
    def setUp(self):
        self.Y = array([[ 0, 1, 0, 2, 1, 0], [0, 0, 1, 0, 0.01, 2 ]])
        self.ind_Y = array([0, 1, 2, 3, 4, 5])
        self.fY = array([2, 5, 1, 3, 2, 6])
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()

        pass

    def test_sqpdfo_swap_in_Y(self):
        
        """
            Matlab test
        """
        QZ,RZ,xbase,scale = bcdfo_build_QR_of_Y_( self.Y, 0, 0, 1, 1, 1e15 )
        Z = QZ.dot(RZ)
#        
        ceY = array([])
        ciY = array([])
        QZ,RZ,Y,ind_Y,fY,ciY,ceY,xbase,scale = sqpdfo_swap_in_Y_( 0, 2, QZ, RZ, self.Y, self.ind_Y,
        self.fY, ciY, ceY, xbase, 0, scale, 0, 1, 1, 1e15 )
        QZ,RZ,Y,ind_Y,fY,ciY,ceY,xbase,scale = sqpdfo_swap_in_Y_( 0, 2, QZ, RZ, Y, ind_Y,
        fY, ciY, ceY, xbase, 0, scale, 0, 1, 1, 1e15 )
        ret = norm_( Z - QZ.dot(RZ))  
        self.assertTrue(ret < 1e-15)
#        
#        #Scaled version :
        QZ,RZ,xbase,scale = bcdfo_build_QR_of_Y_( self.Y, 0, 1, 1, 1, 1e15 )
        Z = QZ.dot(RZ)
        ceY = array([])
        ciY = array([])
        QZ,RZ,Y,ind_Y,fY,ciY,ceY,xbase,scale = sqpdfo_swap_in_Y_( 0, 2, QZ, RZ, self.Y, self.ind_Y,
        self.fY, ciY, ceY, xbase, 0, scale, 1, 1, 1, 1e15 )
        QZ,RZ,Y,ind_Y,fY,ciY,ceY,xbase,scale = sqpdfo_swap_in_Y_( 0, 2, QZ, RZ, Y, ind_Y,
        fY, ciY, ceY, xbase, 0, scale, 1, 1, 1, 1e15 )
        ret = norm_( Z - QZ.dot(RZ))  
        self.assertTrue(ret < 1e-15)

if __name__ == '__main__':
    unittest.main()