# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 11:59:20 2014

@author: jaco_da
"""

import unittest
from bcdfo_gradP import bcdfo_gradP_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from runtime import matlabarray, compare_matlabarray
from random import random
import numpy as np
#import helper

class Test_bcdfo_gradP(unittest.TestCase):
    """
      Reminder :
      This class is a test for gradP which computes the gradient of P at X
    """ 
    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-13;
        self.rel_tol=1e-13;
        pass

    def test_bcdfo_gradP_1(self):
        """
           Easy tests, with models hand-given
        """
        for i in range(0, 20):
            x1=(random()-0.5)*100
            x2=(random()-0.5)*100
            model = matlabarray([0, x1, x2, 0,0,0, ])
            ans = bcdfo_gradP_( model, matlabarray([[(random()-0.5)*100],[(random()-0.5)*100]]), matlabarray([[(random()-0.5)*100],[(random()-0.5)*100]]), 1, 0 )
            correctans = matlabarray([ x1, x2]).T
        
            self.assertTrue(compare_matlabarray(correctans, ans, self.abs_tol, self.rel_tol))
            
        for i in range(0, 20):
            x1=(random()-0.5)*100
            x2=(random()-0.5)*100
            x3=(random()-0.5)*100
            x4=(random()-0.5)*100
            model = matlabarray([0, 0, 0, x1,x2,0, ])
            ans = bcdfo_gradP_( model, matlabarray([[x3],[x4]]), matlabarray([[(random()-0.5)*100],[(random()-0.5)*100]]), 1, 0 )
            correctans = matlabarray([ x1*x3, x2*x4]).T
        
            self.assertTrue(compare_matlabarray(correctans, ans, self.abs_tol, self.rel_tol))

    def test_bcdfo_gradP_2(self):
        """
           Little bit more complicated tests initally written in the matlab code
        """
        Y = matlabarray([[ 1, 2, 1, 3, 3, 1],[1, 2, 2, 1, 2, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 0, 1,1, 1e15 )
        model = ( QZ * np.linalg.solve( RZ.T , matlabarray([1, 2, 3, 4, 5, 6 ]).T ) ).T
        ans = bcdfo_gradP_( model, matlabarray([[0],[0]]), xbase, scale, 0 )
        correctans = matlabarray([ -6.0000, 1.0000]).T
        
        self.assertTrue(compare_matlabarray(correctans, ans, self.abs_tol, self.rel_tol))
        
        #Same test as above but with the shift in interpolation points
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1,1, 1e15 )
        model = ( QZ * np.linalg.solve( RZ.T , matlabarray([1, 2, 3, 4, 5, 6 ]).T ) ).T
        ans = bcdfo_gradP_( model, matlabarray([[0],[0]]), xbase, scale, 1 )
        correctans = matlabarray([ -6.0000, 1.0000]).T
        
        self.assertTrue(compare_matlabarray(correctans, ans, self.abs_tol, self.rel_tol))


if __name__ == '__main__':
    unittest.main()