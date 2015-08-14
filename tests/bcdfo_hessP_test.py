# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 12:24:11 2014

@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from runtime import matlabarray, compare_matlabarray
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from bcdfo_hessP import bcdfo_hessP_
import numpy as np
from random import random
#import helper

class Test_bcdfo_hessP(unittest.TestCase):
    """
      Reminder :
      This class is a test for HessP which computes the Hessian of P at X
    """ 
    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-13;
        self.rel_tol=1e-13;
        pass

    def test_bcdfo_hessP_1(self):
        """
           Easy tests, with models hand-given
        """
        for i in range(0, 20):
            x1=(random()-0.5)*100
            x2=(random()-0.5)*100
            model = matlabarray([(random()-0.5)*100, x1, x2, 0,0,0, ])
            ans = bcdfo_hessP_( model, matlabarray([[(random()-0.5)*100],[(random()-0.5)*100]]), matlabarray([[(random()-0.5)*100],[(random()-0.5)*100]]), 1, 0 )
            correctans = matlabarray([[ 0,0],[0,0]]).T
        
            self.assertTrue(compare_matlabarray(correctans, ans, self.abs_tol, self.rel_tol))
            
        for i in range(0, 20):
            x1=(random()-0.5)*100
            x2=(random()-0.5)*100
            x3=(random()-0.5)*100
            model = matlabarray([(random()-0.5)*100, (random()-0.5)*100, (random()-0.5)*100, x1,x2,x3, ])
            ans = bcdfo_hessP_( model, matlabarray([[(random()-0.5)*100],[(random()-0.5)*100]]), matlabarray([[(random()-0.5)*100],[(random()-0.5)*100]]), 1, 0 )
            correctans = matlabarray([[ x1,x3],[x3,x2]]).T
        
            self.assertTrue(compare_matlabarray(correctans, ans, self.abs_tol, self.rel_tol))

    def test_bcdfo_hessP_2(self):
        """
           Little bit more complicated tests initally written in the matlab code
        """
        Y = matlabarray([[ 1, 2, 1, 3, 3, 1],[ 1, 2, 2, 1, 2, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 0, 1, 1, 1e15 );
        model = ( QZ * np.linalg.solve( RZ.T , matlabarray([1, 2, 3, 4, 5, 6 ]).T ) ).T
        ans = bcdfo_hessP_( model, matlabarray([[0],[0]]) ,xbase, scale, 0 )
        
        correctans = matlabarray([[    4.0000,   -0.5000],
                        [-0.5000,    1.0000]])
        
        self.assertTrue(compare_matlabarray(correctans, ans, self.abs_tol, self.rel_tol))
        
        #Same test as above but with the shift in interpolation points
        Y = matlabarray([[ 1, 2, 1, 3, 3, 1],[ 1, 2, 2, 1, 2, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 );
        model = ( QZ * np.linalg.solve( RZ.T , matlabarray([1, 2, 3, 4, 5, 6 ]).T ) ).T
        ans = bcdfo_hessP_( model, matlabarray([[0],[0]]) ,xbase, scale, 1 )
        
        correctans = matlabarray([[    4.0000,   -0.5000],
                        [-0.5000,    1.0000]])
        
        self.assertTrue(compare_matlabarray(correctans, ans, self.abs_tol, self.rel_tol))
        
  
#  ans =
#
#    4.0000   -0.5000
#   -0.5000    1.0000

if __name__ == '__main__':
    unittest.main()