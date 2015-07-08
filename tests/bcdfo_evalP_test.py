# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 15:04:23 2014

@author: jaco_da
"""

import unittest
from bcdfo_evalP import bcdfo_evalP_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from runtime import matlabarray, compare_matlabarray
import numpy as np
from random import random

class Test_bcdfo_evalP(unittest.TestCase):

    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-15;
        self.rel_tol=1e-15;
        pass

    def test_bcdfo_evalP(self):
#        TEST:
        Y = matlabarray([[ 1, 2, 1, 3, 3, 1],[1, 2, 2, 1, 1.01, 3 ]]) 
        QZ, RZ, xbase, scale  = bcdfo_build_QR_of_Y_( Y, 0, 1 , 1,1, 1e15)
        model = ( QZ * np.linalg.solve( RZ.T , matlabarray([[1], [2], [3], [4], [5], [6] ]) ) ).T
        res = bcdfo_evalP_( model, matlabarray([[1],[3]]), xbase, scale, 1 )
    
        self.assertAlmostEqual(float(res), 6.0)
  
    def test_bcdfo_build_QR_of_Y_2(self):
        Y = matlabarray([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        model = ( QZ * np.linalg.solve( RZ.T , matlabarray([random(),random(),random(),random(),random(), 6 ]).T ) ).T
        res = bcdfo_evalP_( model, matlabarray([[1],[3]]), xbase, scale, 1 )
        
        self.assertTrue(compare_matlabarray(xbase, matlabarray([1,1]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(scale, matlabarray([1,0.499993750117185,0.499993750117185,0.249993750156246,0.249993750156246,0.249993750156246]), self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(float(res), 6)
        
    def test_bcdfo_build_QR_of_Y_3(self):
        Y = matlabarray([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        model = ( QZ * np.linalg.solve( RZ.T , matlabarray([0,0,0,0,0, 6 ]).T ) ).T
        res = res = bcdfo_evalP_( model, matlabarray([[1],[3]]), xbase, scale, 1 )

        self.assertAlmostEqual(float(res), 6)
        self.assertTrue(compare_matlabarray(xbase, matlabarray([1,1]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(scale, matlabarray([1,0.499993750117185,0.499993750117185,0.249993750156246,0.249993750156246,0.249993750156246]), self.abs_tol, self.rel_tol))
        
    def test_bcdfo_build_QR_of_Y_4(self):
        Y = matlabarray([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        model = ( QZ * np.linalg.solve( RZ.T , matlabarray([6,0,0,03,0,0 ]).T ) ).T
        res = bcdfo_evalP_( model, matlabarray([[1],[1]]), xbase, scale, 1 )

        self.assertAlmostEqual(float(res), 6)
        self.assertTrue(compare_matlabarray(xbase, matlabarray([1,1]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(scale, matlabarray([1,0.499993750117185,0.499993750117185,0.249993750156246,0.249993750156246,0.249993750156246]), self.abs_tol, self.rel_tol))


if __name__ == '__main__':
    unittest.main()