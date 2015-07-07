

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 11:58:39 2014

@author: jaco_da
"""

import unittest
from bcdfo_build_QR_of_Y import *
import numpy as np
from runtime import *
import helper
from random import random
from runtime import matlabarray, compare_matlabarray

class Test_bcdfo_build_QR_of_Y(unittest.TestCase):

    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-15;
        self.rel_tol=1e-15;
        pass

#Tests with whichmodel=0 and an active shift

    def test_bcdfo_build_QR_of_Y(self):
        Y = matlabarray([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        model = ( QZ * np.linalg.solve( RZ.T , matlabarray([1, 2, 3, 4, 5, 6 ]).T ) ).T
        res = np.dot(helper.convert(model) , bcdfo_evalZ( (np.array([[1,3]]).T-helper.convert(xbase))*helper.convert(scale[2]),6))
        
        self.assertTrue(compare_matlabarray(xbase, matlabarray([1,1]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(scale, matlabarray([1,0.499993750117185,0.499993750117185,0.249993750156246,0.249993750156246,0.249993750156246]), self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(float(res), 6)
        
    def test_bcdfo_build_QR_of_Y_2(self):
        Y = matlabarray([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        model = ( QZ * np.linalg.solve( RZ.T , matlabarray([random(),random(),random(),random(),random(), 6 ]).T ) ).T
        res = np.dot(helper.convert(model) , bcdfo_evalZ( (np.array([[1,3]]).T-helper.convert(xbase))*helper.convert(scale[2]),6))
        
        self.assertTrue(compare_matlabarray(xbase, matlabarray([1,1]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(scale, matlabarray([1,0.499993750117185,0.499993750117185,0.249993750156246,0.249993750156246,0.249993750156246]), self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(float(res), 6)
        
    def test_bcdfo_build_QR_of_Y_3(self):
        Y = matlabarray([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        model = ( QZ * np.linalg.solve( RZ.T , matlabarray([0,0,0,0,0, 6 ]).T ) ).T
        res = np.dot(helper.convert(model) , bcdfo_evalZ( (np.array([[1,3]]).T-helper.convert(xbase))*helper.convert(scale[2]),6))

        self.assertAlmostEqual(float(res), 6)
        self.assertTrue(compare_matlabarray(xbase, matlabarray([1,1]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(scale, matlabarray([1,0.499993750117185,0.499993750117185,0.249993750156246,0.249993750156246,0.249993750156246]), self.abs_tol, self.rel_tol))
        
    def test_bcdfo_build_QR_of_Y_4(self):
        Y = matlabarray([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        model = ( QZ * np.linalg.solve( RZ.T , matlabarray([6,0,0,03,0,0 ]).T ) ).T
        res = np.dot(helper.convert(model) , bcdfo_evalZ( (np.array([[1,1]]).T-helper.convert(xbase))*helper.convert(scale[2]),6))

        self.assertAlmostEqual(float(res), 6)
        self.assertTrue(compare_matlabarray(xbase, matlabarray([1,1]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(scale, matlabarray([1,0.499993750117185,0.499993750117185,0.249993750156246,0.249993750156246,0.249993750156246]), self.abs_tol, self.rel_tol))

if __name__ == '__main__':
    unittest.main()

