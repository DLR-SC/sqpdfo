# -*- coding: utf-8 -*-
"""
Created on Mon Dec 01 11:12:48 2014

@author: jaco_da
"""
import sys
sys.path.append("../")
sys.path.append("tests/")

import unittest
from bcdfo_augment_Y import bcdfo_augment_Y_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from runtime import matlabarray, compare_matlabarray
#import numpy as np
#import helper

class Test_bcdfo_augment_Y(unittest.TestCase):
    """
      Reminder :
      This class is a test for bcdfo_include_in_Y which attempts to include x in the interpolation set by replacing an existing
      interpolation point
    """ 
    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-14;
        self.rel_tol=1e-14;
        pass

    def test_bcdfo_augment_Y(self):  
        """
            This is the test written in the Matlab Code. Results are the same except for a few signs due to a non-unique QR decomposition
        """
        Y = matlabarray([[ 0, 1, 0, 2, 0], [0, 0, 1, 0, 2 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(  Y, 0, 0, 1, 1, 1e15 )
        p1, QZ, RZ, Y, xbase, scale = bcdfo_augment_Y_( matlabarray([1, 0.01]).T, Y, 0, 0, 1, 1, 1e15 )
        #print "p1, QZ, RZ, Y, xbase, scale", p1, QZ, RZ, Y, xbase, scale
        
        correctQZ=matlabarray([
        [-1,	0,	0,	0,	0,	0],
[0,	-0.894427190999916,	0,	0.447213595499958,	0,	0],
[0,	0,	-0.894427190999916,	0,	0.447213595499958,	0],
[0,	-0.447213595499958,	0,	-0.894427190999916,	0,	0],
[0,	0,	-0.447213595499958,	0,	-0.894427190999916, 	0],
[0,	0,	0,	0,	0,	-1]])

        correctRZ=matlabarray([
[-1,	-1,	-1,	-1,	-1,	-1],
[0,	-1.11803398874990	, 0,	-2.68328157299975, 	0,	-1.11803398874990],
[0,	0,	-1.11803398874990,	0,	-2.68328157299975,	-0.00896663258977416],
[0,	0,	0,	-0.894427190999916	, 0,	0],
[0,	0,	0,	0,	-0.894427190999916,	0.00442741459544958],
[0,	0,	0,	0,	0,	-0.0100000000000000]])

        
        correctY = matlabarray([
         [0,    1.0000,         0,    2.0000,         0,    1.0000],
         [0,         0,    1.0000,         0,    2.0000,    0.0100]])
        
        self.assertEqual(p1, 6)
        self.assertTrue(compare_matlabarray(correctQZ, QZ, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctRZ, RZ, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctY, Y, self.abs_tol, self.rel_tol))

        Y = matlabarray([[ 0, 1, 0, 2, 0], [0, 0, 1, 0, 2 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(  Y, 0, 1, 1, 1, 1e15 )
        p1, QZ, RZ, Y, xbase, scale = bcdfo_augment_Y_( matlabarray([1, 0.01]).T, Y, 0, 1, 1, 1, 1e15 )
        #print "p1, QZ, RZ, Y, xbase, scale", p1, QZ, RZ, Y, xbase, scale
        
        correctQZ=matlabarray([
        [-1,	0,	0,	0,	0,	0],
[0,	-9.701425001453321e-01,	0,	2.425356250363330e-01,	0,	0],
[0,	0,	-9.701425001453321e-01,	0,	2.425356250363330e-01,	0],
[0,	-2.425356250363330e-01,	0,	-9.701425001453319e-01,	0,	0],
[0,	0,	-2.425356250363330e-01,	0,	-9.701425001453319e-01, 	0],
[0,	0,	0,	0,	0,	-1]])

        correctRZ=matlabarray([
[-1,	-1,	-1,	-1,	-1,	-1],
[0,	-5.153882032022076e-01, 0,	-1.091410312663498, 	0,	-0.515388203202208],
[0,	0,	-5.153882032022076e-01,	0,	-1.091410312663498,	-0.00485374419603962],
[0,	0,	0,	-2.425356250363330e-01	, 0,	0],
[0,	0,	0,	0,	-2.425356250363330e-01,	0.00120055134392985],
[0,	0,	0,	0,	0,	-0.00250000000000000]])

        
        correctY = matlabarray([
         [0,    1.0000,         0,    2.0000,         0,    1.0000],
         [0,         0,    1.0000,         0,    2.0000,    0.0100]])
        
        self.assertEqual(p1, 6)
        self.assertTrue(compare_matlabarray(correctQZ, QZ, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctRZ, RZ, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctY, Y, self.abs_tol, self.rel_tol))


if __name__ == '__main__':
    unittest.main()