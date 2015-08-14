# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 14:41:23 2014

@author: jaco_da
"""
import sys
sys.path.append("../")

import unittest
from bcdfo_include_in_Y import bcdfo_include_in_Y_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from runtime import matlabarray, compare_matlabarray, char, ones_
from bcdfo_evalZ import *
import numpy as np
import helper

class Test_bcdfo_include_in_Y(unittest.TestCase):
    """
      Reminder :
      This class is a test for bcdfo_include_in_Y which attempts to include x in the interpolation set by replacing an existing
      interpolation point
    """ 
    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-15;
        self.rel_tol=1e-15;
        pass

    def test_bcdfo_include_in_Y(self):
#        TEST:
        """
            This is the test written in the Matlab Code. Results are the same except for a few signs due to a non-unique QR decomposition
        """
        Y = matlabarray([[ 0, 1, 0, 2, 1, 0],[0, 0, 1, 0, 0.01, 2 ]])
        whichmodel = 0
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, whichmodel, 1, 1, 1, 1e15 )
        QZplus, RZplus, Yplus, pos, xbase, scale = bcdfo_include_in_Y_(matlabarray([-1,1]).T, QZ, RZ, copy_(Y), matlabarray(range(2,6)), 0.01, 'weighted', xbase, whichmodel, 0, scale, 1, 1, 1, 1e15 )
        #print QZplus, RZplus, Yplus, pos, xbase, scale
        
        correctQZplus = matlabarray([
   [1.0000,         0,         0,         0,         0,         0],
   [     0,    -9.701425001453321e-01,    0.0000 ,  -2.425356250363330e-01,    0.0000,   -0.0000],
   [     0,    0.0000,   -9.701425001453321e-01,   -0.0000,    0.0000,   -2.425356250363330e-01],
   [     0,    -2.425356250363330e-01,   -0.0000,    9.701425001453319e-01,   -0.0000,    0.0000],
   [     0,    0.0000,   -2.425356250363330e-01,    0.0000,   -0.0000,    9.701425001453319e-01],
   [     0,   -0.0000,    0.0000,    0.0000,    1.0000,   0.0000]])

        correctRZplus = matlabarray([
   [1.0000,   1.0000,   1.0000,   1.0000,   1.0000,   1.0000],
   [      0,    -5.153882032022076e-01,         0,    -1.091410312663498e+00,   4.547542969431244e-01,    0.0000],
   [      0,         0,   -5.153882032022076e-01,   -0.0000,   -5.153882032022077e-01,   -1.091410312663498e+00],
   [      0,         0,         0,    2.425356250363330e-01,    2.425356250363330e-01,   -0.0000],
   [      0,         0,         0,         0,   -0.2500,   -0.0000],
   [      0,         0,         0,         0,         0,    2.425356250363330e-01]])

        correctYplus = matlabarray([

     [0,     1,     0,     2,    -1,     0],
     [0,     0,     1,     0,     1,     2]])

        correctpos = 5
        correctxbase = matlabarray([[0],[0]])
        correctscale = matlabarray([ 
    [1.0000],
    [0.5000],
    [0.5000],
    [0.2500],
    [0.2500],
    [0.2500],])
                
        self.assertTrue(compare_matlabarray(correctQZplus, QZplus, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctRZplus, RZplus, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctYplus, Yplus, self.abs_tol, self.rel_tol))
        self.assertEqual(correctpos, pos)
        self.assertAlmostEqual(correctxbase, xbase)
        self.assertAlmostEqual(correctscale, scale)
        var=np.dot(QZplus.T, bcdfo_evalZ_((Yplus-np.dot(Y[:,1], ones_(1,6)))*scale[2],6))
        ident=np.linalg.solve(var,RZplus)
        self.assertTrue(compare_matlabarray(matlabarray(ident), matlabarray(np.eye(6)),self.abs_tol, self.rel_tol))

       #  RZplus\(QZplus'*bcdfo_evalZ((Yplus-Y(:,1)*ones(1,6))*scale(2),6))
       #  should gives the identity matrix.

if __name__ == '__main__':
    unittest.main()