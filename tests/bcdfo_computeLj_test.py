# -*- coding: utf-8 -*-
"""
Created on Mon Dec 01 13:41:32 2014

@author: jaco_da
"""

import unittest
from bcdfo_computeLj import *#bcdfo_computeLj_
from bcdfo_build_QR_of_Y import *#bcdfo_build_QR_of_Y_
from runtime import matlabarray, compare_matlabarray
from bcdfo_evalP import bcdfo_evalP_
import numpy as np
from random import random
#import helper

class Test_bcdfo_computeLj(unittest.TestCase):
    """
       Reminder :
      This class is a test for computeLj which computes the coefficients of the j-th Lagrange polynomial
    """ 
    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-14;
        self.rel_tol=1e-14;
        pass

    def test_bcdfo_computeLj_1(self):
        """
        This test verify the basic properties of the lagrange polynomials : We verify here that for x being a point of the interpolant set, we have l_j(x) = 1 if the index j corresponds to the index of x in the interpolant set
        and 0 otherwise
        """
        for k in range(0,20):
            x1_coord=[(random()-0.5)*100 for p in range(0,5)]
            x2_coord=[(random()-0.5)*100 for p in range(0,5)]
            Y = matlabarray([x1_coord, x2_coord])
            whichmodel = 0
            QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(  Y, whichmodel, 1, 1, 1, 1e15 )
            for i in range(1,6):
                L= bcdfo_computeLj_( QZ, RZ, i, Y, whichmodel, scale, 1 )
                for j in range(0,5):
                    value=bcdfo_evalP_( L, matlabarray([[x1_coord[j]],[x2_coord[j]]]), xbase, scale, 1 )
                    if i==j+1: # we have to be careful with the indices
                        self.assertAlmostEqual(value, 1.0,places=10)
                    else:
                        self.assertAlmostEqual(value, 0,places=10)

    def test_bcdfo_computeLj_2(self):
        """
        More complicated test initally written in the matlab code
        """
        Y = matlabarray([[ 0, 1, 0, 2, 0], [0, 0, 1, 0, 2 ]])
        whichmodel = 0
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(  Y, whichmodel, 1, 1, 1, 1e15 )
        L1 = bcdfo_computeLj_( QZ, RZ, 1, Y, whichmodel, scale, 1 )
        correctL1 = matlabarray([1.0000,  -3.0000,  -3.0000,   4.0000,   4.0000])
        self.assertTrue(compare_matlabarray(correctL1,L1, self.abs_tol, self.rel_tol))
        


if __name__ == '__main__':
    unittest.main()