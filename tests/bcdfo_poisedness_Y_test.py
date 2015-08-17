# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 13:54:53 2014

@author: jaco_da
"""
import sys
sys.path.append("../")
sys.path.append("tests/")

import unittest
from bcdfo_poisedness_Y import *#bcdfo_poisedness_Y_
from bcdfo_build_QR_of_Y import *#bcdfo_build_QR_of_Y
import numpy as np
from numpy import array
#import helper

class Test_bcdfo_poisedness_Y(unittest.TestCase):
    """
      Reminder :
      This class is a test for bcdfo_replace_in_Y which omputes the poisedness of the interpolation set Y in a ball of radius Delta centered at Y(:,1), 
      assuming that a QR factorization of the  interpolation system is known (and given by QZ, RZ).
    """ 
    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()

        pass
    
    def test_bcdfo_poisedness_Y(self):
        """
        This is the test written in the Matlab Code. Results are the same.
        """
        Y = array([ [ 0, 1, 0, 2, 1, 0 ], [ 0, 0, 1, 0, 0.01, 2] ])
        #print Y
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )

        lSolver = 1
        whichmodel = 0
        hardcons = 0 
        stratLam = 1
        
        #print "Warning: Not enough input arguments"
        lambd ,Y_radius  = bcdfo_poisedness_Y_( QZ, RZ, Y, 0.001, xbase, lSolver, whichmodel, hardcons, None, None, None, stratLam,scale, 1)
        #print "lambd = ", lambd
        #print "Y_radius = ", Y_radius
        
        correctlambda = 2.048585522856004e+02
        correctY_radius = 2
        self.assertAlmostEqual(Y_radius, correctY_radius, places=15)
        self.assertAlmostEqual(lambd, correctlambda, places=11)
        
        #Same test without the shift  in interpolation points
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 0, 1, 1, 1e15 )
        lambd ,Y_radius  = bcdfo_poisedness_Y_( QZ, RZ, Y, 0.001, xbase, lSolver, whichmodel, hardcons, None, None, None, stratLam,scale, 0)
        self.assertAlmostEqual(Y_radius, correctY_radius, places=15)
        self.assertAlmostEqual(lambd, correctlambda, places=11)
        

if __name__ == '__main__':
    unittest.main()