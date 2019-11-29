# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 15:04:23 2014

@author: jaco_da
"""

import unittest
from sqpdfo.bcdfo_evalP import bcdfo_evalP_
from sqpdfo.bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
import numpy as np
from random import random


class Test_bcdfo_evalP(unittest.TestCase):
    """
 Reminder :
 This class is a test for evalP which computes the value of the model at x
    """ 

    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-15;
        self.rel_tol=1e-15;
        pass

    
    def test_bcdfo_evalP_1(self):
        """
        We check that the model on the interpolant set equals the given values of the initial function on the interpolant set, here we verify that the sixth point is interpolated
        """
        Y = np.array([[ 1, 2, 1, 3, 3, 1],[1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale  = bcdfo_build_QR_of_Y_( Y, 0, 1 , 1,1, 1e15)
        
        #we verify that the sixth point is interpolated
        model = ( QZ.dot(  np.linalg.solve( RZ.T , np.array([[1], [2], [3], [4], [5], [6] ]) ) )).T
        res = bcdfo_evalP_( model, np.array([[1],[3]]), xbase, scale, 1 )
        self.assertAlmostEqual(float(res), 6.0, places=12)
        
  
    def test_bcdfo_evalP_2(self):
        """
        We check that the model on the interpolant set equals the given values of the initial function on the interpolant set, here we verify that the sixth point is interpolated
        """
        Y = np.array([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        
        #we verify that the sixth point is interpolated
        model = ( QZ.dot(  np.linalg.solve( RZ.T , np.array([[random(),random(),random(),random(),random(), 6 ]]).T ) )).T
        res = bcdfo_evalP_( model, np.array([[1],[3]]), xbase, scale, 1 )
        self.assertAlmostEqual(float(res), 6, places=13) #NB : the random numbers forces us to reduce the precision
      
        
    def test_bcdfo_evalP_3(self):
        """
        We check that the model on the interpolant set equals the given values of the initial function on the interpolant set, here we verify that the sixth point is interpolated
        """
        Y = np.array([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        
        #we verify that the sixth point is interpolated
        model = ( QZ.dot(  np.linalg.solve( RZ.T , np.array([[0,0,0,0,0, 6 ]]).T ) )).T
        res = res = bcdfo_evalP_( model, np.array([[1],[3]]), xbase, scale, 1 )
        self.assertAlmostEqual(float(res), 6, places=14)
        
    def test_bcdfo_evalP_4(self):
        """
        We check that the model on the interpolant set equals the given values of the initial function on the interpolant set, here we verify that the first point is interpolated
        """
        Y = np.array([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        
        #we verify that the first point is interpolated
        model = ( QZ.dot(  np.linalg.solve( RZ.T , np.array([[6,0,0,0o3,0,0 ]]).T ) )).T
        res = bcdfo_evalP_( model, np.array([[1],[1]]), xbase, scale, 1 )
        self.assertAlmostEqual(float(res), 6, places=15)
        
    
    def test_bcdfo_evalP_5(self):
        """
        We check that in the very easy case where the fonction to interpolate is always equal to 1 on the interpolant set, then the model value on any
        x point will logically be 1
        """
        Y = np.array([[ 1, 2, 1, 3, 3, 1],[1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale  = bcdfo_build_QR_of_Y_( Y, 0, 1 , 1,1, 1e15)
        model = ( QZ.dot( np.linalg.solve( RZ.T , np.array([[1], [1], [1], [1], [1], [1] ]) ) )).T
        for i in range(0,50):
            res = bcdfo_evalP_( model, np.array([[(random()-0.5)*100],[(random()-0.5)*100]]), xbase, scale, 1 )
            self.assertAlmostEqual(float(res), 1.0, places=15)
        
        
if __name__ == '__main__':
    unittest.main()