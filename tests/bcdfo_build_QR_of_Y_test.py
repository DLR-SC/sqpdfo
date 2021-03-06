

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 11:58:39 2014

@author: jaco_da
"""

import unittest
from sqpdfo.bcdfo_build_QR_of_Y import *
from sqpdfo.bcdfo_evalZ import *
from numpy import *
from sqpdfo.runtime import *
from random import random
from sqpdfo.runtime import compare_array

class Test_bcdfo_build_QR_of_Y(unittest.TestCase):
    """
    Reminder :
    This class is a test for build_QR_of_Y which gives us not only the QR and RZ matrices, but only xbase and scale
    """ 
    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-15;
        self.rel_tol=1e-15;
        pass

#Tests with whichmodel=0 and an active shift
    def test_bcdfo_build_QR_of_Y_0(self):
        """
        Tests with whichmodel=0 and an active shift : we verifiy that xbase and scale are OK
        """
        Y = array([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        self.assertTrue(compare_array(xbase, array([1,1]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(scale, array([1,0.499993750117185,0.499993750117185,0.249993750156246,0.249993750156246,0.249993750156246]), self.abs_tol, self.rel_tol))
        Y = array([[ 0, 1, 0, 2, 1, 0],[0, 0, 1, 0, 0.01, 2 ]]) 
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1,1,1e15 )
        self.assertTrue(compare_array(xbase, array([0,0]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(scale, array([1,0.5, 0.5, 0.25, 0.25, 0.25]), self.abs_tol, self.rel_tol))
        
    def test_bcdfo_build_QR_of_Y_1(self):
        """
        Tests with whichmodel=0 and an active shift : we verifiy that the model value calculated with QZ and RZ is OK
        """
        Y = array([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        model = ( QZ.dot(np.linalg.solve( RZ.T , array([[1, 2, 3, 4, 5, 6 ]]).T ) )).T
        res =  model.dot(bcdfo_evalZ_( (array([[1],[3]])-xbase)*scale[1],6))
        self.assertTrue(compare_array(scale, array([1,0.499993750117185,0.499993750117185,0.249993750156246,0.249993750156246,0.249993750156246]), self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(float(res), 6,places=12)
        
    def test_bcdfo_build_QR_of_Y_2(self):
        """
        Tests with whichmodel=0 and an active shift : we verifiy that the model value calculated with QZ and RZ is OK
        """
        Y = array([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        for i in range(0,100):
            model = ( QZ .dot(np.linalg.solve( RZ.T , array([[random(),random(),random(),random(),random(), 6 ]]).T )) ).T
            res = np.dot(model, bcdfo_evalZ_( (array([[1],[3]])-xbase)*scale[1],6))
            self.assertAlmostEqual(float(res), 6,places=12)
        array([1,3])
    def test_bcdfo_build_QR_of_Y_3(self):       
        """
        Tests with whichmodel=0 and an active shift : we verifiy that the model value calculated with QZ and RZ is OK
        """
        Y = array([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        model = ( QZ.dot( np.linalg.solve( RZ.T , array([[0,0,0,0,0, 6 ]]).T ) )).T
        res = np.dot(model, bcdfo_evalZ_( (array([[1],[3]])-xbase)*scale[1],6))

        self.assertAlmostEqual(float(res), 6,places=14)
        
    def test_bcdfo_build_QR_of_Y_4(self):
        """
        Tests with whichmodel=0 and an active shift : we verifiy that the model value calculated with QZ and RZ is OK
        """
        Y = array([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
        model = ( QZ.dot( np.linalg.solve( RZ.T , array([[6,0,0,0o3,0,0 ]]).T ))).T
        res = np.dot(model, bcdfo_evalZ_( (array([[1],[1]])-xbase)*scale[1],6))

        self.assertAlmostEqual(float(res), 6,places=15)

if __name__ == '__main__':
    unittest.main()

