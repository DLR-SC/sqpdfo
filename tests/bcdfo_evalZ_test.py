# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 13:35:53 2014

@author: jaco_da
"""

import unittest
from sqpdfo.bcdfo_evalZ import bcdfo_evalZ_
from numpy import *
from sqpdfo.runtime import  compare_array


class Test_bcdfo_evalZ(unittest.TestCase):
    """
    Reminder :
    This class is a test for evalZ which computes the matrix of the values of given points on the monomials basis 1, x1, x2,...,0.5*x1^2, 0.5*x2^2...x1*x2....
    """ 
    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-15;
        self.rel_tol=1e-15;
        pass

    
    def test_bcdfo_evalZ_1(self):
        """
        evaluation of the  6 following points on the monomial basis [1, x1, x2, 0.5*x1^2, 0.5*x2^2, x1*x2] :
        [ 0, 1, 0, 2, 1, 0], [0, 0, 1, 0, 0.01, 2 ]
        """
        res = bcdfo_evalZ_(array([[ 0, 1, 0, 2, 1, 0], [0, 0, 1, 0, 0.01, 2 ]]), 6)
        correctres = array([
    [1.0000,    1.0000,    1.0000,    1.0000,    1.0000,    1.0000],
     [    0,    1.0000,         0,    2.0000,    1.0000,         0],
      [   0,         0,    1.0000,         0,    0.0100,    2.0000],
       [  0,    0.5000,         0,    2.0000,    0.5000,         0],
        [ 0,         0,    0.5000,         0,    0.00005,    2.0000],
         [0,         0,         0,         0,    0.0100,         0]])
            
        self.assertTrue(compare_array(res, correctres, self.abs_tol, self.rel_tol))
        
    def test_bcdfo_evalZ_2(self):
        """
        evaluation of the  6 following points on the monomial basis [1, x1, x2, 0.5*x1^2, 0.5*x2^2, x1*x2] :
        [ 1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6 ]
        """
        res = bcdfo_evalZ_(array([[ 1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6 ]]), 6)
        correctres = array([
    [1.0000,    1.0000,    1.0000,    1.0000,    1.0000,    1.0000],
     [    1,         2,         3,         4,         5,         6],
      [   1,         2,         3,         4,         5,         6],
       [  0.5,         2,         4.5,         8,       12.5,         18],
        [ 0.5,         2,         4.5,         8,       12.5,         18],
         [1,         4,         9,         16,       25,         36]])
        self.assertTrue(compare_array(res, correctres, self.abs_tol, self.rel_tol))
        
    def test_bcdfo_evalZ_3(self):  
        """
        evaluation of the  6 following points on the monomial basis [1, x1, x2, 0.5*x1^2, 0.5*x2^2] :
        [ 0, 1, 0, 2, 1, 0], [0, 0, 1, 0, 0.01, 2 ]
        """
        res = bcdfo_evalZ_(array([[ 0, 1, 0, 2, 1, 0], [0, 0, 1, 0, 0.01, 2 ]]), 5)
        correctres = array([
    [1.0000,    1.0000,    1.0000,    1.0000,    1.0000,    1.0000],
     [    0,    1.0000,         0,    2.0000,    1.0000,         0],
      [   0,         0,    1.0000,         0,    0.0100,    2.0000],
       [  0,    0.5000,         0,    2.0000,    0.5000,         0],
        [ 0,         0,    0.5000,         0,    0.00005,    2.0000]]);
            
        self.assertTrue(compare_array(res, correctres, self.abs_tol, self.rel_tol))
        
    def test_bcdfo_evalZ_4(self):
        """
        evaluation of the  3 following points on the monoomial basis [1, x1, x2, 0.5*x1^2, 0.5*x2^2] :
        [ 0, 1, 0], [0, 0, 1]
        """
        res = bcdfo_evalZ_(array([[ 0, 1, 0], [0, 0, 1]]), 5)
        correctres = array([
    [1.0000,    1.0000,    1.0000],
     [    0,    1.0000,         0],
      [   0,         0,    1.0000],
       [  0,    0.5000,         0],
        [ 0,         0,    0.5000]]);
            
        self.assertTrue(compare_array(res, correctres, self.abs_tol, self.rel_tol))
        

if __name__ == '__main__':
    unittest.main()