# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 13:35:53 2014

@author: jaco_da
"""
import sys
sys.path.append("../")
 
import unittest
from bcdfo_evalZ import bcdfo_evalZ_
import numpy as np
from runtime import matlabarray, compare_matlabarray
#import helper

class Test_bcdfo_evalZ(unittest.TestCase):

    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-15;
        self.rel_tol=1e-15;
        pass

    #evaluation of the  6 following points on the polynomial basis [1, x1, x2, 0.5*x1^2, 0.5*x2^2, x1*x2]
    def test_bcdfo_evalZ_1(self):
        res = bcdfo_evalZ_(matlabarray([[ 0, 1, 0, 2, 1, 0], [0, 0, 1, 0, 0.01, 2 ]]), 6)
        correctres = matlabarray([
    [1.0000,    1.0000,    1.0000,    1.0000,    1.0000,    1.0000],
     [    0,    1.0000,         0,    2.0000,    1.0000,         0],
      [   0,         0,    1.0000,         0,    0.0100,    2.0000],
       [  0,    0.5000,         0,    2.0000,    0.5000,         0],
        [ 0,         0,    0.5000,         0,    0.00005,    2.0000],
         [0,         0,         0,         0,    0.0100,         0]])
            
        self.assertTrue(compare_matlabarray(res, correctres, self.abs_tol, self.rel_tol))
    #evaluation of the  6 following points on the polynomial basis [1, x1, x2, 0.5*x1^2, 0.5*x2^2, x1*x2]
    def test_bcdfo_evalZ_2(self):
        res = bcdfo_evalZ_(matlabarray([[ 1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6 ]]), 6)
        correctres = np.array([
    [1.0000,    1.0000,    1.0000,    1.0000,    1.0000,    1.0000],
     [    1,         2,         3,         4,         5,         6],
      [   1,         2,         3,         4,         5,         6],
       [  0.5,         2,         4.5,         8,       12.5,         18],
        [ 0.5,         2,         4.5,         8,       12.5,         18],
         [1,         4,         9,         16,       25,         36]])
        self.assertTrue(compare_matlabarray(res, correctres, self.abs_tol, self.rel_tol))
        
    #evaluation of the  6 following points on the polynomial basis [1, x1, x2, 0.5*x1^2, 0.5*x2^2]
    def test_bcdfo_evalZ_3(self):
        res = bcdfo_evalZ_(matlabarray([[ 0, 1, 0, 2, 1, 0], [0, 0, 1, 0, 0.01, 2 ]]), 5)
        correctres = matlabarray([
    [1.0000,    1.0000,    1.0000,    1.0000,    1.0000,    1.0000],
     [    0,    1.0000,         0,    2.0000,    1.0000,         0],
      [   0,         0,    1.0000,         0,    0.0100,    2.0000],
       [  0,    0.5000,         0,    2.0000,    0.5000,         0],
        [ 0,         0,    0.5000,         0,    0.00005,    2.0000]]);
            
        self.assertTrue(compare_matlabarray(res, correctres, self.abs_tol, self.rel_tol))
        
   #evaluation of the  3 following points on the polynomial basis [1, x1, x2, 0.5*x1^2, 0.5*x2^2]
    def test_bcdfo_evalZ_4(self):
        res = bcdfo_evalZ_(matlabarray([[ 0, 1, 0], [0, 0, 1]]), 5)
        correctres = matlabarray([
    [1.0000,    1.0000,    1.0000],
     [    0,    1.0000,         0],
      [   0,         0,    1.0000],
       [  0,    0.5000,         0],
        [ 0,         0,    0.5000]]);
            
        self.assertTrue(compare_matlabarray(res, correctres, self.abs_tol, self.rel_tol))
        

if __name__ == '__main__':
    unittest.main()