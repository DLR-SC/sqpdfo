# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 15:13:04 2014

@author: jaco_da
"""

import sys
sys.path.append("../")
import unittest
from bcdfo_projgrad import bcdfo_projgrad_
from runtime import matlabarray, compare_matlabarray
#import numpy as np
#import helper

class Test_bcdfo_projgrad(unittest.TestCase):
    """
      Reminder :
      This class is a test for bcdfo_projgrad which computes the projected gradient and its infinity norm.
    """ 
    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-14;
        self.rel_tol=1e-14;
        pass

    def test_bcdfo_projgrad(self):
        """
            Small tests using the logic answer of the projected gradient to see if the function works correctly
        """
        #print "No test specified..."
#          n     : dimension
#  x     : current iterate
#  g     : gradient at x
#  bl    : lower bounds
#  bu    : upper bounds
        
        #TEST 1 WITH UNDISTURBING BOUNDS 
        n = 2
        x = matlabarray([0.5,0.5]).T
        g = matlabarray([-1,-1]).T
        bl = matlabarray([0,0]).T
        bu = matlabarray([2,2]).T
        
        gnorm, gn = bcdfo_projgrad_(n,x,g,bl,bu)
        
        #print "gnorm", gnorm
        #print "gn", gn
        
        correctgn = matlabarray([-1,-1]).T
        correctgnorm = 1
        
        self.assertTrue(compare_matlabarray(correctgn, gn, self.abs_tol, self.rel_tol))
        self.assertEqual(gnorm, correctgnorm)

        #TEST 2 WITH A DISTURBING UPPER BOUND
        n = 2
        x = matlabarray([0.5,0.5]).T
        g = matlabarray([-1,-1]).T
        bl = matlabarray([0,0]).T
        bu = matlabarray([1,2]).T
        
        gnorm, gn = bcdfo_projgrad_(n,x,g,bl,bu)
        
        #print "gnorm", gnorm
        #print "gn", gn
        
        correctgn = matlabarray([-0.5,-1]).T
        correctgnorm = 1
        
        self.assertTrue(compare_matlabarray(correctgn, gn, self.abs_tol, self.rel_tol))
        self.assertEqual(gnorm, correctgnorm)
        
        
        #1 TEST 3 WITH A DISTURBING LOWER BOUND
        n= 2
        x = matlabarray([0.5,0.5]).T
        g = matlabarray([-1,1]).T
        bl = matlabarray([0,0]).T
        bu = matlabarray([2,2]).T
        
        gnorm, gn = bcdfo_projgrad_(n,x,g,bl,bu)
        
        #print "gnorm", gnorm
        #print "gn", gn
        
        correctgn = matlabarray([-1,0.5]).T
        correctgnorm = 1
        
        self.assertTrue(compare_matlabarray(correctgn, gn, self.abs_tol, self.rel_tol))
        self.assertEqual(gnorm, correctgnorm)

if __name__ == '__main__':
    unittest.main()