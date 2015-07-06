# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 13:35:53 2014

@author: jaco_da
"""
import sys
sys.path.append("../")
 
import unittest
from bcdfo_evalZ import *
import numpy as np
#import helper

class Test_bcdfo_evalZ(unittest.TestCase):

    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()

        pass

    #@unittest.expectedFailure
    def test_bcdfo_evalZ(self):
        res = bcdfo_evalZ_(matlabarray([[ 0, 1, 0, 2, 1, 0], [0, 0, 1, 0, 0.01, 2 ]]), 6)
        #print "res", res
        correctres = matlabarray([
    [1.0000,    1.0000,    1.0000,    1.0000,    1.0000,    1.0000],
     [    0,    1.0000,         0,    2.0000,    1.0000,         0],
      [   0,         0,    1.0000,         0,    0.0100,    2.0000],
       [  0,    0.5000,         0,    2.0000,    0.5000,         0],
        [ 0,         0,    0.5000,         0,    0.0001,    2.0000],
         [0,         0,         0,         0,    0.0100,         0]])
            
        #print "abs:\n\n", abs(res - correctres)
        #print "5e-5  -> 1e-1 ?"
        # Just a numerical Error, Tested also with the below test case and correct result was obtained
        self.assertTrue((abs(res - correctres) < 1e-4).all())
        
    def test_bcdfo_evalZ_Original(self):
        #ret =    bcdfo_evalZ( np.array([[ 0.0, 1.0, 0.0, 2.0, 1.0, 0.0], [0.0, 0.0, 1.0, 0.0, 0.01, 2.0 ]]), 6 )
        res = bcdfo_evalZ_(matlabarray([[ 1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6 ]]), 6)
        correctres = np.array([
    [1.0000,    1.0000,    1.0000,    1.0000,    1.0000,    1.0000],
     [    1,         2,         3,         4,         5,         6],
      [   1,         2,         3,         4,         5,         6],
       [  0.5,         2,         4.5,         8,       12.5,         18],
        [ 0.5,         2,         4.5,         8,       12.5,         18],
         [1,         4,         9,         16,       25,         36]])
        self.assertTrue((abs(res - correctres) < 1e-4).all())
        
        
        


#  should give
# ans =
#
#    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000
#         0    1.0000         0    2.0000    1.0000         0
#         0         0    1.0000         0    0.0100    2.0000
#         0    0.5000         0    2.0000    0.5000         0
#         0         0    0.5000         0    0.0001    2.0000
#         0         0         0         0    0.0100         0
#

if __name__ == '__main__':
    unittest.main()