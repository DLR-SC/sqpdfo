# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 16:28:06 2014
INPUT VALUES
key,xy,lm

key =

     2


xy =

    -2
     2
     2
     1
     1

lm = none (Not enough input arguments.)
OUTPUT VALUES
 outdic,fvalue,info.ci,info.ce

outdic =

     0


fvalue =

     3.354626279025119e-04


ans =

     []


ans =

     4
    -1
     1

@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from evalfgh import *
from ecdfo_func import *
from ecdfo_global_variables import set_prob
#import numpy as np
#import helper

class Test_evalfgh(unittest.TestCase):
    """
      Reminder :
      This class is a test for evalfgh
    """ 
    def setUp(self):
        self.key = 2
        self.xy = matlabarray([-2, 2, 2, 1, 1]).T
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        
        set_simul_not_initialized(0)
        set_prob(5)

    def test_evalfgh(self):
        """
        Test with some values, results compared with matlab
        """
        msg,out2,out3,out4 = evalfgh_(self.key,self.xy,lm=None)

        #print "msg,out2,out3,out4", msg,out2,out3,out4
        
        correctmsg = 0
        correctout2 = 3.354626279025119e-04
        correctout3 =  matlabarray([])
        correctout4 =    matlabarray([4, -1, 1]).T
        
        self.assertEqual(correctmsg, msg)
        self.assertAlmostEqual(correctout2, out2, places=15)
        self.assertEqual(correctout3, out3)
        
        self.assertTrue(compare_matlabarray(correctout4, out4, 1e-15, 1e-15))
        
        #self.assertEqual(correctout4, out4)


if __name__ == '__main__':
    unittest.main()
    