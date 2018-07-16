# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 16:28:06 2014

@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from evalfgh import *
from ecdfo_func import *
from ecdfo_global_variables import set_prob
from numpy import array, double


class Test_evalfgh(unittest.TestCase):
    """
      Reminder :
      This class is a test for evalfgh
    """ 
    def setUp(self):
        self.key = 2
        self.xy = array([[-2, 2, 2, 1, 1]]).T
        set_simul_not_initialized(0)
        set_prob(5)

    def test_evalfgh(self):
        """
        Test with some values, results compared with matlab
        """
        msg,out2,out3,out4 = evalfgh_(self.key,self.xy,lm=None)

        correctmsg = 0
        correctout2 = 3.354626279025119e-04
        correctout4 =    array([4, -1, 1])
        
        self.assertEqual(correctmsg, msg)
        self.assertAlmostEqual(correctout2, double(out2), places=15)
        self.assertTrue(isempty_(out3))
        self.assertTrue(compare_array(correctout4, out4, 1e-15, 1e-15))
        

if __name__ == '__main__':
    unittest.main()
    