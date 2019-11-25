# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 14:16:23 2014
INPUT VALUES

prob=5

OUTPUT VALUES

 [x0,lx,ux,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info] = sqpdfo_init_prob(5)

x0 =

    -2
     2
     2
     1
     1


lx =

  -Inf
  -Inf
  -Inf
  -Inf
  -Inf


ux =

   Inf
   Inf
   Inf
   Inf
   Inf


dxmin =

   1.0000e-06


li =

     []


ui =

     []


dcimin =

   1.4901e-08


infb =

   1.0000e+20


n =

     5


nb =

     0


mi =

     0


me =

     3


info =

     0


@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from sqpdfo_init_prob import *
from numpy import array
#import numpy as np
#import helper

class Test_sqpdfo_init_prob(unittest.TestCase):
    """
        Reminder :
        This class is a test for sqpdfo_init_prob which initializes the problem
    """ 
    #def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()

    #    pass

    def test_sqpdfo_init_prob(self):
        """
            Test with some values, results compared with matlab
        """
        x0,lx,ux,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info = sqpdfo_init_prob_(5)
        
        correctx0 = array([[-2, 2, 2, 1, 1]]).T
        correctlx = array( [[ -np.Inf, -np.Inf, -np.Inf, -np.Inf, -np.Inf]]).T
        correctux = array([[ np.Inf, np.Inf, np.Inf, np.Inf, np.Inf]]).T
        correctdxmin =  1.0000e-06
        correctdcimin = 1.4901e-08
        correctinfb = 1.0000e+20
        correctn =  5
        correctnb =    0
        correctmi =     0
        correctme =   3
        correctinfo = 0
        
        self.assertEqual(correctdxmin, dxmin)
        self.assertTrue(isempty_( li))
        self.assertTrue(isempty_( ui))
        self.assertAlmostEqual(correctdcimin, dcimin, 4)
        self.assertEqual(correctinfb, infb)
        self.assertEqual(correctn, n)
        self.assertEqual(correctnb, nb)
        self.assertEqual(correctmi, mi)
        self.assertEqual(correctme, me)
        self.assertEqual(correctinfo,info)

        #print "thru onon vec"   
        self.assertTrue(compare_array(correctx0, x0, 1e-15, 1e-15))
        self.assertTrue(compare_array(correctux, ux, 1e-15, 1e-15))
        self.assertTrue(compare_array(correctlx, lx, 1e-15, 1e-15))     
          
if __name__ == '__main__':
    unittest.main()
    