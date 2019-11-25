# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 14:50:56 2014
INPUT VALUES
 x,lm,lb,ub,

x =

   0.500000000000000
   1.000000000000000
   0.500000000000000


lm =

                   0
                   0
                   0
  -0.333333332763891
  -0.000000000249999


lb =

  -0.500000000000000
                   0
                -Inf


ub =

   Inf
   Inf
   Inf
            
    info

info = 

         g: [3x1 double]
        ai: []
        ae: [2x3 double]
        hl: []
     niter: 0
      flag: 0
    nsimul: [0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        ci: []
        ce: [2x1 double]
         f: 1.500000000000000
                                    
                                    
info.g

ans =

     0
     1
     0

K>> info.ae

ans =

     1     1     1
     1     2     3

K>> info.ce

ans =

     2
     3
                                    
            
OUTPUT VALUES

feas =

     0
     0
     0
     2
     3


compl =

     0
     0
     0


info = 

         g: [3x1 double]
        ai: []
        ae: [2x3 double]
        hl: []
     niter: 0
      flag: 0
    nsimul: [0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        ci: []
        ce: [2x1 double]
         f: 1.500000000000000
      glag: [3x1 double]
                        
                        
                        
                        
K>> info.g

ans =

     0
     1
     0

K>> info.ae

ans =

     1     1     1
     1     2     3

K>> info.ce

ans =

     2
     3

K>> info.glag

ans =

  -0.333333333013890
   0.666666666736111
  -0.333333333513889
@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from sqpdfo_optimality import *
#import numpy as np
import helper
from numpy import array

class dummyInfo():
    def __init__(self):
        self.glag = None
        
        
        self.g =  array([[0,   1,  0]]).T

        self.ai = array( [])
        self.ae = array([[1,    1,    1], [ 1,    2,    3]])
        self.hl = array( [])
        self.niter =  0
        self.flag = 0
        self.nsimul = array([0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.ci = array([])
        self.ce = array( [[     2,  3]]).T
        self.f = 1.500000000000000        

class Test_sqpdfo_optimality(unittest.TestCase):
    """
    #Reminder :
    #This class is a test for sqpdfo_optimality
    """ 
    def setUp(self):
        self.options = helper.dummyOptions()
        self.info = dummyInfo()
        
        self.x =array([[  0.500000000000000,  1.000000000000000,   0.500000000000000]]).T
        self.lm = array(   [ [0,0,0,  -0.333333332763891,  -0.000000000249999]]).T
        #print "lm", self.lm
        self.lb = array([[-0.500000000000000, 0, -np.Inf]]).T
        self.ub = array([[ np.Inf, np.Inf, np.Inf]]).T
        #self.values = helper.dummyValues()

    def test_sqpdfo_optimality(self):
        """
         Test with some values, results compared with matlab
         """
        feas,compl,info = sqpdfo_optimality_(self.x,self.lm,self.lb,self.ub,self.info,self.options)
        
        correctfeas = array([ 0, 0, 0, 2, 3]).T
        correctcompl = array([0, 0, 0]).T
        correctglag = array([-0.333333333013890, 0.666666666736111, -0.333333333513889]).T
        
        self.assertTrue(compare_array(feas, correctfeas, 1e-15,1e-15))
        self.assertTrue(compare_array(compl, correctcompl,1e-15,1e-15))
        self.assertTrue(compare_array(info.glag, correctglag, 1e-14, 1e-14))

if __name__ == '__main__':
    unittest.main()
    