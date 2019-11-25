# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 11:35:42 2014
TEST:
%  X = [ 0 1 0 ; 0 0 1 ];
%  fX = [ 1 100 101 ];
%  y = [ 2; 4 ];
%
%  To put y on 4-th position in the set X call:
%  [X, fX, neval, xstatus, sstatus, dstatus] = sqpdfo_augmX_evalf(@banana, ...
%     y, 4, X, fX, 0, [0;0], [], [1 2 3], 1e25, 3, 100, 1, ...
%     0, [1 1 1], 1, [1 1 1], [ 0 0 0], 0, [ 1, 1])
%  where
%     function fx = banana( x )
%     fx  = 100 * ( x(2) - x(1)^2 ) ^2 + (1-x(1))^2;
%
%  which gives
%  X =
%     0     1     0     2
%     0     0     1     4
%  fX =
%     1     100     101   1
%  neval =
%     4
%  xstatus =
%     1     1     1     1
%  sstatus =
%     1     1     1     1
%  dstatus =
%     0     0     0     0
@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from sqpdfo_augmX_evalf import *
#import numpy as np
import helper
from numpy import array
from copy import copy
from runtime import compare_array

class dummyInfo():
    def __init__(self):
        self.nsimul = array([0,0,0])
        self.ci = array([])
        self.ce = array([])

class Test_sqpdfo_augmX_evalf(unittest.TestCase):
    """
          Reminder :
          This class is a test for sqpdfo_augmX_evalf which adds a new points y to the set of all points X.
    """
    def setUp(self):
        self.X = array([[ 0, 1, 0],[0, 0, 1] ])
        self.fX = array([ 1, 100, 101 ])
        self.y = array([[ 2, 4 ]]).T
        
        self.ciX = copy(self.X)
        self.ceX = copy(self.X)
        
        self.info = dummyInfo()
        self.options = helper.dummyOptions()
        self.values = helper.dummyValues()
        self.abs_tol=1e-15
        self.rel_tol=1e-15
    #def banana(self, x, y):
#        #return outdic,fvalue,info.ci,info.c 
#        print "x", x
#        print "y", y
#        return None,100 * ( y[2] - y[1]^2 ) ^2 + (1-y[1])^2,array([]),array([])
        
    def banana(self, somenumber, x):

        #return outdic,fvalue,info.ci,info.c 
        #print "x", x
        #print "y", y
        #print "banana x", x
        #print "summand 1:", 100 * ( x[2] - x[1]**2 )**2
        #print "summand 2:", (1-x[1])**2
        #print "[[0]] + [[-3]]", array([0]) + array([-3])
        #print "[[-3]]^2 = ", array([-3])**2
        fx = 100 * ( x[1] - x[0]**2 )**2 + (1-x[0])**2
        #print "banana fx", fx
    #outdic = 0
     #fvalue = fx
     #infoci = fx
     #infoce = fx
        return 0,fx,array([[1.0],[0.0]]), array([[0.0],[1.0]])#fx,fx


    #@unittest.expectedFailure
    def test_sqpdfo_augmX_evalf(self):
        """
          This is the matlab test.
        """
        #(f=None,y=None,m=None,X=None,fX=None,ciX=None,ceX=None,nfix=None,xfix=None,indfix=None,
        #indfree=None,fxmax=None,neval=None,xstatus=None,xstatus_val=None,sstatus=None,dstatus=None,
        #scaleX=None,scalefacX=None,info=None,options=None,values=None)
        #returns X, fX, neval, xstatus, sstatus, dstatus
    
        X,fX,ciX,ceX,neval,xstatus,sstatus,dstatus,info,outdic = sqpdfo_augmX_evalf_(
        self.banana,self.y, 3,self.X,
        self.fX,self.ciX, self.ceX, 0, array([[0],[0]]), array([]), array([[0, 1, 2]]), 
        1e25, 3, array([1, 1, 1]), 1,
        array([1, 1, 1]), array([[ 0, 0, 0]]), 0, 
        array([[ 1, 1]]), self.info,self.options, self.values)
        
        #sstatus=1, dstatus=array([1, 1, 1]), scaleX=array([ 0, 0, 0]), scalefacX=0, 
        #info=array([ 1, 1]))
        correctX = array([[0,     1,     0,     2],[0,     0,     1,     4]])  
        correctfX = array([1,     100,     101,   1])
        correctneval = 4
        correctxstatus =array([ 1,     1,     1,     1])
        correctsstatus =array([ 1,     1,     1,     1])
        correctdstatus = array([0,     0,     0,     0])
        
        self.assertTrue(compare_array(X, correctX,self.abs_tol,self.rel_tol))
        self.assertTrue(compare_array(sstatus, correctsstatus,self.abs_tol,self.rel_tol))
        self.assertTrue(compare_array(dstatus, correctdstatus,self.abs_tol,self.rel_tol))
        self.assertTrue(compare_array(xstatus, correctxstatus,self.abs_tol,self.rel_tol))
        self.assertTrue(compare_array(fX, correctfX,self.abs_tol,self.rel_tol))
        self.assertEqual(neval, correctneval)

if __name__ == '__main__':
    unittest.main()