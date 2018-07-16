# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 13:39:12 2014

@author: jaco_da
"""


import sys
sys.path.append("../")
import unittest
from bcdfo_computeP import *
from bcdfo_build_QR_of_Y import *
from bcdfo_evalP import *
from random import random
from runtime import compare_array
from numpy import array


class Test_bcdfo_computeP(unittest.TestCase):
    """
    Reminder :
    This class is a test for computeP which computes the coefficients of the polynomial model
    """ 
    def setUp(self):      
        self.abs_tol=1e-14;
        self.rel_tol=1e-14;
        pass
    
    def test_bcdfo_computeP_0(self):
        """
        Easy test which verify that the model computed is indeed an interpolant model, i.e. that it interpolates all the initial points used to create the model,
        at least when there are enough point to create a complete quadratic model (for instance, with points with 2 coordinates, we need 6 points to 
        have a fully determined quadratic model (with random points, poisedness is assumed to be verified))
        """
        for i in range(0,10):
            x1_coord=[(random()-0.5)*100 for p in range(0,6)]
            x2_coord=[(random()-0.5)*100 for p in range(0,6)]
            y_coord=[(random()-0.5)*100 for p in range(0,6)]
            Y = array([x1_coord,x2_coord])
            fY =array([y_coord])
            whichmodel = 0
        
            QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(Y, whichmodel, 1, 1, 1, 1e15 )
            P = bcdfo_computeP_( QZ, RZ, Y, fY, whichmodel, array([0, 0, 0, 0, 0,0]), 0,
                                1, 1, array([0, 0]), scale, 1, 1)
            for j in range(0,6):
                value=bcdfo_evalP_( P, array([[x1_coord[j]],[x2_coord[j]]]), xbase, scale, 1 )
                self.assertAlmostEqual(double(value), y_coord[j],places=8)
    
    def test_bcdfo_computeP_1(self):
        """
        More complicated test initally written in the matlab code
        """
        Y = array([[ 0, 1, 0, 2, 0],[0, 0, 1, 0, 2]])
        fY =array([[ 1, 2, 3, 4, 5 ]])
        whichmodel = 0
        
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(Y, whichmodel, 1, 1, 1, 1e15 )
        #QZ = array(QZ)
        #RZ = array(RZ)
        
        #xbase = array(xbase)
        #scale = array(scale)
        P = bcdfo_computeP_( QZ, RZ, Y, fY, whichmodel, array([0, 0, 0, 0, 0]), 0,
            1, 1, array([0, 0]), scale, 1, 1)
        self.assertTrue(compare_array(P, array([1, 1, 4, 4, 0]), self.abs_tol, self.rel_tol))

    def test_bcdfo_computeP_2(self):
        """
        See test_bcdfo_computeP_0 with whichmodel = 2.
        """
        for i in range(0, 10):
            x1_coord = [(random() - 0.5) * 100 for p in range(0, 6)]
            x2_coord = [(random() - 0.5) * 100 for p in range(0, 6)]
            y_coord = [(random() - 0.5) * 100 for p in range(0, 6)]
            Y = array([x1_coord, x2_coord])
            fY = array([y_coord])
            whichmodel = 2

            QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(Y, whichmodel, 1, 1, 1, 1e15)
            P = bcdfo_computeP_(QZ, RZ, Y, fY, whichmodel, array([0, 0, 0, 0, 0, 0]), 0,
                                1, 1, array([0, 0]), scale, 1, 1)
            for j in range(0, 6):
                value = bcdfo_evalP_(P, array([[x1_coord[j]], [x2_coord[j]]]), xbase, scale, 1)
                self.assertAlmostEqual(double(value), y_coord[j], places=8)


if __name__ == '__main__':
    unittest.main()