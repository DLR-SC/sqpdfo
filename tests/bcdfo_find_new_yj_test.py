# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 11:17:28 2014

@author: jaco_da
"""
import sys
sys.path.append("../")

import unittest
from bcdfo_find_new_yj import bcdfo_find_new_yj_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from runtime import compare_array
from numpy import array, double


class Test_bcdfo_find_new_yj(unittest.TestCase):
    """
      Reminder :
      This class is a test for bcdfo_find_new_yj which compute the best point to replace yj among the current set of interpolation point Y,
      and also yields the improvement of poisedness corresponding to |L_j(new_y)|.
    """ 
    def setUp(self):
        self.abs_tol=1e-15;
        self.rel_tol=1e-15;
        pass

    def test_bcdfo_find_new_yj(self):
        """
            This is the test written in the Matlab Code. By Looking carefully, the result is slightly different because the sign of ynew[1] is opposite
            to the sign computed in matlab. This, after investigation, comes from a difference in the computation of pvalue & mvalue on approximately the 
            fiftheen significant digit, which leads to mvalue < pvalue in python instead of mvalue > pvalue in matlab. It is interesting to keep this in
            memory, but this is not really a bug.
        """
        Y = array([[ 3.0, 1.0, 0, 2.0, 1.0, 0.0],[0.0, 0.0, 1.0, 0.0, 0.01, 2.0 ]])
        whichmodel = 0
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y , whichmodel, 1, 1,1, 1e15 )
        ynew, improvement,msgTR = bcdfo_find_new_yj_( QZ, RZ, Y, 4, 1.0, 0.001, xbase, 1, whichmodel, scale, 1)

        correctynew = array([[3.280776988023534, -0.959774286524791]]).T
        correctimprovement =  314.8805392927235
        self.assertAlmostEqual(correctimprovement, double(improvement), 11)
        self.assertTrue(compare_array(correctynew, ynew, self.abs_tol, self.rel_tol))
        
        #Same test as above but without the shifting in the interpolation points
        Y = array([[ 3.0, 1.0, 0, 2.0, 1.0, 0.0],[0.0, 0.0, 1.0, 0.0, 0.01, 2.0 ]])
        whichmodel = 0
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y , whichmodel, 0, 1,1, 1e15 )
        ynew, improvement, msgTR = bcdfo_find_new_yj_( QZ, RZ, Y, 4, 1.0, 0.001, xbase, 1, whichmodel, scale, 0 )

        correctynew = array([[3.280776988023534, -0.959774286524791]]).T
        correctimprovement =  314.8805392927235
        self.assertAlmostEqual(correctimprovement, double(improvement), places=11)
        self.assertTrue(compare_array(correctynew, ynew, self.abs_tol, self.rel_tol))


if __name__ == '__main__':
    unittest.main()