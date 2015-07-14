# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 15:23:24 2014

@author: jaco_da
"""
import sys
sys.path.append("../")
sys.path.append("tests/")

import unittest
from bcdfo_replace_in_Y import *#bcdfo_replace_in_Y_
from bcdfo_build_QR_of_Y import *#bcdfo_build_QR_of_Y_
from runtime import matlabarray, compare_matlabarray
import numpy as np
#import helper

class Test_bcdfo_replace_in_Y(unittest.TestCase):
    """
      Reminder :
      This class is a test for bcdfo_replace_in_Y which updates the interpolation set for a transformation of Y in Yplus where the
      vector ynew replaces Y(:,j).  Also update the factorization of Z(Y)
      accordingly.
    """ 
    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-10;
        self.rel_tol=1e-10;
        pass

    #@unittest.expectedFailure
    def test_bcdfo_replace_in_Y(self):
        """
        This is the test written in the Matlab Code. Results are the same except for maybe some signs which may be different due to a slightly different QR factorization.
        """
        Y = matlabarray([[ 1.0, 1.0, 0.0, 0.0, 0.0, -1.0], [1.0, 0.0, -1.0, 1.0, 0.0, 0.0, ]])
        whichmodel = 0
        QZ,RZ,xbase,scale = bcdfo_build_QR_of_Y_( Y , 0, 0, 1, 1, 1e15 )
        #print "QZ\n", QZ
        #print "RZ\n", RZ

        ynew = 0.25*Y[:,3].T
        #print "ynew", ynew
        QZplus, RZplus, Yplus, xbase, scale = bcdfo_replace_in_Y_( QZ, RZ, ynew, Y, 3, xbase, whichmodel, scale, 0, 1, 1, 1e15 )
        
        #print "QZplus:\n", QZplus, "\nRZplus:\n",  RZplus, "\nYplus:\n", Yplus
        
        correctYplus = matlabarray([
    [1.0000,    1.0000,         0,         0,         0,   -1.0000],
    [1.0000,         0,   -0.2500,    1.0000,         0,         0]])
                  
        correctQZplus = matlabarray([
   [-4.714045207910318e-01,   -4.714045207910318e-01,                    0.720457635944399,    0.152655619705795	, 0.114859096884849	,3.42976275268113e-17],
   [-4.714045207910317e-01,   -4.714045207910317e-01,   -0.576366108755519,   -0.122124495764636,	-0.0918872775078794	,-0.447213595499958],
   [-4.714045207910317e-01,   4.714045207910317e-01,   -0.189120129435405,    0.633289525385556,	0.344577290654547,	1.09881792804699e-17],
   [-2.357022603955159e-01,   -2.357022603955158e-01,   -0.288183054377760,   -0.0610622478823181,	-0.0459436387539394,	0.894427190999916],
   [-2.357022603955159e-01,    2.357022603955158e-01,    0.108068645391660,   0.181336372499005,	-0.918872775078793,	-2.55294270180895e-16],
   [-4.714045207910317e-01,    4.714045207910317e-01,    0.135085806739575,   -0.723957711635059,	0.114859096884849	,-1.00034123182046e-16]])

        correctRZplus=matlabarray([
[-2.12132034355964,	-1.06066017177982,	-0.360919086230634,	-1.06066017177982	,-0.471404520791032,	-0.117851130197758],
[0,	-1.06066017177982,	-0.581889955351430,	0.117851130197758,	-0.471404520791032,	-0.117851130197758],
[0,	0,	0.771114813471740,	0.585371829204824,	0.720457635944399	,1.15273221751104],
[0,	0,	0,	0.876613331340854,	0.152655619705795,	0.244248991529272],
[0,	0,	0,	0,	0.114859096884849,	0.183774555015759],
[0,	0,	0,	0,	0,	0.894427190999916]])
         
#        print QZplus
        self.assertTrue(compare_matlabarray(correctQZplus, QZplus,self.abs_tol, self.rel_tol)) 
        self.assertTrue(compare_matlabarray(correctRZplus, RZplus,self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctYplus, Yplus,self.abs_tol, self.rel_tol))

        #the same test with scale version :
        Y = matlabarray([[ 1.0, 1.0, 0.0, 0.0, 0.0, -1.0], [1.0, 0.0, -1.0, 1.0, 0.0, 0.0, ]])
        whichmodel = 0
        QZ,RZ,xbase,scale = bcdfo_build_QR_of_Y_( Y , 0, 1, 1, 1, 1e15 )

        ynew = 0.25*Y[:,3].T
        
        QZplus, RZplus, Yplus, xbase, scale = bcdfo_replace_in_Y_( QZ, RZ, ynew, Y, 3, xbase, whichmodel, scale, 1, 1, 1, 1e15 )
        
        correctYplus = matlabarray([
    [1.0000,    1.0000,         0,         0,         0,   -1.0000],
    [1.0000,         0,   -0.2500,    1.0000,         0,         0]])
    
        correctQZplus=matlabarray([
    [-1,	0,	0,	0,	0,	0],
[0,	0,	-0.855247610220823,	0.470034549360494,	3.88578058618805e-16,	-0.218217890235992],
[0,	0.975900072948533	, 0.0127268989616193,	0.0231571273585230,	-0.216612144429553,	-4.44089209850063e-16],
[0,	0,	0.191239179404800,	-0.105102920414355,	1.91166527052644e-15,	-0.975900072948533],
[0,	-0.218217890235992,	0.0569164224419050,	0.103561821874555	  ,-0.968718959392966,	-1.85962356624714e-15],
[0,	0,	0.478097948512001, 	0.869919303746259	,0.121089869924121,	6.03683769639929e-16]])
        
        correctRZplus=matlabarray([
[-1,	-1,	-1,	-1,	-1,	-1],
[0,	-0.458257569495584,	-0.579641270939355,	0,	-0.458257569495584,	-0.458257569495584],
[0,	0,	0.522905402079391	, 0.401602276750081,	0.497221866452481	, 1.03269156878592],
[0,	0,	0,	-0.220716132870144,	-0.0467322721208927,	-0.114485128324656],
[0,	0,	0,	0,	0.0242179739848241,	0.0484359479696486],
[0,	0,	0,	0,	0,	-0.195180014589707]])

        self.assertTrue(compare_matlabarray(correctQZplus, QZplus,self.abs_tol, self.rel_tol)) 
        self.assertTrue(compare_matlabarray(correctRZplus, RZplus,self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctYplus, Yplus,self.abs_tol, self.rel_tol))

if __name__ == '__main__':
    unittest.main()