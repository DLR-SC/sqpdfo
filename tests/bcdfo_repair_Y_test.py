# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 15:16:08 2014

@author: jaco_da
"""

import sys
sys.path.append("../")

import unittest
from bcdfo_repair_Y import bcdfo_repair_Y_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from runtime import matlabarray, compare_matlabarray
#import numpy as np
#import helper

class Test_bcdfo_repair_Y(unittest.TestCase):
    """
      Reminder :
      This class is a test for bcdfo_repair_Y  which epairs the interpolation set Y by first replacing the interpolation points
      at a distance larger that farfact*Delta from Y(:,1), and then replaces 
      interpolation points with the best possible point obtained by taking the best
      optimal replacement, until this improvement no longer causes a relative simplex
      volume increase by at least threshold.  For conventions on how polynomials
      are represented, see the documentation of evalZ.

    """ 
    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-10
        self.rel_tol=1e-10
        pass

    def test_bcdfo_repair_Y(self):
#        TEST:
        """
        This is the test written in the Matlab Code. Results are the same except for a few signs due to a non-unique QR decomposition.
        """
        Y = matlabarray([[ 0, 1, 0, 2, 1, 0],[0, 0, 1, 0, 0.01, 2 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(  Y, 0, 0, 1, 1, 1e15 )
        QZplus, RZplus, Yplus, replaced, maximprove, Y_radius, xbase, scale = bcdfo_repair_Y_( QZ, RZ, Y, 0.7, 10, 1.0e-10, 1.1, 0.001, xbase, 1, 0, 0, matlabarray([-10,-10]), matlabarray([10,10]), matlabarray([1,2]), 1, scale, 0, 1, 1e15 )
        #print QZplus, RZplus, Yplus, replaced, maximprove, Y_radius, xbase, scale
        
        correctQZplus = matlabarray([

   [1 ,	0,	0,	0,	0,	0],
   [0,	-0.659566419393763,	0.583075869829704,	-0.259014354637790,	0.383735140707789,	0.103216153328673],
   [0	,-0.640135734527690,	-0.703005089235313,	-0.116047869763477,	-0.172250116473354,	0.229941025457003],
   [0,	0.165662796629888,	0.0580981831826901,	-0.900707614268715,	-0.383735140707789,	-0.103216153328673],
   [0,	0.156045791159224,	0.217778643032798,	0.0539460623491271,	-0.274281054964081,	0.921998860598591],
   [0,	0.321564812538689	,-0.339121568587501,	-0.324438086497736,	0.774980311901311	, 0.275205518205477]])


        correctRZplus = matlabarray([

  [1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000],
  [ 0,	0.761619249046136,	0.120705803297469,	-0.987807245527749,	0.141209675627013,	-0.968179886736932],
  [      0,        0,  0.747318893145360,	1.28234810602479,	-0.504495501740099,	-0.970452892405031],
  [      0,        0,        0,   -2.31944393781301	,0.0241326579145779,	-0.124203614828700],
  [      0,        0,        0,        0,  -0.545860227607381,	-0.893062342874869],
  [      0,        0,        0,        0,        0,  2.30387977211119 ]])

        correctYplus = matlabarray([
[0,	-0.502338481034726,	0.356130119179943	,  2,	-0.603012769702487,	0],
[0,	-0.487539697418576,	-0.602637083218470,	0,	0.355493490030523,	2]])


        correctreplaced = matlabarray([2,    3,    5]) 

        correctmaximprove =  1.000153015137598

        correctY_radius =  2

        self.assertTrue(compare_matlabarray(correctRZplus, RZplus, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctQZplus, QZplus, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctYplus, Yplus, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctreplaced, replaced, self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(correctmaximprove, maximprove, places=13)
        self.assertEqual(correctY_radius, Y_radius)
        
        #The scaled version (i.e. with shift in interpolation points)
        Y = matlabarray([[ 0, 1, 0, 2, 1, 0],[0, 0, 1, 0, 0.01, 2 ]])
        QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(  Y, 0, 1, 1, 1, 1e15 )
        QZplus, RZplus, Yplus, replaced, maximprove, Y_radius, xbase, scale = bcdfo_repair_Y_( QZ, RZ, Y, 0.7, 10, 1.0e-10, 1.1, 0.001, xbase, 1, 0, 0, matlabarray([-10,-10]), matlabarray([10,10]), matlabarray([1,2]), 1, scale, 1, 1, 1e15 )
        #print QZplus, RZplus, Yplus, replaced, maximprove, Y_radius, xbase, scale
        
        
        correctQZplus = matlabarray([

  [1.0000,         0,         0,         0,         0 ,        0],
  [0,	-0.701665921132481,	0.663289830649876,	-0.149898529361902,	0.206022638530179,	0.0528832580876775
],
  [0,	-0.680994993999292,	-0.714381653387058,	-0.0589306015912921,	-0.0924789168850424,	0.117811313462405
],
  [  0,	0.0881184482538805,	0.0261672408486609,	-0.900323579771228,	-0.412045277060358,	-0.105766516175355
],
  [  0,	0.0830030233294948,	0.111907596144591,	0.0351764352940761,	-0.294515944525440,	0.944780485014161
],
  [0,	0.171044995438928,	-0.191042977604362,	-0.402787394538547,	0.832154655348015,	0.282005557794224]])

        correctRZplus = matlabarray([

 [ 1.0000,  1.0000,   1.0000,   1.0000,   1.0000,   1.0000],
 [ 0,	0.357961293192034	,0.0762419683750912,	-0.657606697005541,	0.0866619532593967,	-0.639493482334544
],
 [       0,         0,   0.349110509250466,	0.676373451074206	,-0.313769619685537,	-0.658427855314762
],
 [       0,         0,         0,    -0.600060319247516,	0.0159399616505447,	-0.0413423839442541
],
 [       0,         0,         0,         0,    -0.146532784415969,	-0.239736889147762
],
 [       0,         0,         0,         0,         0,   0.590201555969485]])

        correctYplus = matlabarray([
[0,	-0.502338481034726,	0.356130119179943	,  2,	-0.603012769702487,	0],
[0,	-0.487539697418576,	-0.602637083218470,	0,	0.355493490030523,	2]])


        correctreplaced = matlabarray([2,    3,    5]) # changed to python indices

        correctmaximprove =  1.000153015137598

        correctY_radius =  2

        self.assertTrue(compare_matlabarray(correctRZplus, RZplus, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctQZplus, QZplus, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctYplus, Yplus, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctreplaced, replaced, self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(correctmaximprove, maximprove, places=13)
        self.assertEqual(correctY_radius, Y_radius)
  

if __name__ == '__main__':
    unittest.main()