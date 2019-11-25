# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 11:47:43 2014
INPUT VALUES

simul = 

    @sqpdfo_evalfgh


x =

  -0.003054026564966
                   0
   0.315521468125839


null_step =

     0


constrained_pbl =

     2


lm =

                   0
  -0.999999999499998
                   0
  -0.000000000375002
   0.000000000083334


M =

     1     0     0
     0     1     0
     0     0     1


n =

     3


me =

     2


mi =

     0


s =

  -0.503054026564966
  -1.000000000000000
  -0.184478531874161


gx =

  -0.009783923878659
   1.000000000000000
                   0


gci =

     []


gce =

   1.000000000000000   0.999999999999999   1.000000000000000
   1.000000000000000   1.999999999999999   3.000000000000000


info = 

         g: [3x1 double]
        ai: []
        ae: [2x3 double]
        hl: []
     niter: 1
      flag: 0
    nsimul: [0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        ci: []
        ce: [2x1 double]
         f: 0.099563123926544
      glag: [3x1 double]
     glagn: 0.666666666736111
     feasn: 3
     compl: 0



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

   0.312467441560872
  -0.056489622187450

K>> info.glag

ans =

  -0.333333333013890
   0.666666666736111
  -0.333333333513889


options = 

           algo_method: 101
    algo_globalization: 112
           hess_approx: 131
          bfgs_restart: 0
          algo_descent: 120
                   tol: [1.000000000000000e-05 1.000000000000000e-05 1.000000000000000e-05]
                 dxmin: 1.000000000000000e-06
                 miter: 500
                msimul: 500
               verbose: 2
                  fout: 1
                   inf: Inf
                   df1: 0


values = 

                       success: 0
              fail_on_argument: 1
               fail_on_problem: 2
                 fail_on_simul: 3
                 stop_on_simul: 4
              stop_on_max_iter: 5
             stop_on_max_simul: 6
                 stop_on_dxmin: 7
          fail_on_non_decrease: 8
            fail_on_ascent_dir: 9
           fail_on_max_ls_iter: 10
              fail_on_ill_cond: 11
    stop_on_small_trust_region: 15
             fail_on_null_step: 20
         fail_on_infeasible_QP: 21
          fail_on_unbounded_QP: 22
                  fail_strange: 99
                    nsimultype: 16
                max_null_steps: 1
                        newton: 100
                  quasi_newton: 101
            cheap_quasi_newton: 102
                 unit_stepsize: 110
                    linesearch: 111
                 trust_regions: 112
                        powell: 120
                         wolfe: 121
                          bfgs: 130
                         model: 131
                         dline: '--------------------------------------------------------------------------------------'
                         eline: '======================================================================================'
                         sline: '**************************************************************************************'


fcmodel =

   0.099563123926544  -0.013550429460611   1.384968815034240                   0   6.144976749236061
   0.312467441560872   1.384968815034241   1.384968815034240   1.384968815034241   0.000000000000006
  -0.056489622187450   1.384968815034241   2.769937630068480   4.154906445102721   0.000000000000014


Y =

  -0.003054026564966  -0.500000000000000   0.500000000000000   0.500000000000000   0.500000000000000
                   0   1.000000000000000                   0   1.000000000000000   1.000000000000000
   0.315521468125839   0.500000000000000   0.500000000000000  -0.500000000000000   0.500000000000000


fY =

   0.099563123926544   1.500000000000000   0.500000000000000   1.500000000000000   1.500000000000000


ciY =

     []


ceY =

   0.312467441560872   1.000000000000000   1.000000000000000   1.000000000000000   2.000000000000000
  -0.056489622187450   2.000000000000000   1.000000000000000                   0   3.000000000000000


sigma =

     1


scale =

   1.000000000000000
   0.722037918214987
   0.722037918214987
   0.722037918214987
   0.521338755340232


shift_Y =

     1


QZ =

   1.000000000000000                   0                   0                   0                   0
                   0  -0.437717028004984  -0.826365739060244   0.329535025867438  -0.130115853870500
                   0   0.880814115424577  -0.315023389192614   0.329471408091712  -0.127966204837389
                   0   0.162491294867564  -0.418566952172131  -0.884320607054428  -0.127966204837389
                   0   0.078529462977709  -0.206595343893201  -0.028849989870035   0.974843149132506


RZ =

   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000
                   0   0.819739267991799  -0.132165180682231   0.386491133279071   0.503816009553305
                   0                   0  -0.369537503265736  -0.294775124607443  -0.596996335387400
                   0                   0                   0   0.876403859762092   0.237890849609900
                   0                   0                   0                   0  -0.092396452142661


whichmodel =

     0


ind_Y =

     5     2     3     4     1


i_xbest =

     5


m =

     5
                    
OUPUT VALUES:

K>> M,pc,info

M =

   2.620301230585867                   0                   0
                   0                   0                   0
                   0                   0                   0


pc =

     1


info = 

         g: [3x1 double]
        ai: []
        ae: [2x3 double]
        hl: []
     niter: 1
      flag: 0
    nsimul: [0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        ci: []
        ce: [2x1 double]
         f: 0.099563123926544
      glag: [3x1 double]
     glagn: 0.666666666736111
     feasn: 3
     compl: 0

K>> info.g

ans =

  -0.009783923878659
   1.000000000000000
                   0

K>> info.ae

ans =

   1.000000000000000   0.999999999999999   1.000000000000000
   1.000000000000000   1.999999999999999   3.000000000000000

K>> info.ce

ans =

   0.312467441560872
  -0.056489622187450

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
from sqpdfo_computeHessian import *
from sqpdfo_evalfgh import *
#import numpy as np
import helper
from numpy import array

class dummyInfo():
    def __init__(self):
        self.g  = array([[ 0, 1, 0]]).T
        self.ai = array( [])
        self.ae = array([[1,     1,     1],[1,     2,     3]])
        self.hl = array( [])
        self.niter = 1
        self.flag = 0
        self.nsimul = array( [0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.ci = array( [])
        self.ce = array([[0.312467441560872, -0.056489622187450]]).T
        self.f = 0.099563123926544
        self.glag = array([[-0.333333333013890, 0.666666666736111,-0.333333333513889]]).T
        self.glagn =  0.666666666736111
        self.feasn = 3
        self.compl = 0
        
class Test_sqpdfo_computeHessian(unittest.TestCase):
    """
      Reminder :
      This class is a test for sqpdfo_computeHessian which computes an approximation of the Hessian of the Lagrangian
    """ 
    def setUp(self):
        self.simul = sqpdfo_evalfgh_
        self.x = array([[-0.003054026564966, 0, 0.315521468125839]]).T
        self.null_step = 0
        self.constrained_pbl = 2
        self.lm = array([[0, -0.999999999499998, 0,-0.000000000375002, 0.000000000083334]]).T
        self.M = array([[1,     0,     0],[0,     1,     0],[0,     0,     1]])
        self.n = 3
        self.me = 2
        self.mi =0
        self.s =array([ [ -0.503054026564966,-1.000000000000000,-0.184478531874161]]).T
        self.gx = array([[-0.009783923878659, 1.000000000000000, 0]]).T
        self.gci =array( [])
        self.gce = array([[1.000000000000000,   0.999999999999999,   1.000000000000000], [1.000000000000000 ,  1.999999999999999,   3.000000000000000]])

        self.info = dummyInfo()
        self.options = helper.dummyOptions()
        self.values = helper.dummyValues()
        
        self.fcmodel = array([[0.099563123926544,  -0.013550429460611,   1.384968815034240,  0,   6.144976749236061],
                        [0.312467441560872,   1.384968815034241,   1.384968815034240,   1.384968815034241,   0.000000000000006],
                          [-0.056489622187450,   1.384968815034241,   2.769937630068480 ,  4.154906445102721,   0.000000000000014]])


        self.Y = array([
            [  -0.003054026564966,  -0.500000000000000,   0.500000000000000,   0.500000000000000,   0.500000000000000],
            [               0,   1.000000000000000,                   0,   1.000000000000000,   1.000000000000000],
           [0.315521468125839,   0.500000000000000,   0.500000000000000,  -0.500000000000000,   0.500000000000000]])
        self.fY = array([   0.099563123926544,   1.500000000000000,   0.500000000000000,   1.500000000000000,   1.500000000000000])
        self.ciY = array(  [])
        self.ceY = array([ 
        [  0.312467441560872,   1.000000000000000 ,  1.000000000000000,   1.000000000000000 ,  2.000000000000000],
        [  -0.056489622187450,   2.000000000000000,   1.000000000000000,                   0,   3.000000000000000]])
        self.sigma = 1
        self.scale = array([[   1.000000000000000, 0.722037918214987, 0.722037918214987, 0.722037918214987, 0.521338755340232]]).T
        self.shift_Y = 1
        self.QZ = array([
               [1.000000000000000,                   0,                   0,                   0,                   0],
                [               0,  -0.437717028004984,  -0.826365739060244,   0.329535025867438,  -0.130115853870500],
                [               0,   0.880814115424577,  -0.315023389192614,   0.329471408091712,  -0.127966204837389],
                [              0 ,  0.162491294867564 , -0.418566952172131 , -0.884320607054428 , -0.127966204837389],
                [               0,   0.078529462977709,  -0.206595343893201,  -0.028849989870035,   0.974843149132506]])
        self.RZ = array([
               [1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000],
                   [               0,   0.819739267991799,  -0.132165180682231,   0.386491133279071,   0.503816009553305],
                [               0,                   0,  -0.369537503265736,  -0.294775124607443,  -0.596996335387400],
                [               0,                   0,                   0,   0.876403859762092,   0.237890849609900],
                [               0,                   0,                   0,                   0,  -0.092396452142661]])
        self.whichmodel = 0
        self.ind_Y = array([  4,     1,     2,     3,     0])
        self.i_xbest = 4
        self.m = 4
        
        self.abs_tol=1e-14
        self.rel_tol=1e-14

    def test_sqpdfo_computeHessian(self):
        """ 
         Test with some values, results compared with matlab 
        """
        M,pc,info = sqpdfo_computeHessian_(self.simul,self.x,self.null_step,self.constrained_pbl,
            self.lm,self.M,self.n,self.me,self.mi,self.s,self.gx,self.gci,self.gce,self.info,
            self.options,self.values,self.fcmodel,self.Y,self.fY,self.ciY,self.ceY,
            self.sigma,self.scale,self.shift_Y,self.QZ,self.RZ,self.whichmodel,
            self.ind_Y,self.i_xbest,self.m)
            
        #print("M:\n", M, "\npc:\n", pc, "\ninfo.g:\n", info.g)

        correctM = array([
            [          3.2036145300413947,                   0,                   0],
            [                  0,                   0,                   0],
            [                 0,                   0,                   0]])
        correctpc = 1
        correctinfog = array([
            [-0.009783923878659],
            [ 1.000000000000000],
            [0]])
        
        self.assertEqual(pc, correctpc)
        self.assertTrue(compare_array(M, correctM, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(correctinfog, info.g, self.abs_tol, self.rel_tol))
        


if __name__ == '__main__':
    unittest.main()