# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 17:27:08 2014
-------------------INPUT VALUES--------------------------


func = 

    @evalfgh


n =

     3


nb =

     2


mi =

     0


me =

     2


lm =

                   0
                   0
                   0
  -0.333333332763891
  -0.000000000249999


nitold =

     0


nit =

     0


i_xbest =

     1


lb =

  -0.500000000000000
                   0
                -Inf


ub =

   Inf
   Inf
   Inf


m =

     4


X =

   0.500000000000000  -0.500000000000000   0.500000000000000   0.500000000000000
   1.000000000000000   1.000000000000000                   0   1.000000000000000
   0.500000000000000   0.500000000000000   0.500000000000000  -0.500000000000000


fX =

   1.500000000000000   1.500000000000000   0.500000000000000   1.500000000000000


ciX =

     []


ceX =

     2     1     1     1
     3     2     1     0


ind_Y =

     1     2     3     4


QZ =

     1     0     0     0
     0     1     0     0
     0     0     1     0
     0     0     0     1


RZ =

     1     1     1     1
     0    -1     0     0
     0     0    -1     0
     0     0     0    -1


delta =

     1


cur_degree =

     4


neval =

     0


maxeval =

   600


maxit =

   600


fcmodel =

   1.500000000000000                   0   1.000000000000000                   0
   2.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000
   3.000000000000000   1.000000000000000   2.000000000000000   3.000000000000000


gx =

     0
     1
     0


normgx =

     1


show_errg =

     0


pquad =

    10


pdiag =

     7


plin =

     4


stallfact =

     2.220446049250313e-15


eps_rho =

     1.000000000000000e-14


Deltamax =

      100000


rep_degree =

     4


epsilon =

     1.000000000000000e-05


verbose =

     1


eta1 =

     1.000000000000000e-04


eta2 =

   0.900000000000000


gamma1 =

   0.010000000000000


gamma2 =

   0.500000000000000


gamma3 =

     2


interpol_TR =

     1


factor_CV =

   100


Lambda_XN =

     1.000000000000000e-10


Lambda_CP =

   1.200000000000000


factor_FPU =

     1


factor_FPR =

    10


Lambda_FP =

     1.000000000000000e-10


criterion_S =

distance


criterion_FP =

distance


criterion_CP =

standard


mu =

     0


theta =

     1


eps_TR =

     1.000000000000000e-04


eps_L =

     1.000000000000000e-03


lSolver =

     1


stratLam =

     1


eps_current =

     1.000000000000000e-05


vstatus =

     0
     0
     0


xstatus =

     1
     1
     1
     1


sstatus =

     1
     1
     1
     1


dstatus =

     0
     0
     0
     0


ndummyY =

     0


sspace_save =

     []


xspace_save =

     []


xfix =

     0
     0
     0


fxmax =

     1.000000000000000e+25


poised_model =

     1


M =

     1     0     0
     0     1     0
     0     0     1


kappa_ill =

     1.000000000000000e+15


kappa_th =

        2000


eps_bnd =

     1.000000000000000e-06


poised =

     1


Y_radius =

     1


c = 

           free: 0
          fixed: 1
    alwaysfixed: 2
             in: 1
            out: 0
         unused: 0
            inY: 1
          dummy: 1
        nodummy: 0


level =

toplevel


whichmodel =

     0


hardcons =

     0


noisy =

     0


scaleX =

     0


scalefacX =

     1     1     1


CNTsin =

     0


shrink_Delta =

     1


scale =

     1
     1
     1
     1


shift_Y =

     1


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
        
        
-----------------------------OUTPUT VALUES -------------------------------------

nit =

     0


i_xbest =

     6


x =

  -0.500000000000000
                   0
   0.500000000000000


fx =

   0.500000000000000


m =

    11


X =

  Columns 1 through 6

   0.500000000000000  -0.500000000000000   0.500000000000000   0.500000000000000  -0.003054026564966  -0.500000000000000
   1.000000000000000   1.000000000000000                   0   1.000000000000000                   0                   0
   0.500000000000000   0.500000000000000   0.500000000000000  -0.500000000000000   0.315521468125839   0.500000000000000

  Columns 7 through 11

  -0.500002253056410  -0.500000286192314  -0.499999900588771  -0.500009741475754  -0.499990015903692
   0.000000442034118  -0.000009988577082   0.000009999638990   0.000000127297985  -0.000000045756443
   0.499990267148911   0.499999612175087   0.500000134714239   0.500002260696260   0.499999438103251


fX =

  Columns 1 through 6

   1.500000000000000   1.500000000000000   0.500000000000000   1.500000000000000   0.099563123926544   0.500000000000000

  Columns 7 through 11

   0.499992520305321   0.499999898467405   0.500000035403030   0.500012002272037   0.499989454106944


ciX =

     []


ceX =

  Columns 1 through 6

   2.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   0.312467441560872   0.000000000000000
   3.000000000000000   2.000000000000000   1.000000000000000                   0  -0.056489622187450  -0.000000000000000

  Columns 7 through 11

  -0.000011543873381  -0.000010662594309   0.000010233764458  -0.000007353481509   0.000009376443116
  -0.000030567541440  -0.000021426821217   0.000020502831926  -0.000002704791005   0.000008206893176


ind_Y =

     6     9    10     7     8    11


Delta =

   5.672511970824019


eps_current =

     1.000000000000000e-05


cur_degree =

     6


fcmodel =

   0.500000000000000  -0.000010001181020   0.000000000002677   0.000010001072100   0.000000000222233   0.000000000197548
   0.000000000000000   0.000010001165044   0.000010001165044   0.000010001165044   0.000000000000000  -0.000000000000000
  -0.000000000000000   0.000010001165043   0.000020002330087   0.000030003495131   0.000000000000000   0.000000000000000


gx =

  -1.000001597463464
   0.000000267687146
   0.999990706694728


normgx =

   1.000001597463464


vstatus =

     0
     0
     0


xstatus =

     0
     0
     0
     0
     0
     1
     1
     1
     1
     1
     1


sstatus =

     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1


dstatus =

     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0


M =

   1.0e+05 *

   3.981667193893546                   0                   0
                   0   6.903934759857267                   0
                   0                   0                   0


ndummyY =

     0


sspace_save =

     []


xspace_save =

     []


msg =

 Convergence in 7 evaluations of the objective function.


CNTsin =

     0


neval =

     7


lm =

                   0
  -0.000005713064576
                   0
   1.999997749517402
  -0.999996152071198


info = 

         g: [3x1 double]
        ai: []
        ae: [2x3 double]
        hl: []
     niter: 4
      flag: 0
    nsimul: [0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        ci: []
        ce: [2x1 double]
         f: 0.500000000000000
      glag: [3x1 double]
     glagn: 6.661338147750939e-16
     feasn: 2.220446049250313e-16
     compl: 0

K>> info.g

ans =

  -1.000001597463464
   0.000000267687146
   0.999990706694728

K>> info.ae

ans =

   1.000000000001365   1.000000000001130   0.999999999995923
   0.999999999985469   1.999999999999835   2.999999999990382

K>> info.ce

ans =

   1.0e-15 *

   0.222044604925031
  -0.111022302462516

K>> info.glag

ans =

   1.0e-15 *

   0.666133814775094
   0.661787688722732
  -0.555111512312578        

@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from evalfgh import evalfgh_
from runtime import matlabarray
from ecdfo_main import *
from ecdfo_global_variables import set_prob, set_threshold
import numpy as np
import helper



class dummyInfo():
    def __init__(self):        
        self.g = matlabarray([
     [0.0],
     [1.0],
     [0.0]])

        self.ai = matlabarray([])
        self.ae = matlabarray([
         [1.0,     1.0,     1.0],
        [1.0 ,    2.0,     3.0]])
        self.hl = matlabarray( [])
        self.niter = 0
        self.flag = 0
        self.nsimul = matlabarray( [0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.ci = matlabarray( [])
        self.ce =  matlabarray([
         [2],
        [3]])
        self.f = 1.500000000000000
        self.glag =  matlabarray([
[-0.333333333013890],
[   0.666666666736111],
[  -0.333333333513889]])

class correctInfo():
    def __init__(self):

        self.g = matlabarray([-1.000001597463464, 0.000000267687146, 0.999990706694728]).T
        self.ai = matlabarray( [])
        self.ae = matlabarray([[1.000000000001365,   1.000000000001130,   0.999999999995923],
   [0.999999999985469,   1.999999999999835,   2.999999999990382]])
        self.hl = matlabarray( [])
        self.niter = 4
        self.flag = 0
        self.nsimul = matlabarray([0, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.ci = matlabarray( [])
        self.ce =  1.0e-15 * matlabarray([ 0.222044604925031, -0.111022302462516]).T
        self.f = 0.500000000000000
        self.glag =    1.0e-15 * matlabarray([0.666133814775094, 0.661787688722732, -0.555111512312578]).T
        self.glagn = 6.661338147750939e-16
        self.feasn = 2.220446049250313e-16
        self.compl = 0

class Test_ecdfo_main(unittest.TestCase):
    """
      Reminder :
      This class is a test for ecdfo_main realizes the optimization loop
    """
    def setUp(self):
        
        set_threshold(1.000000000000000e-08)
        set_prob(3)        

        self.func = evalfgh_
        self.n = 3
        self.nb = 2
        self.mi = 0
        self.me = 2
        self.lm = matlabarray([ 
[                  0.0],
[                   0.0],
[                   0.0],
[  -0.333333332763891],
[  -0.000000000249999]])
        self.nitold = 0
        self.nit = 0
        self.i_xbest = 1
        self.lb = matlabarray([
  [-0.500000000000000],
   [                0.0],
    [            -np.Inf]])
        self.ub = matlabarray([
   [np.Inf],
   [np.Inf],
   [np.Inf]])
        self.m = 4
        self.X = matlabarray([
   [0.500000000000000,  -0.500000000000000,   0.500000000000000,   0.500000000000000],
   [1.000000000000000,   1.000000000000000,                   0.0,   1.000000000000000],
   [0.500000000000000,   0.500000000000000,   0.500000000000000,  -0.500000000000000]])
        self.fX = matlabarray([1.500000000000000,   1.500000000000000,   0.500000000000000,   1.500000000000000])
        self.ciX = matlabarray([])
        self.ceX = matlabarray([
     [2.0,     1.0,     1.0,     1.0],
     [3.0,     2.0,     1.0,     0.0]])
        self.ind_Y = matlabarray([1,     2,     3,     4])
        self.QZ = matlabarray([
     [1.0,     0.0,     0.0,     0.0],
     [0.0,     1.0,     0.0,     0.0],
     [0.0,     0.0,     1.0,     0.0],
     [0.0,     0.0,     0.0,     1.0]])


        self.RZ = matlabarray([
     [1.0,     1.0,     1.0,     1.0],
     [0.0,    -1.0,     0.0,     0.0],
     [0.0 ,    0.0,    -1.0,     0.0],
     [0.0,     0.0,     0.0,    -1.0]])
        self.delta = 1
        self.cur_degree = 4
        self.neval = 0
        self.maxeval = 600
        self.maxit = 600
        self.fcmodel = matlabarray([
   [1.500000000000000,                   0.0,   1.000000000000000,                   0.0],
   [2.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000],
   [3.000000000000000,   1.000000000000000,   2.000000000000000,   3.000000000000000]])
        self.gx = matlabarray([
     [0.0],
     [1.0],
     [0.0]])
        self.normgx = 1.0
        self.show_errg = 0
        self.pquad = 10
        self.pdiag = 7
        self.plin = 4
        self.stallfact = 2.220446049250313e-15
        self.eps_rho = 1.000000000000000e-14
        self.Deltamax = 100000.0
        self.rep_degree = 4
        self.epsilon = 1.000000000000000e-05
        self.verbose = 1
        self.eta1 = 1.000000000000000e-04
        self.eta2 = 0.900000000000000
        self.gamma1 = 0.010000000000000
        self.gamma2 = 0.500000000000000
        self.gamma3 = 2.0
        self.interpol_TR = 1
        self.factor_CV = 100
        self.Lambda_XN = 1.000000000000000e-10
        self.Lambda_CP = 1.200000000000000
        self.factor_FPU = 1.0
        self.factor_FPR = 10.0
        self.Lambda_FP = 1.000000000000000e-10
        self.criterion_S = 'distance'
        self.criterion_FP = 'distance'
        self.criterion_CP = 'standard'
        self.mu = 0.0
        self.theta = 1.0
        self.eps_TR = 1.000000000000000e-04
        self.eps_L = 1.000000000000000e-03
        self.lSolver = 1
        self.stratLam = 1
        self.eps_current = 1.000000000000000e-05
        self.vstatus = matlabarray([
     [0],
     [0],
     [0]])
        self.xstatus = matlabarray([
     [1],
     [1],
     [1],
     [1]])
        self.sstatus = matlabarray([
     [1],
     [1],
     [1],
     [1]])
        self.dstatus = matlabarray([
     [0],
     [0],
     [0],
     [0]])
        self.ndummyY = 0
        self.sspace_save = matlabarray([])
        self.xspace_save = matlabarray([])
        self.xfix = matlabarray([
     [0],
     [0],
     [0]])
        self.fxmax = 1.000000000000000e+25
        self.poised_model = 1
        self.M = matlabarray([
     [1.0,     0.0,     0.0],
     [0.0,     1.0,     0.0],
     [0.0,     0.0,     1.0]])
        self.kappa_ill = 1.000000000000000e+15
        self.kappa_th = 2000
        self.eps_bnd = 1.000000000000000e-06
        self.poised = 1
        self.Y_radius = 1
        self.c = helper.dummyUnionStruct()
        self.c.free = 0
        self.c.fixed = 1
        self.c.alwaysfixed = 2
        self.c.in_ = 1
        self.c.out = 0
        self.c.unused = 0
        self.c.inY = 1
        self.c.dummy = 1
        self.c.nodummy = 0
        self.level = 'toplevel'
        self.whichmodel = 0
        self.hardcons = 0
        self.noisy = 0
        self.scaleX = 0
        self.scalefacX = matlabarray([ 1,     1,     1])
        self.CNTsin = 0
        self.shrink_Delta = 1
        self.scale = matlabarray([
     [1],
     [1],
     [1],
     [1]])
        self.shift_Y = 1

        self.info = dummyInfo()
        
        self.options = helper.dummyOptions()
        self.values = helper.dummyValues()
        
        self.abs_tol=1e-5
        self.rel_tol=1e-5


    def test_ecdfo_main(self):
             
        """
            Test with some values, results compared with matlab
        """
        nit, i_xbest, x, fx, m, X, fX, ciX, ceX, ind_Y, Delta, eps_current, cur_degree, fcmodel, gx, normgx, vstatus, xstatus, sstatus, dstatus, M, ndummyY, sspace_save, xspace_save, msg, CNTsin, neval, lm, info = ecdfo_main_(self.func,self.n,self.nb,self.mi,self.me,self.lm,self.nitold, self.nit, self.i_xbest, self.lb, self.ub, self.m, self.X, self.fX, self.ciX, self.ceX, self.ind_Y, self.QZ, self.RZ, self.delta, self.cur_degree, self.neval, self.maxeval, self.maxit, self.fcmodel, self.gx, self.normgx, self.show_errg, self.pquad, self.pdiag, self.plin, self.stallfact, self.eps_rho, self.Deltamax, self.rep_degree, self.epsilon, self.verbose, self.eta1, self.eta2, self.gamma1, self.gamma2, self.gamma3, self.interpol_TR, self.factor_CV, self.Lambda_XN, self.Lambda_CP, self.factor_FPU, self.factor_FPR, self.Lambda_FP, self.criterion_S, self.criterion_FP, self.criterion_CP, self.mu, self.theta, self.eps_TR, self.eps_L, self.lSolver, self.stratLam, self.eps_current, self.vstatus, self.xstatus, self.sstatus, self.dstatus, self.ndummyY, self.sspace_save, self.xspace_save, self.xfix, self.fxmax, self.poised_model, self.M, self.kappa_ill, self.kappa_th, self.eps_bnd, self.poised, self.Y_radius, self.c, self.level, self.whichmodel, self.hardcons, self.noisy, self.scaleX, self.scalefacX, self.CNTsin, self.shrink_Delta, self.scale, self.shift_Y,self.info,self.options,self.values)
        
        
        correctnit = 0
        correcti_xbest = 6
        correctx = matlabarray([  -0.500000000000000, 0,  0.500000000000000]).T
        correctfx = matlabarray([0.500000000000000])
        correctm = 11
        correctX = matlabarray([
   [0.500000000000000,  -0.500000000000000,   0.500000000000000,   0.500000000000000,  -0.003054026564966,  -0.500000000000000, -0.500002253056410,  -0.500000286192314,  -0.500009741475754,  -0.499999900588771,    -0.499990015903692],
   [1.000000000000000,   1.000000000000000,                   0,   1.000000000000000,                   0,                   0, 0.000000442034118 , -0.000009988577082 ,  0.000000127297985 ,  0.000009999638990,    -0.000000045756443],
   [0.500000000000000,   0.500000000000000,   0.500000000000000,  -0.500000000000000,   0.315521468125839,   0.500000000000000, 0.499990267148911 ,  0.499999612175087 ,  0.500002260696260 ,  0.500000134714239,     0.499999438103251]])
        correctfX = matlabarray([   1.500000000000000,   1.500000000000000,   0.500000000000000,   1.500000000000000,   0.099563123926544,   0.500000000000000, 0.499992520305321,   0.499999898467405,   0.500000035403030,   0.500012002272037,   0.499989454106944])
        correctciX = matlabarray( [])
        correctceX = matlabarray([
  [  2.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,   0.312467441560872,   0.000000000000000, -0.000011543873381,  -0.000010662594309,   -0.000007353481509, 0.000010233764458,     0.000009376443116],
  [ 3.000000000000000 ,  2.000000000000000,   1.000000000000000 ,                  0,  -0.056489622187450,  -0.000000000000000, -0.000030567541440,  -0.000021426821217,  -0.000002704791005,  0.000020502831926 ,    0.000008206893176]])
        correctind_Y = matlabarray([   6,     10, 9,     7,     8,    11])
        correctDelta = 5.672511970824019
        correcteps_current =  1.000000000000000e-05
        correctcur_degree =  6
        correctfcmodel = matlabarray([
   [0.500000000000000,  -0.000010001181020,   0.000000000002677,   0.000010001072100,   0.000000000222233,   0.000000000197548],
   [0.000000000000000,   0.000010001165044,   0.000010001165044,   0.000010001165044,   0.000000000000000,  -0.000000000000000],
  [-0.000000000000000,   0.000010001165043,   0.000020002330087,   0.000030003495131,   0.000000000000000,   0.000000000000000]])
        correctgx = matlabarray([ -1.000001597463464, 0.000000267687146, 0.999990706694728]).T
        correctnormgx = 1.000001597463464
        correctvstatus = matlabarray([  0,  0,  0]).T
        correctxstatus = matlabarray([ 0,  0, 0, 0, 0, 1, 1, 1, 1,  1, 1]).T
        correctsstatus = matlabarray([ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]).T
        correctdstatus = matlabarray([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).T
        correctM =   1.0e+05 * matlabarray([
   [3.981667193893546,                   0,                   0],
   [                0,   6.903934759857267,                   0],
   [                0,                  0 ,                  0]])

        correctndummyY =  0
        correctsspace_save = matlabarray(  [])
        correctxspace_save = matlabarray(  [])
        correctmsg = "Convergence in 7 evaluations of the objective function."
        correctCNTsin =  0
        correctneval = 7
        correctlm = matlabarray([0, -0.000005713064576, 0, 1.999997749517402, -0.999996152071198]).T
        
        correctinfo = correctInfo()
        
        self.assertEqual(nit, correctnit)
        self.assertTrue(compare_matlabarray(correctx, x, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctfx, fx, self.abs_tol, self.rel_tol))
        
        self.assertEqual(correcti_xbest, i_xbest)
        
        #Difference due to different convergence paths in auskommentierten Zeilen
        self.assertTrue(compare_matlabarray(correctm, m, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctX, X, self.abs_tol, 1e-3))
        self.assertEqual(correctciX, ciX)
        self.assertTrue(compare_matlabarray(correctceX, ceX, self.abs_tol, 1e-3))
        self.assertTrue(compare_matlabarray(correctind_Y, ind_Y, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctDelta, Delta, self.abs_tol, self.rel_tol))
        self.assertEqual(correcteps_current, eps_current)
        self.assertEqual(correctcur_degree, cur_degree)
        self.assertTrue(compare_matlabarray(correctfcmodel,fcmodel, self.abs_tol, 1e-4))
        self.assertTrue(compare_matlabarray(correctgx, gx, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctnormgx, normgx, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctvstatus, vstatus, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctxstatus, xstatus, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctdstatus, dstatus, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctsstatus, sstatus, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctM, M, 2000, 1e-2))
        self.assertTrue(compare_matlabarray(correctndummyY, ndummyY, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctsspace_save, sspace_save, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctxspace_save, xspace_save, self.abs_tol, self.rel_tol))
        self.assertEqual(str(correctmsg),str(msg))
        self.assertTrue(compare_matlabarray(correctCNTsin, CNTsin, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctneval, neval, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(correctlm, lm, self.abs_tol, self.rel_tol))

        
        

if __name__ == '__main__':
    unittest.main()
    