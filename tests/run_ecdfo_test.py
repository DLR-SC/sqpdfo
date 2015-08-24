
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 17:32:15 2014

@author: jaco_da
"""
import sys
sys.path.append("../")
import helper
from runtime import *
from ecdfo_init_prob import ecdfo_init_prob_
from ecdfo_global_variables import set_prob, set_threshold,get_prob, set_check_condition
from ecdfo import ecdfo_
from evalfgh import evalfgh_
import unittest
from numpy import array, zeros, arange
#from run_ecdfo import *
#import numpy as np
#import helper

class Test_run_ecdfo(unittest.TestCase):

    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-5
        self.rel_tol=1e-8
        pass

    def test_run_ecdfo_prob1(self):
        set_prob(1)
        set_check_condition(0)
        prob=get_prob()
        options = helper.dummyUnionStruct()
        options.tol=zeros(3)

        x,lx,ux,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(prob,nargout=13)
        lb=zeros_(n,1)
        ub=zeros_(n,1)
        lb[arange(0,n)]=lx
        ub[arange(0,n)]=ux
        if mi:
            lb[arange(n,n + mi)]=li
            ub[arange(n,n + mi)]=ui
        set_threshold(1e-08)
        options.algo_method='quasi-Newton'
        options.algo_globalization='trust regions'
        options.hess_approx='model'
        options.bfgs_restart=0
        options.algo_descent='Powell'
        if nb + mi + me == 0:
            options.algo_descent='Wolfe'
        options.tol[0]=1e-05
        options.tol[1]=1e-05
        options.tol[2]=1e-05
        options.dxmin=dxmin
        options.miter=500
        options.msimul=500
        options.verbose=0
        lm=array([])
        x,lm,info=ecdfo_(evalfgh_,x,lm,lb,ub,options,nargout=3)
        
        self.assertTrue(compare_array(x, array([[1.950000000000000,0.262499999999991]]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(lm, array([[-0.637499999995624,0,  0.737500000002291]]), self.abs_tol, self.rel_tol))

        self.assertTrue(compare_array(info.g, array([[  -0.100000000005919, -2.950000000008655]]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(info.ae, array([[ 0.999999999998983, 3.999999999999313]]), self.abs_tol, self.rel_tol))
        self.assertEqual(info.niter,4)
        self.assertAlmostEqual(info.ce, -3.819167204710539e-14,places=10)
        self.assertEqual(info.flag,0)
        self.assertTrue(compare_array(info.nsimul, array([[0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0 ]]), self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(info.f,-3.909687499999972,places=10)
        self.assertEqual(info.compl,0)
        
    def test_run_ecdfo_prob2(self):
        set_prob(2) 
        set_check_condition(0)
        prob=get_prob()
        options = helper.dummyUnionStruct()
        options.tol=zeros(3)

        x,lx,ux,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(prob,nargout=13)
        lb=zeros_(n,1)
        ub=zeros_(n,1)
        lb[arange(0,n)]=lx
        ub[arange(0,n)]=ux
        if mi:
            lb[arange(n ,n + mi)]=li
            ub[arange(n,n + mi)]=ui
        set_threshold(1e-08)
        options.algo_method='quasi-Newton'
        options.algo_globalization='trust regions'
        options.hess_approx='model'
        options.bfgs_restart=0
        options.algo_descent='Powell'
        if nb + mi + me == 0:
            options.algo_descent='Wolfe'
        options.tol[0]=1e-05
        options.tol[1]=1e-05
        options.tol[2]=1e-05
        options.dxmin=dxmin
        options.miter=500
        options.msimul=500
        options.verbose=0
        lm=array([])
        x,lm,info=ecdfo_(evalfgh_,x,lm,lb,ub,options,nargout=3)
        
        self.assertTrue(compare_array(x, array([[   0.333326758778846,  0.666659126169760]]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(lm, array([[0,0,    -1.333312643708242]]), self.abs_tol, self.rel_tol))

        self.assertTrue(compare_array(info.g, array([[    1.333307035124744,    1.333318252334031]]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(info.ae, array([[ 1.000000000014289,  1.000000000017430]]), self.abs_tol, self.rel_tol))
        self.assertEqual(info.niter,6)
        self.assertAlmostEqual(info.ce, -1.411505139448099e-05,places=10)
        self.assertEqual(info.flag,0)
        self.assertTrue(compare_array(info.nsimul, array([[0, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0 ]]), self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(info.f, 0.666647846741449,places=10)
        self.assertEqual(info.compl,0)
#        
    def test_run_ecdfo_prob3(self):
        set_check_condition(0)
        set_prob(3) 
        prob=get_prob()
        options = helper.dummyUnionStruct()
        options.tol=zeros(3)

        x,lx,ux,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(prob,nargout=13)
        lb=zeros_(n,1)
        ub=zeros_(n,1)
        lb[arange(0,n)]=lx
        ub[arange(0,n)]=ux
        if mi:
            lb[arange(n,n + mi)]=li
            ub[arange(n,n + mi)]=ui
        set_threshold(1e-08)
        options.algo_method='quasi-Newton'
        options.algo_globalization='trust regions'
        options.hess_approx='model'
        options.bfgs_restart=0
        options.algo_descent='Powell'
        if nb + mi + me == 0:
            options.algo_descent='Wolfe'
        options.tol[0]=1e-05
        options.tol[1]=1e-05
        options.tol[2]=1e-05
        options.dxmin=dxmin
        options.miter=500
        options.msimul=500
        options.verbose=0
        lm=array([])
        x,lm,info=ecdfo_(evalfgh_,x,lm,lb,ub,options,nargout=3)
        
        self.assertTrue(compare_array(x, array([[-0.5,0,0.5]]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(lm, array([[ 0, -0.000005713064576 ,0,   1.999997749517402,  -0.999996152071198]]), self.abs_tol, self.rel_tol))

        self.assertTrue(compare_array(info.g, array([[ -1.000001597463464,  0.000000267687146 ,  0.999990706694728]]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(info.ae, array([[ 1.000000000001365,  1.000000000001130, 0.999999999995923],[   0.999999999985469 , 1.999999999999835,    2.999999999990382 ]]), self.abs_tol, self.rel_tol))
        self.assertEqual(info.niter,4)
        self.assertTrue(compare_array(info.ce, array([[0.222044604925031e-15, -0.111022302462516e-15]]),self.abs_tol, self.rel_tol))
        self.assertEqual(info.flag,0)
        self.assertTrue(compare_array(info.nsimul, array([[0, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0 ]]), self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(info.f,   0.500000000000000,places=10)
        self.assertEqual(info.compl,0)

if __name__ == '__main__':
    unittest.main()