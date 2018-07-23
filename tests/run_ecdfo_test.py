
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
import unittest
from numpy import array, zeros, arange, double


class Test_run_ecdfo(unittest.TestCase):
    """
      Reminder :
      This class is a test for run_ecdfo
      these tests are run with whichmodel=0 and pquad=(n+1)(n+2)/2 in ecdfo_() !!
    """
    def setUp(self):
        self.abs_tol=1e-5
        self.rel_tol=1e-8
        pass

    def test_run_ecdfo_prob1(self):
        """
         Test which compare python and matlab results
        """
        set_prob(1)
        set_check_condition(0)
        prob=get_prob()
        options = helper.dummyUnionStruct()
        options.tol=zeros(3)

        x,lb,ub,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(prob,nargout=13)
        set_threshold(1e-08)
        options.hess_approx='model'
        options.bfgs_restart=0
        options.algo_descent='Powell'
        options.tol[0]=1e-05
        options.tol[1]=1e-05
        options.tol[2]=1e-05
        options.dxmin=dxmin
        options.miter=500
        options.msimul=500
        options.verbose=0
        lm=array([])
        x,lm,info=ecdfo_(x,lm,lb,ub,options)

        self.assertTrue(compare_array(x, array([[1.950000000000000,0.262499999999991]]), 1e-5, 1e-5))
        self.assertTrue(compare_array(lm, array([[-0.6375,0,  0.7375]]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(info.g, array([[-0.1, -2.95]]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(info.ae, array([[ 1., 4.]]), self.abs_tol, self.rel_tol))
        self.assertEqual(info.niter,4)
        self.assertAlmostEqual(double(info.ce), -3.819167204710539e-14,places=5)
        self.assertEqual(info.flag,0)
        self.assertTrue(compare_array(info.nsimul, array([[0, 9, 0, 0 ]]), self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(info.f,-3.909687499999972,places=10)
        self.assertEqual(info.compl,0)
        self.assertAlmostEqual(dcimin, 1.490116119384766e-08,places=10)
        self.assertTrue(compare_array(info.glag, 1e-10*array([0.76387451919401,-0.381947806715743]), self.abs_tol, self.rel_tol))
        
    def test_run_ecdfo_prob2(self):
        """
         Test which compare python and matlab results
        """
        set_prob(2) 
        set_check_condition(0)
        prob=get_prob()
        options = helper.dummyUnionStruct()
        options.tol=zeros(3)

        x,lb,ub,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(prob,nargout=13)
        set_threshold(1e-08)
        options.hess_approx='model'
        options.bfgs_restart=0
        options.algo_descent='Powell'
        options.tol[0]=1e-05
        options.tol[1]=1e-05
        options.tol[2]=1e-05
        options.dxmin=dxmin
        options.miter=500
        options.msimul=500
        options.verbose=0
        lm=array([])
        x,lm,info=ecdfo_(x,lm,lb,ub,options)
        self.assertTrue(compare_array(x, array([[   0.333326758778846,  0.666659126169760]]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(lm, array([[0,0,    -1.333312643708242]]), self.abs_tol, self.rel_tol))

        self.assertTrue(compare_array(info.g, array([[    1.333307035124744,    1.333318252334031]]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(info.ae, array([[ 1.000000000014289,  1.000000000017430]]), self.abs_tol, self.rel_tol))
        self.assertEqual(info.niter,5)
        self.assertAlmostEqual(double(info.ce), -1.411505139448099e-05,places=10)
        self.assertEqual(info.flag,0)
        self.assertTrue(compare_array(info.nsimul, array([[0, 11, 0, 0]]), self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(info.f, 0.666647846741449,places=10)
        self.assertEqual(info.compl,0)
        self.assertAlmostEqual(dcimin, 1.490116119384766e-08,places=10)
        self.assertTrue(compare_array(info.glag, 1e-05*array([  -0.560839309260430,0.560883753109032]), self.abs_tol, self.rel_tol))
#        
    def test_run_ecdfo_prob3(self):
        """
         Test which compare python and matlab results
        """
        set_check_condition(0)
        set_prob(3) 
        prob=get_prob()
        options = helper.dummyUnionStruct()
        options.tol=zeros(3)

        x,lb,ub,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(prob,nargout=13)
        set_threshold(1e-08)
        options.hess_approx='model'
        options.bfgs_restart=0
        options.algo_descent='Powell'
        options.tol[0]=1e-05
        options.tol[1]=1e-05
        options.tol[2]=1e-05
        options.dxmin=dxmin
        options.miter=500
        options.msimul=500
        options.verbose=0
        lm=array([])
        x,lm,info=ecdfo_(x,lm,lb,ub,options)

        self.assertTrue(compare_array(x, array([[-0.5,0,0.5]]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(lm, array([[ 0, -0.000005713064576 ,0,   1.999997749517402,  -0.999996152071198]]), self.abs_tol, self.abs_tol))
        self.assertTrue(compare_array(info.g, array([[ -1.000001597463464,  0.000000267687146 ,  0.999990706694728]]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(info.ae, array([[ 1.000000000001365,  1.000000000001130, 0.999999999995923],[   0.999999999985469 , 1.999999999999835,    2.999999999990382 ]]), self.abs_tol, self.rel_tol))
        self.assertEqual(info.niter,4)
        self.assertTrue(compare_array(info.ce, array([[0.222044604925031e-15, -0.111022302462516e-15]]),self.abs_tol, self.rel_tol))
        self.assertEqual(info.flag,0)
        self.assertTrue(compare_array(info.nsimul, array([[0, 11, 0, 0]]), self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(info.f,   0.500000000000000,places=10)
        self.assertEqual(info.compl,0)
        self.assertAlmostEqual(dcimin, 1.490116119384766e-08,places=10)
        self.assertTrue(compare_array(info.glag, 1e-07*array([0.062859997207454,  0.188580686583393,-0.188580611126810]), self.abs_tol, self.rel_tol))

    # def test_run_ecdfo_prob4(self):
    #     """
    #      Test which compare python and matlab results
    #     """
    #     set_check_condition(1)
    #     set_prob(4)
    #     prob=get_prob()
    #     options = helper.dummyUnionStruct()
    #     options.tol=zeros(3)
    #
    #     x,lb,ub,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(prob,nargout=13)
    #     set_threshold(1e-08)
    #     options.hess_approx='model'
    #     options.bfgs_restart=0
    #     options.algo_descent='Powell'
    #     options.tol[0]=1e-05
    #     options.tol[1]=1e-05
    #     options.tol[2]=1e-05
    #     options.dxmin=dxmin
    #     options.miter=500
    #     options.msimul=500
    #     options.verbose=1
    #     lm=array([])
    #     x,lm,info=ecdfo_(x,lm,lb,ub,options)
    #
    #     self.assertTrue(compare_array(x, array([[ -0.499998511434003,  -0.000002977131994,   0.500001488565997,   0.999999998348743]]), self.abs_tol, 1e-6))
    #     self.assertTrue(compare_array(lm, array([[ 0,0,0,0,1.999999758015728,-0.999999892175830,-0.333333335490867]]), self.abs_tol, 1e-6))
    #     self.assertTrue(compare_array(info.g, array([[   -0.999997022863534,   -0.000005954272170 ,    1.000002977136998,   0.999999999997250]]), self.abs_tol, 1e-6))
    #     self.assertTrue(compare_array(info.ae, array([[1.00000000000119,	1.00000000000153,	1.00000000000242,	1.09566597368092e-12],[1.00000000000310,	1.99999999999671,	3.00000000000154,	5.36746539189640e-12],[1.43689077548211e-12,	3.68227391508966e-12,	-3.71857711253322e-12,	2.99999999015165]]), self.abs_tol, 1e-6))
    #     self.assertEqual(info.niter,21)
    #     self.assertTrue(compare_array(info.ce,    1.0e-08 *array([ 0 , 0.000000022204460, -0.495377106002337]),self.abs_tol, self.rel_tol))
    #     self.assertEqual(info.flag,0)
    #     self.assertTrue(compare_array(info.nsimul, array([[0, 38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0 ]]), self.abs_tol, self.rel_tol))
    #     self.assertAlmostEqual(info.f,   1.499999998362038,places=9)
    #     self.assertEqual(info.compl,0)
    #     self.assertAlmostEqual(dcimin, 1.490116119384766e-08,places=10)
    #     self.assertTrue(compare_array(info.glag, 1e-05*array([0.284297517494370,-0.598060297747812,0.305862976301974,-0.000319574466889]), self.abs_tol, self.rel_tol))
#
    def test_run_ecdfo_prob5(self):
        """
         Test which compare python and matlab results
        """
        set_check_condition(1)
        set_prob(5) 
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
        x,lm,info=ecdfo_(x,lm,lb,ub,options)
        
        self.assertTrue(compare_array(x, array([[   -1.717135373541669,1.595700197732079,1.827260995778992,0.763703525309065,0.763584463389690]]), 1e-3, 1e-4))
        self.assertTrue(compare_array(lm, array([[ 0,0,0,0,0, 0.040162755804678,-0.037957678618516, 0.005222725990309]]), 1e-3, 1e-4))
        self.assertTrue(compare_array(info.g, array([[ 0.091732656086263,  -0.098713648653029,  -0.086204099493362,     -0.206254630841206,  -0.206286791111375]]), 1e-4, 1e-4))
        self.assertTrue(compare_array(info.ae, array([[-3.43427074714240,	3.19140039544429,	3.65452199155318,	1.52740705050108	,1.52716892687853],
[-1.50249321407029e-11,	1.82726099577928	,1.59570019773225	,-3.81792231695518	,-3.81851762654697],
[8.84566167319830,	7.63877736315284,	2.39648693851626e-11	,3.91573621289396e-11,	1.95315487528657e-11]]), 1,1e-4))
        self.assertTrue(compare_array(info.ce,    1.0e-07 *array([  0.661051284822634,  -0.005370299760443, -0.042228007757217]),self.abs_tol, self.rel_tol))
        self.assertEqual(info.flag,0)
        self.assertAlmostEqual(info.f,      0.053949845718415,places=7)
        self.assertEqual(info.compl,0)
        self.assertAlmostEqual(dcimin, 1.490116119384766e-08,places=10)
        self.assertTrue(compare_array(info.glag, 1e-05*array([0.134582429518748,-0.155845706949209,0.249955931404255,0.971384006251408,-0.941356094893986]), self.abs_tol, self.rel_tol))
    
if __name__ == '__main__':
    unittest.main()