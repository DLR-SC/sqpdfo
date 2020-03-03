
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 17:32:15 2014

@author: jaco_da
"""

import unittest
from sqpdfo.runtime import *
from sqpdfo.sqpdfo_global_variables import set_threshold, set_check_condition
from sqpdfo.sqpdfo import sqpdfo_
from numpy import array, zeros, arange, double
from numpy.testing import assert_almost_equal
from sqpdfo import helper
import tests.benchmarks


class Test_run_sqpdfo(unittest.TestCase):
    """
      Reminder :
      This class is a test for run_sqpdfo
      these tests are run with whichmodel=0 and pquad=(n+1)(n+2)/2 in sqpdfo_() !!
    """
    def setUp(self):
        self.abs_tol=1e-5
        self.rel_tol=1e-5
        pass

    def test_run_sqpdfo_prob1(self):
        """
         Test which compare python and matlab results
        """
        args = tests.benchmarks.get(1)
        set_check_condition(0)
        set_threshold(1e-08)

        options = helper.dummyUnionStruct()
        options.hess_approx='model'
        options.bfgs_restart=0
        options.algo_descent='Powell'
        options.tol_grad=1e-05
        options.tol_feas=1e-05
        options.tol_bnds=1e-05
        options.dxmin=1e-6
        options.miter=500
        options.msimul=500
        options.verbose=0
        options.whichmodel = 'subbasis'
        options.final_degree = 'quadratic'

        x,lm,info=sqpdfo_(options, *args)

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
        self.assertTrue(compare_array(info.glag, 1e-10*array([0.76387451919401,-0.381947806715743]), self.abs_tol, self.rel_tol))
        
    def test_run_sqpdfo_prob2(self):
        """
         Test which compare python and matlab results
        """
        args = tests.benchmarks.get(2)
        set_check_condition(0)
        set_threshold(1e-08)

        options = helper.dummyUnionStruct()
        options.hess_approx='model'
        options.bfgs_restart=0
        options.algo_descent='Powell'
        options.tol_grad=1e-05
        options.tol_feas=1e-05
        options.tol_bnds=1e-05
        options.dxmin=1e-6
        options.miter=500
        options.msimul=500
        options.verbose=0
        options.whichmodel = 'subbasis'
        options.final_degree = 'quadratic'

        x,lm,info=sqpdfo_(options, *args)

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
        self.assertTrue(compare_array(info.glag, 1e-05*array([  -0.560839309260430,0.560883753109032]), self.abs_tol, self.rel_tol))
#        
    def test_run_sqpdfo_prob3(self):
        """
         Test which compare python and matlab results
        """
        set_check_condition(0)
        args = tests.benchmarks.get(3)
        set_threshold(1e-08)

        options = helper.dummyUnionStruct()
        options.hess_approx='model'
        options.bfgs_restart=0
        options.algo_descent='Powell'
        options.tol_grad=1e-05
        options.tol_feas=1e-05
        options.tol_bnds=1e-05
        options.dxmin=1e-6
        options.miter=500
        options.msimul=500
        options.verbose=0
        options.whichmodel = 'subbasis'
        options.final_degree = 'quadratic'

        x,lm,info=sqpdfo_(options, *args)

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
        self.assertTrue(compare_array(info.glag, 1e-07*array([0.062859997207454,  0.188580686583393,-0.188580611126810]), self.abs_tol, self.rel_tol))

    def test_run_sqpdfo_prob4(self):
        """
         Test which compare python and matlab results
        """
        set_check_condition(0)
        args = tests.benchmarks.get(4)
        set_threshold(1e-08)

        options = helper.dummyUnionStruct()
        options.hess_approx='model'
        options.bfgs_restart=0
        options.algo_descent='Powell'
        options.tol_grad=1e-05
        options.tol_feas=1e-05
        options.tol_bnds=1e-05
        options.dxmin=1e-6
        options.miter=500
        options.msimul=500
        options.verbose=0
        options.whichmodel = 'subbasis'
        options.final_degree = 'quadratic'

        x,lm,info=sqpdfo_(options, *args)

        self.assertTrue(compare_array(x, array([[ -0.5,  0.0,   0.5,   1.0]]), self.abs_tol, 1e-5))
        self.assertTrue(compare_array(lm, array([[ 0,0,0,0,1.999999758015728,-0.999999892175830,-0.333333335490867]]), self.abs_tol, 1e-6))
        self.assertTrue(compare_array(info.g, array([[   -1.0,   0.0 ,    1.0,   1.0]]), 1e-4, 1e-4))
        self.assertTrue(compare_array(info.ae, array([[1.00000000000119,	1.00000000000153,	1.00000000000242,	1.09566597368092e-12],[1.00000000000310,	1.99999999999671,	3.00000000000154,	5.36746539189640e-12],[1.43689077548211e-12,	3.68227391508966e-12,	-3.71857711253322e-12,	2.99999999015165]]), self.abs_tol, 1e-6))
        self.assertTrue(compare_array(info.ce,  array([ 0., 0., 0.]), 1e-4, 1e-4))
        self.assertEqual(info.flag,0)
        self.assertAlmostEqual(info.f, 1.5, places=6)
        self.assertEqual(info.compl,0)
        self.assertTrue(compare_array(info.glag, array([0.,0.,0.,0.]), 1e-4, 1e-4))

    def test_run_sqpdfo_prob5(self):
        """
         Test which compare python and matlab results
        """
        set_check_condition(0)
        set_threshold(1e-08)

        options = helper.dummyUnionStruct()
        options.hess_approx='model'
        options.bfgs_restart=0
        options.algo_descent='Powell'
        options.tol_grad=1e-05
        options.tol_feas=1e-05
        options.tol_bnds=1e-05
        options.dxmin=1e-6
        options.miter=500
        options.msimul=500
        options.verbose=0
        options.whichmodel = 'subbasis'
        options.final_degree = 'quadratic'

        x,lm,info=sqpdfo_(options, *tests.benchmarks.get(5))

        assert_almost_equal(x, array([   -1.717135373541669,1.595700197732079,1.827260995778992,0.763703525309065,0.763584463389690]), decimal=3)
        assert_almost_equal(lm, array([[ 0,0,0,0,0, 0.040162755804678,-0.037957678618516, 0.005222725990309]]).T, decimal=3)
        assert_almost_equal(info.g.T, array([[ 0.091732656086263,  -0.098713648653029,  -0.086204099493362,     -0.206254630841206,  -0.206286791111375]]), decimal=4)
        assert_almost_equal(info.ae, array([[-3.43427074714240,	3.19140039544429,	3.65452199155318,	1.52740705050108	,1.52716892687853],
[-1.50249321407029e-11,	1.82726099577928	,1.59570019773225	,-3.81792231695518	,-3.81851762654697],
[8.84566167319830,	7.63877736315284,	2.39648693851626e-11	,3.91573621289396e-11,	1.95315487528657e-11]]), decimal=3)
        assert_almost_equal(info.ce.T,    1.0e-07 *array([[  0.661051284822634,  -0.005370299760443, -0.042228007757217]]), decimal=5)
        self.assertEqual(info.flag,0)
        assert_almost_equal(info.f,      0.053949845718415, decimal=7)
        self.assertEqual(info.compl,0)
        assert_almost_equal(info.glag, array([[0.,0.,0.,0.,0.]]).T, decimal=4)
    
if __name__ == '__main__':
    unittest.main()