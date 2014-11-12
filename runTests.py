# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 14:31:47 2014

@author: jaco_da
"""
import unittest
import ecdfo_check_cond_test
import ecdfo_check_convex_test
import blls_test
import sqplab_badsimul_test
import sqplab_bfgs_test
import sqplab_checkoptions_test
import sqplab_options_test
import sqplab_tcg_test

loader = unittest.TestLoader() 

suite = loader.loadTestsFromTestCase(ecdfo_check_convex_test.Test_ecdfo_check_convex)

suite.addTests(loader.loadTestsFromTestCase(ecdfo_check_cond_test.Test_ecdfo_check_cond))
suite.addTests(loader.loadTestsFromTestCase(blls_test.Test_blls))
suite.addTests(loader.loadTestsFromTestCase(sqplab_badsimul_test.Test_sqplab_badsimul))
suite.addTests(loader.loadTestsFromTestCase(sqplab_bfgs_test.Test_sqplab_bfgs))
suite.addTests(loader.loadTestsFromTestCase(sqplab_checkoptions_test.Test_sqplab_checkoptions))
suite.addTests(loader.loadTestsFromTestCase(sqplab_options_test.Test_sqplab_options))
suite.addTests(loader.loadTestsFromTestCase(sqplab_tcg_test.Test_sqplab_tcg))


unittest.TextTestRunner(verbosity=2).run(suite)