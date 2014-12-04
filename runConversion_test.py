# -*- coding: utf-8 -*-
"""
Created on Wed Dec 03 13:17:11 2014

@author: jaco_da
"""
import unittest
import bcdfo_evalL_test
import bcdfo_find_new_yj_test
import bcdfo_hessP_test
import bcdfo_include_in_Y_test
import bcdfo_repair_Y_test
import bcdfo_replace_in_Y_test
import bcdfo_solve_TR_MS_test

import bcdfo_build_QR_of_Y_test
import bcdfo_solve_TR_MS_bc_test
import bcdfo_evalZ_test
import bcdfo_poisedness_Y_test
import bcdfo_augment_Y_test
import bcdfo_gradP_test
import bcdfo_evalP_test

loader = unittest.TestLoader() 

suite = loader.loadTestsFromTestCase(bcdfo_evalL_test.Test_bcdfo_evalL)

#ECDFO Tests
suite.addTests(loader.loadTestsFromTestCase(bcdfo_find_new_yj_test.Test_bcdfo_find_new_yj))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_hessP_test.Test_bcdfo_hessP))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_include_in_Y_test.Test_bcdfo_include_in_Y))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_repair_Y_test.Test_bcdfo_repair_Y))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_replace_in_Y_test.Test_bcdfo_replace_in_Y))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_solve_TR_MS_test.Test_bcdfo_solve_TR_MS))


#BCDFO Tests
suite.addTests(loader.loadTestsFromTestCase(bcdfo_build_QR_of_Y_test.Test_bcdfo_build_QR_of_Y))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_solve_TR_MS_bc_test.Test_bcdfo_solve_TR_MS_bc))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_evalZ_test.Test_bcdfo_evalZ))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_poisedness_Y_test.Test_bcdfo_poisedness_Y))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_augment_Y_test.Test_bcdfo_augment_Y))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_gradP_test.Test_bcdfo_gradP))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_evalP_test.Test_bcdfo_evalP))

unittest.TextTestRunner(verbosity=2).run(suite)