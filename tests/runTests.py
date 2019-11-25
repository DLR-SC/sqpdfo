import unittest
import sys
import warnings
warnings.simplefilter('ignore', FutureWarning)
sys.path.append("../")

#Import runtime  test
import runtime_test

#Import of BCDFO Tests
import bcdfo_evalZ_test
import bcdfo_build_QR_of_Y_test
import bcdfo_evalP_test
import bcdfo_evalL_test
import bcdfo_computeLj_test
import bcdfo_computeP_test
import bcdfo_solve_TR_MS_test
import bcdfo_augment_Y_test
import bcdfo_gradP_test
import bcdfo_hessP_test
import bcdfo_find_new_yj_test
import bcdfo_replace_in_Y_test
import bcdfo_include_in_Y_test
import bcdfo_poisedness_Y_test
import bcdfo_projgrad_test
import bcdfo_repair_Y_test
import bcdfo_solve_TR_MS_bc_test


#Import of blls, and sqpdfo_check_convex/check_cond Tests

import sqpdfo_check_cond_test
import sqpdfo_check_convex_test
import blls_test
import sqpdfo_bfgs_update_test
import sqpdfo_compute_multiplier_test
import sqpdfo_options_test
import sqpdfo_truncated_cg_test

#Import of sqpdfo ests which do not print anything on the screen

import sqpdfo_augmX_evalf_test
import sqpdfo_func_test
import sqpdfo_solve_TR_bc_test
import sqpdfo_init_prob_test
import sqpdfo_evalfgh_test
import sqpdfo_optimality_test
import sqpdfo_swap_in_Y_test
import sqpdfo_find_smallf_test
import sqpdfo_computeHessian_test

#Import of sqpdfo tests which print things on the screen

import sqpdfo_iter_printout_test
import sqpdfo_prelim_test
import sqpdfo_main_test
import sqpdfo_finish_test
import sqpdfo_test
import run_sqpdfo_test

loader = unittest.TestLoader() 

#NB: tests are added in the suite in a logical order, such that one test is added only after all of its dependencies.

#BCDFO Tests
suite = loader.loadTestsFromTestCase(bcdfo_evalZ_test.Test_bcdfo_evalZ)
suite.addTests(loader.loadTestsFromTestCase(bcdfo_build_QR_of_Y_test.Test_bcdfo_build_QR_of_Y))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_evalP_test.Test_bcdfo_evalP))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_evalL_test.Test_bcdfo_evalL))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_computeLj_test.Test_bcdfo_computeLj))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_computeP_test.Test_bcdfo_computeP))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_solve_TR_MS_test.Test_bcdfo_solve_TR_MS))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_augment_Y_test.Test_bcdfo_augment_Y))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_gradP_test.Test_bcdfo_gradP))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_hessP_test.Test_bcdfo_hessP))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_find_new_yj_test.Test_bcdfo_find_new_yj))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_replace_in_Y_test.Test_bcdfo_replace_in_Y))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_include_in_Y_test.Test_bcdfo_include_in_Y))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_poisedness_Y_test.Test_bcdfo_poisedness_Y))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_projgrad_test.Test_bcdfo_projgrad))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_repair_Y_test.Test_bcdfo_repair_Y))
suite.addTests(loader.loadTestsFromTestCase(bcdfo_solve_TR_MS_bc_test.Test_bcdfo_solve_TR_MS_bc))
#
# blls, and sqpdfo_check_convex/check_cond Tests
#
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_check_convex_test.Test_sqpdfo_check_convex))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_check_cond_test.Test_sqpdfo_check_cond))
suite.addTests(loader.loadTestsFromTestCase(blls_test.Test_blls))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_bfgs_update_test.Test_sqpdfo_bfgs_update))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_compute_multiplier_test.Test_sqpdfo_compute_multiplier))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_options_test.Test_sqpdfo_options))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_truncated_cg_test.Test_sqpdfo_truncated_cg))
#
##sqpdfo Tests which do not print anything on the screen
#
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_augmX_evalf_test.Test_sqpdfo_augmX_evalf))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_func_test.Test_sqpdfo_func))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_solve_TR_bc_test.Test_sqpdfo_solve_TR_bc))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_init_prob_test.Test_sqpdfo_init_prob))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_evalfgh_test.Test_sqpdfo_evalfgh))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_optimality_test.Test_sqpdfo_optimality))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_swap_in_Y_test.Test_sqpdfo_swap_in_Y))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_find_smallf_test.Test_sqpdfo_find_smallf))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_computeHessian_test.Test_sqpdfo_computeHessian))
#
##sqpdfo Tests which print text on the screen (at least with some verbose >0)
#
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_iter_printout_test.Test_sqpdfo_iter_printout))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_prelim_test.Test_sqpdfo_prelim))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_main_test.Test_sqpdfo_main))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_finish_test.Test_sqpdfo_finish))
suite.addTests(loader.loadTestsFromTestCase(sqpdfo_test.Test_sqpdfo))
suite.addTests(loader.loadTestsFromTestCase(run_sqpdfo_test.Test_run_sqpdfo))


unittest.TextTestRunner(verbosity=1).run(suite)
