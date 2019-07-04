import unittest
import sys
import warnings
warnings.simplefilter('ignore', FutureWarning)
sys.path.append("../")

# Import 
import ecdfo_on_cuter_equalities_test
#import ecdfo_main_test

loader = unittest.TestLoader() 

# Tests
suite = loader.loadTestsFromTestCase(ecdfo_on_cuter_equalities_test.Test_ecdfo_on_cuter_equalities)
#suite.addTests(loader.loadTestsFromTestCase(ecdfo_main_test.Test_ecdfo_main))

unittest.TextTestRunner(verbosity=1).run(suite)
