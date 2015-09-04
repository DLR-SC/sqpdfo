# -*- coding: utf-8 -*-
from runtime import *
import numpy as np
	
class dummyUnionStruct():
      """
      Creates an empty object which can for instance later be used options, info in the ECDFO algorithm
      """
      pass
	
class dummyOptions():
      """
      Those are some basics options for the ECDFO algorithm : they are used for the ECDFO tests
      """
      def __init__(self):
		self.algo_method = 101
		self.algo_globalization = 112
		self.hess_approx = 131
		self.bfgs_restart = 0
		self.algo_descent = 120
		self.tol = np.array( [1.0e-05, 1.0e-05, 1.0e-05])
		self.dxmin = 1.0e-06
		self.miter = 500
		self.msimul = 500
		self.verbose = 0#2 set to a low level to avoid messages : matlab is the reference to debug
		self.fout = 1
		self.inf = np.Inf#Inf
		self.df1 = 0

class dummyValues():
      """
      Those are some basics options for the ECDFO algorithm : they are used for the ECDFO tests
      """
      def __init__(self):
           # Define ecdfo constant values
        
           self.success                =   0; # solution found
           self.fail_on_argument       =   1; # an argument is wrong
           self.fail_on_problem        =   2; # unaccepted problem structure
           self.fail_on_simul          =   3; # error on a simulation
           self.stop_on_simul          =   4; # stop required by the simulator
           self.stop_on_max_iter       =   5; # max iterations
           self.stop_on_max_simul      =   6; # max simulations
           self.stop_on_dxmin          =   7; # stop on dxmin
           self.fail_on_non_decrease   =   8; # the merit function no longer decrease
           self.fail_on_ascent_dir     =   9; # nondescent direction in linesearch
           self.fail_on_max_ls_iter    =  10; # too many stepsize trials in linesearch
           self.fail_on_ill_cond       =  11; # ill-conditioning
          
           self.stop_on_small_trust_region = 15;
        
           self.fail_on_null_step      =  20; # null step d is solution of 'values.max_null_steps' QPs
           self.fail_on_infeasible_QP  =  21; # infeasible QP
           self.fail_on_unbounded_QP   =  22; # unbounded QP
        
           self.fail_strange           =  99; # should not have happened, call a guru
         
           self.nsimultype             =  16; # nb of simulation types
           self.max_null_steps         =   1; # maximum nbr of null QP steps

           self.newton = 100
           self.quasi_newton = 101
           self.cheap_quasi_newton = 102
           self.unit_stepsize = 110
           self.linesearch = 111
           self.trust_regions = 112
           self.powell = 120
           self.wolfe = 121
           self.bfgs = 130
           self.model = 131
           self.dline = '--------------------------------------------------------------------------------------'
           self.eline = '======================================================================================'
           self.sline = '**************************************************************************************'
