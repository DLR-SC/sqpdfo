# -*- coding: utf-8 -*-
from runtime import *
import numpy as np
	
class dummyUnionStruct():
	pass
	
class dummyOptions():
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
	def __init__(self):
		self.success = 0
		self.fail_on_argument = 1
		self.fail_on_problem = 2
		self.fail_on_simul = 3
		self.stop_on_simul = 4
		self.stop_on_max_iter = 5
		self.stop_on_max_simul = 6
		self.stop_on_dxmin = 7
		self.fail_on_non_decrease = 8
		self.fail_on_ascent_dir = 9
		self.fail_on_max_ls_iter = 10
		self.fail_on_ill_cond = 11
		self.stop_on_small_trust_region = 15
		self.fail_on_null_step = 20
		self.fail_on_infeasible_QP = 21
		self.fail_on_unbounded_QP = 22
		self.fail_strange = 99
		self.nsimultype = 16
		self.max_null_steps = 1
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
