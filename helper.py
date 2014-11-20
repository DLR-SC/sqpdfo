# -*- coding: utf-8 -*-
"""
Created on Thu Nov 06 14:50:24 2014

@author: jaco_da
"""
from runtime import *

#class dummyOptions:
#	def __init__(self):
#		self.verbose = 4
#		self.fout = 42

#so far for bad sqplab_badsimul Unittest
#class dummyInfo():
#	def __init__(self):
#		self.flag = False
		
#class dummyValues():
	#def __inite__(self):
	#	self.stop_on_simul = False
	#	self.fail_on_simul = False
		
class generalDecorator():
	def convert(self, a):
		if a.shape == (1,1):
			print "shape (1,1)"
			b = np.zeros(1)
			b[0] = a[0]
			return b[0]

		if a.shape[1] == 1:
			a = a.T			
			#print "transposed a", a
		if a.shape[0] == 1:
			b = np.zeros(a.shape)#[0]
			#print "b", b
			#print "b.shape", b.shape
			for i in range(b.shape[1]):#len(b)):
				b[0][i] = a[i+1]
		else:
			b = np.zeros(a.shape)#[0]
			#print "b", b
			#print "b.shape", b.shape
			for i in range(b.shape[0]):
				#print "i", i
				for j in range(b.shape[1]):
					#print "j", j
					b[i][j] = a[i+1,j+1]
					
		return b
		
	def printTypes(self, lst):
		for arg in lst:
			print type(arg)
		
class dummyOptions():
	def __init__(self):
		self.algo_method = 101
		self.algo_globalization = 112
		self.hess_approx = 131
		self.bfgs_restart = 0
		self.algo_descent = 120
		self.tol = matlabarray( [1.0e-05, 1.0e-05, 1.0e-05])
		self.dxmin = 1.0e-06
		self.miter = 500
		self.msimul = 500
		self.verbose = 7#2 set to a high level to get many messages.
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
