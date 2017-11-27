# -*- coding: utf-8 -*-

from runtime import *
from numpy import inf,array
from copy import copy


class ecdfoValues():
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

def sqplab_options_(info_=None,options_=None,*args,**kwargs):
    """
# [info,options,values] = sqplab_options (info,options)
#
# Set the options of the optimization solver 'ecdfo'

#-----------------------------------------------------------------------
#
# Author: Jean Charles Gilbert, INRIA.
#
# Copyright 2008, 2009, INRIA.
#
# ecdfo is distributed under the terms of the Q Public License version
# 1.0.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the Q Public
# License version 1.0 for more details.
#
# You should have received a copy of the Q Public License version 1.0
# along with this program.  If not, see
# <http://doc.trolltech.com/3.0/license.html>.
#
#-----------------------------------------------------------------------
    """
#    varargin = cellarray(args)
#    nargin = 2-[info,options].count(None)+len(args)

    info=copy(info_)
    options=copy(options_)
# Define ecdfo constant values

    values = ecdfoValues()
# Initialization

    info.flag=values.success
# Force the existence of 'options'

    if isempty_(options):
        options.xxx=0
# Set options
    
    if isfield_(options,'fout'):
        if options.fout < 0:
            fprintf_('\n### ecdfo: options.fout = "%0i" is not a valid file identifier (use \'fopen\' to have a valid one)'%(options.fout))
            fprintf_('\n            options.fout is set to 1\n\n')
            options.fout=1
    else:
        options.fout=1
    if isfield_(options,'verbose'):
        if (options.verbose < 0) or (options.verbose > 6):
            fprintf_(options.fout,'\n### ecdfo: options.verbose = "%0i" and should be in [0,6], reset to 1\n\n'%(options.verbose))
            options.verbose=1
    else:
        options.verbose=1
    if isfield_(options,'algo_method'):
        algo_method_str=options.algo_method.lower().replace(" ","")
        if 'newton' == algo_method_str:
            options.algo_method=values.newton
        elif (algo_method_str=='quasi-newton') or (algo_method_str=='quasi newton') or (algo_method_str=='quasinewton'): 
            options.algo_method=values.quasi_newton
        elif (algo_method_str=='cheap quasi-newton')or(algo_method_str=='cheap quasi newton') or(algo_method_str=='cheap quasinewton'):
            options.algo_method=values.cheap_quasi_newton
        else:
            if options.verbose:
                fprintf_(options.fout,'\n### ecdfo: options.algo_method "{}" not recognized\n\n'.format(options.algo_method))
            info.flag=values.fail_on_argument
            return info,options,values
    else:
        options.algo_method=values.quasi_newton
    if isfield_(options,'algo_globalization'):
        algo_globalization_str=options.algo_globalization.lower().replace(" ","")
        if (algo_globalization_str=='unit step-size') or (algo_globalization_str=='unit stepsize') :
            options.algo_globalization=values.unit_stepsize
        elif (algo_globalization_str=='line-search') or (algo_globalization_str=='linesearch'):
            options.algo_globalization=values.linesearch
        elif (algo_globalization_str=='trust regions')or(algo_globalization_str=='trust-regions')or(algo_globalization_str=='trustregions'):
            options.algo_globalization=values.trust_regions
        else:
            if options.verbose:
                fprintf_(options.fout,'\n### ecdfo: options.algo_globalization "%s" not recognized\n\n'%(options.algo_globalization))
            info.flag=values.fail_on_argument
            return info,options,values
    else:
        options.algo_globalization=values.linesearch
    if isfield_(options,'algo_descent'):
        algo_descent_str=options.algo_descent.lower().replace(" ","")
        if (algo_descent_str=='powell'):
            options.algo_descent=values.powell
        elif (algo_descent_str=='wolfe'):
            options.algo_descent=values.wolfe
        else:
            if options.verbose:
                fprintf_(options.fout,'\n### ecdfo: options.algo_descent "%s" not recognized\n\n'%(options.algo_descent))
            info.flag=values.fail_on_argument
            return info,options,values
    if isfield_(options,'dxmin'):
        if (options.dxmin <= 0):
            if options.verbose:
                fprintf_(options.fout,'\n### ecdfo: options.dxmin = %g must be > 0\n\n'%(options.dxmin))
            info.flag=values.fail_on_argument
            return info,options,values
    else:
        options.dxmin=1e-08
    if isfield_(options,'inf'):
        if options.inf <= 0:
            if options.verbose:
                fprintf_('\n### ecdfo: incorrect value of options.inf %g (should be > 0)\n\n'%(options.inf))
            info.flag=values.fail_on_argument
            return info,options,values
    else:
        options.inf=inf
    if isfield_(options,'miter'):
        if options.miter <= 0:
            if options.verbose:
                fprintf_('\n### ecdfo: incorrect value of options.miter %g (should be > 0)\n\n'%(options.miter))
            info.flag=values.fail_on_argument
            return info,options,values
    else:
        options.miter=1000

#    if isfield_(options,'msimul'):
#       if options.msimul <= 0:
#         if options.verbose:
#            fprintf_('\n### ecdfo: incorrect value of options.msimul %g (should be > 0)\n\n',options.msimul)
#         info.flag = values.fail_on_argument;
#         return
#    else:
#       options.msimul = 1000;
    
    if isfield_(options,'tol'):
        if any_(options.tol <= 0):
            if options.verbose:
                fprintf_('\n### ecdfo: incorrect value of some options.tol (should be > 0)\n\n')
            info.flag=values.fail_on_argument
            return info,options,values
    else:
        options.tol=array([[1e-06],[1e-06],[1e-06]])
    if not isfield_(options,'df1'):
        options.df1=0
    if isfield_(options,'hess_approx'):
        if 'bfgs' == options.hess_approx.lower().replace(" ",""):
            options.hess_approx=values.bfgs
        elif 'model' == options.hess_approx.lower().replace(" ",""):
            options.hess_approx=values.model
        else:
            if options.verbose:
                fprintf_(options.fout,'\n### ecdfo: options.hess_approx "%s" not recognized\n\n'%(options.hess_approx))
            info.flag=values.fail_on_argument
            return info,options,values
    return info,options,values
