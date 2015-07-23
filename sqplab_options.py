# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 10:42:12 2014
%%
% [info,options,values] = sqplab_options (info,options)
%
% Set the options of the optimization solver 'ecdfo'

%-----------------------------------------------------------------------
%
% Author: Jean Charles Gilbert, INRIA.
%
% Copyright 2008, 2009, INRIA.
%
% ecdfo is distributed under the terms of the Q Public License version
% 1.0.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the Q Public
% License version 1.0 for more details.
%
% You should have received a copy of the Q Public License version 1.0
% along with this program.  If not, see
% <http://doc.trolltech.com/3.0/license.html>.
%
%-----------------------------------------------------------------------


@author: jaco_da
"""

# Autogenerated with SMOP version 
# c:\Users\jaco_da\AppData\Local\Continuum\Anaconda\Scripts\smop-script.py sqplab_options.m

from __future__ import division
#try:
from runtime import *
from numpy import inf
#except ImportError:
    #from smop.runtime import *

def sqplab_options_(info=None,options=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 2-[info,options].count(None)+len(args)

    values.success=0
    values.fail_on_argument=1
    values.fail_on_problem=2
    values.fail_on_simul=3
    values.stop_on_simul=4
    values.stop_on_max_iter=5
    values.stop_on_max_simul=6
    values.stop_on_dxmin=7
    values.fail_on_non_decrease=8
    values.fail_on_ascent_dir=9
    values.fail_on_max_ls_iter=10
    values.fail_on_ill_cond=11
    values.stop_on_small_trust_region=15
    values.fail_on_null_step=20
    values.fail_on_infeasible_QP=21
    values.fail_on_unbounded_QP=22
    values.fail_strange=99
    values.nsimultype=16
    values.max_null_steps=1
    values.newton=100
    values.quasi_newton=101
    values.cheap_quasi_newton=102
    values.unit_stepsize=110
    values.linesearch=111
    values.trust_regions=112
    values.powell=120
    values.wolfe=121
    values.bfgs=130
    values.model=131
    info.flag=values.success
    if isempty_(options):
        options.xxx=0
    if isfield_(options,char('fout')):
        if options.fout < 0:
            fprintf_(char('\\n### ecdfo: options.fout = "%0i" is not a valid file identifier (use \'fopen\' to have a valid one)'),options.fout)
            fprintf_(char('\\n            options.fout is set to 1\\n\\n'))
            options.fout=1
    else:
        options.fout=1
    if isfield_(options,char('verbose')):
        if (options.verbose < 0) or (options.verbose > 6):
            fprintf_(options.fout,char('\\n### ecdfo: options.verbose = "%0i" and should be in [0,6], reset to 1\\n\\n'),options.verbose)
            options.verbose=1
    else:
        options.verbose=1
    if isfield_(options,char('algo_method')):
        if char('newton') == lower___(regexprep___(strtrim___(options.algo_method),char('  *'),char(' '))):
            options.algo_method=values.newton
        else:
            if [char('quasi-newton'),char('quasi newton'),char('quasinewton')] == lower___(regexprep___(strtrim___(options.algo_method),char('  *'),char(' '))):
                options.algo_method=values.quasi_newton
            else:
                if [char('cheap quasi-newton'),char('cheap quasi newton'),char('cheap quasinewton')] == lower___(regexprep___(strtrim___(options.algo_method),char('  *'),char(' '))):
                    options.algo_method=values.cheap_quasi_newton
                else:
                    if options.verbose:
                        fprintf_(options.fout,char('\\n### ecdfo: options.algo_method "%s" not recognized\\n\\n'),options.algo_method)
                    info.flag=values.fail_on_argument
                    return info,options,values
    else:
        options.algo_method=values.quasi_newton
    if isfield_(options,char('algo_globalization')):
        if [char('unit step-size'),char('unit stepsize')] == lower___(regexprep___(strtrim___(options.algo_globalization),char('  *'),char(' '))):
            options.algo_globalization=values.unit_stepsize
        else:
            if [char('line-search'),char('linesearch')] == lower___(regexprep___(strtrim___(options.algo_globalization),char('  *'),char(' '))):
                options.algo_globalization=values.linesearch
            else:
                if [char('trust regions'),char('trust-regions'),char('trustregions')] == lower___(regexprep___(strtrim___(options.algo_globalization),char('  *'),char(' '))):
                    options.algo_globalization=values.trust_regions
                else:
                    if options.verbose:
                        fprintf_(options.fout,char('\\n### ecdfo: options.algo_globalization "%s" not recognized\\n\\n'),options.algo_globalization)
                    info.flag=values.fail_on_argument
                    return info,options,values
    else:
        options.algo_globalization=values.linesearch
    if isfield_(options,char('algo_descent')):
        if char('powell') == lower__(regexprep__(strtrim__(options.algo_descent),char('  *'),char(' '))):
            options.algo_descent=values.powell
        else:
            if char('wolfe') == lower__(regexprep__(strtrim__(options.algo_descent),char('  *'),char(' '))):
                options.algo_descent=values.wolfe
            else:
                if options.verbose:
                    fprintf_(options.fout,char('\\n### ecdfo: options.algo_descent "%s" not recognized\\n\\n'),options.algo_descent)
                info.flag=values.fail_on_argument
                return info,options,values
    if isfield_(options,char('dxmin')):
        if (options.dxmin <= 0):
            if options.verbose:
                fprintf_(options.fout,char('\\n### ecdfo: options.dxmin = %g must be > 0\\n\\n'),options.dxmin)
            info.flag=values.fail_on_argument
            return info,options,values
    else:
        options.dxmin=1e-08
    if isfield_(options,char('inf')):
        if options.inf <= 0:
            if options.verbose:
                fprintf_(char('\\n### ecdfo: incorrect value of options.inf %g (should be > 0)\\n\\n'),options.inf)
            info.flag=values.fail_on_argument
            return info,options,values
    else:
        options.inf=inf
    if isfield_(options,char('miter')):
        if options.miter <= 0:
            if options.verbose:
                fprintf_(char('\\n### ecdfo: incorrect value of options.miter %g (should be > 0)\\n\\n'),options.miter)
            info.flag=values.fail_on_argument
            return info,options,values
    else:
        options.miter=1000
    if isfield_(options,char('tol')):
        if any_(options.tol <= 0):
            if options.verbose:
                fprintf_(char('\\n### ecdfo: incorrect value of some options.tol (should be > 0)\\n\\n'))
            info.flag=values.fail_on_argument
            return info,options,values
    else:
        options.tol=[[1e-06],[1e-06],[1e-06]]
    if not isfield_(options,char('df1')):
        options.df1=0
    if isfield_(options,char('hess_approx')):
        if char('bfgs') == lower__(regexprep__(strtrim__(options.hess_approx),char('  *'),char(' '))):
            options.hess_approx=values.bfgs
        else:
            if char('model') == lower__(regexprep__(strtrim__(options.hess_approx),char('  *'),char(' '))):
                options.hess_approx=values.model
            else:
                if options.verbose:
                    fprintf_(options.fout,char('\\n### ecdfo: options.hess_approx "%s" not recognized\\n\\n'),options.hess_approx)
                info.flag=values.fail_on_argument
                return info,options,values
    return info,options,values

    
#class ecdfoValues():
#	def __init__(self):
#		self.success=0
#		self.fail_on_argument=1
#		self.fail_on_problem=2
#		self.fail_on_simul=3
#		self.stop_on_simul=4
#		self.stop_on_max_iter=5
#		self.stop_on_max_simul=6
#		self.stop_on_dxmin=7
#		self.fail_on_non_decrease=8
#		self.fail_on_ascent_dir=9
#		self.fail_on_max_ls_iter=10
#		self.fail_on_ill_cond=11
#		self.stop_on_small_trust_region=15
#		self.fail_on_null_step=20
#		self.fail_on_infeasible_QP=21
#		self.fail_on_unbounded_QP=22
#		self.fail_strange=99
#		self.nsimultype=16
#		self.max_null_steps=1
#		self.newton=100
#		self.quasi_newton=101
#		self.cheap_quasi_newton=102
#		self.unit_stepsize=110
#		self.linesearch=111
#		self.trust_regions=112
#		self.powell=120
#		self.wolfe=121
#		self.bfgs=130
#		self.model=131

#def sqplab_options_(info=None,options=None,*args,**kwargs):
#    #varargin = cellarray(args)
#    #nargin = 2-[info,options].count(None)+len(args)
#
#    values = ecdfoValues()
#    
#    info.flag=values.success
#    if isempty_(options):
#        options.xxx=0
#    if isfield_(options,char('fout')):
#        if options.fout < 0:
#            fprintf_(char('\\n### ecdfo: options.fout = "%0i" is not a valid file identifier (use \'fopen\' to have a valid one)'),options.fout)
#            fprintf_(char('\\n            options.fout is set to 1\\n\\n'))
#            options.fout=1
#    else:
#        options.fout=1
#    if isfield_(options,char('verbose')):
#        if (options.verbose < 0) or (options.verbose > 6):
#            fprintf_(options.fout,char('\\n### ecdfo: options.verbose = "%0i" and should be in [0,6], reset to 1\\n\\n'),options.verbose)
#            options.verbose=1
#    else:
#        options.verbose=1
#    
#    if isfield_(options,char('algo_method')):
#        #print "type char('newton')", type(char('newton'))
#        #print "type lower___(...)", type(lower___(regexprep___(strtrim___(options.algo_method),char('  *'),char(' '))))
#        if char('newton') == lower___(regexprep___(strtrim___(options.algo_method),char('  *'),char(' '))):
#            options.algo_method=values.newton
#        else:
#            if [char('quasi-newton'),char('quasi newton'),char('quasinewton')] == lower___(regexprep___(strtrim___(options.algo_method),char('  *'),char(' '))):
#                options.algo_method=values.quasi_newton
#            else:
#                if [char('cheap quasi-newton'),char('cheap quasi newton'),char('cheap quasinewton')] == lower___(regexprep___(strtrim___(options.algo_method),char('  *'),char(' '))):
#                    options.algo_method=values.cheap_quasi_newton
#                else:
#                    if options.verbose:
#                        fprintf_(options.fout,char('\\n### ecdfo: options.algo_method "%s" not recognized\\n\\n'),options.algo_method)
#                    info.flag=values.fail_on_argument
#                    #print "return 1"																				
#                    return info,options,values
#    else:
#        options.algo_method=values.quasi_newton
#    if isfield_(options,char('algo_globalization')):
#        if [char('unit step-size'),char('unit stepsize')] == lower___(regexprep___(strtrim___(options.algo_globalization),char('  *'),char(' '))):
#            options.algo_globalization=values.unit_stepsize
#        else:
#            if [char('line-search'),char('linesearch')] == lower___(regexprep___(strtrim___(options.algo_globalization),char('  *'),char(' '))):
#                options.algo_globalization=values.linesearch
#            else:
#                if [char('trust regions'),char('trust-regions'),char('trustregions')] == lower___(regexprep___(strtrim___(options.algo_globalization),char('  *'),char(' '))):
#                    options.algo_globalization=values.trust_regions
#                else:
#                    if options.verbose:
#                        fprintf_(options.fout,char('\\n### ecdfo: options.algo_globalization "%s" not recognized\\n\\n'),options.algo_globalization)
#                    info.flag=values.fail_on_argument
#                    #print "return 2"																				
#                    return info,options,values
#    else:
#        options.algo_globalization=values.linesearch
#    if isfield_(options,char('algo_descent')):
#        if char('powell') == lower__(regexprep__(strtrim__(options.algo_descent),char('  *'),char(' '))):
#            options.algo_descent=values.powell
#        else:
#            if char('wolfe') == lower__(regexprep__(strtrim__(options.algo_descent),char('  *'),char(' '))):
#                options.algo_descent=values.wolfe
#            else:
#                if options.verbose:
#                    fprintf_(options.fout,char('\\n### ecdfo: options.algo_descent "%s" not recognized\\n\\n'),options.algo_descent)
#                info.flag=values.fail_on_argument
#                #print "return 3"																
#                return info,options,values
#    if isfield_(options,char('dxmin')):
#        if (options.dxmin <= 0):
#            if options.verbose:
#                fprintf_(options.fout,char('\\n### ecdfo: options.dxmin = %g must be > 0\\n\\n'),options.dxmin)
#            info.flag=values.fail_on_argument
#            #print "return 4"												
#            return info,options,values
#    else:
#        options.dxmin=1e-08
#    if isfield_(options,char('inf')):
#        if options.inf <= 0:
#            if options.verbose:
#                fprintf_(char('\\n### ecdfo: incorrect value of options.inf %g (should be > 0)\\n\\n'),options.inf)
#            info.flag=values.fail_on_argument
#            #print "return 5"
#            return info,options,values
#    else:
#        options.inf=inf
#    if isfield_(options,char('miter')):
#        if options.miter <= 0:
#            if options.verbose:
#                fprintf_(char('\\n### ecdfo: incorrect value of options.miter %g (should be > 0)\\n\\n'),options.miter)
#            info.flag=values.fail_on_argument
#            #print "return 6"												
#            return info,options,values
#    else:
#        options.miter=1000
#    #print "isfield optrions tol", isfield_(options,char('tol'))								
#    if isfield_(options,char('tol')):
#        #print "hellooo any"					
#        if any_(options.tol <= 0):
#            if options.verbose:
#                fprintf_(char('\\n### ecdfo: incorrect value of some options.tol (should be > 0)\\n\\n'))
#            info.flag=values.fail_on_argument
#            #print "return 7"																								
#            return info,options,values
#    else:
#        options.tol=[[1e-06],[1e-06],[1e-06]]
#    if not isfield_(options,char('df1')):
#        options.df1=0
#    if isfield_(options,char('hess_approx')):
#        if char('bfgs') == lower__(regexprep__(strtrim__(options.hess_approx),char('  *'),char(' '))):
#            options.hess_approx=values.bfgs
#        else:
#            if char('model') == lower__(regexprep__(strtrim__(options.hess_approx),char('  *'),char(' '))):
#                options.hess_approx=values.model
#            else:
#                if options.verbose:
#                    fprintf_(options.fout,char('\\n### ecdfo: options.hess_approx "%s" not recognized\\n\\n'),options.hess_approx)
#                info.flag=values.fail_on_argument
#                #print "return 8"																												
#                return info,options,values
#    return info,options,values
