# -*- coding: utf-8 -*-

from runtime import *
from numpy import inf, array
from copy import copy


class ecdfoValues():
    def __init__(self):
        # Define ecdfo constant values

        self.success = 0;  # solution found
        self.fail_on_argument = 1;  # an argument is wrong
        self.fail_on_problem = 2;  # unaccepted problem structure
        self.fail_on_simul = 3;  # error on a simulation
        self.stop_on_simul = 4; # simulator wants to stop
        self.stop_on_max_iter = 5;  # max iterations
        self.stop_on_max_simul = 6;  # max simulations
        self.stop_on_dxmin = 7;  # stop on dxmin
        self.fail_on_non_decrease = 8;  # the merit function no longer decrease
        self.fail_on_ascent_dir = 9;  # nondescent direction in linesearch
        self.fail_on_ill_cond = 11;  # ill-conditioning
        self.stop_on_small_trust_region = 15; # small trust-region radius
        self.fail_on_null_step = 20;  # null step d is solution of 'values.max_null_steps' QPs
        self.fail_on_infeasible_QP = 21;  # infeasible QP
        self.fail_on_unbounded_QP = 22;  # unbounded QP
        self.fail_unexpected = 30;  # should not have happened

        self.nsimultype = 4;  # nb of simulation types
        self.max_null_steps = 1;  # maximum nbr of null QP steps

        self.powell = 120
        self.wolfe = 121
        self.bfgs = 130
        self.model = 131

def ecdfo_options_(info_=None,options_=None,*args,**kwargs):
    """
# [info,options,values] = ecdfo_options (info,options)
#
# Set the options of the optimizer 'ecdfo'
    """

    info = copy(info_)
    options = copy(options_)

    # Define ecdfo constant values

    values = ecdfoValues()

    # Initialization

    info.flag = values.success

    # Force the existence of 'options'

    if isempty_(options):
        options.xyz = 0

    # Set options

    if isfield_(options, 'fout'):
        if options.fout < 0:
            fprintf_(
                '\n### ecdfo: options.fout = "%0i" is not a valid file identifier (use \'fopen\' to have a valid one)' % (
                    options.fout))
            fprintf_('\n            options.fout is set to 1\n\n')
            options.fout = 1
    else:
        options.fout = 1
    if isfield_(options, 'verbose'):
        if (options.verbose < 0) or (options.verbose > 6):
            fprintf_(options.fout,
                     '\n### ecdfo: options.verbose = "%0i" and should be in [0,6], reset to 1\n\n' % (options.verbose))
            options.verbose = 1
    else:
        options.verbose = 1
    if isfield_(options, 'algo_descent'):
        algo_descent_str = options.algo_descent.lower().replace(" ", "")
        if (algo_descent_str == 'powell'):
            options.algo_descent = values.powell
        elif (algo_descent_str == 'wolfe'):
            options.algo_descent = values.wolfe
        else:
            if options.verbose:
                fprintf_(options.fout,
                         '\n### ecdfo: options.algo_descent "%s" not recognized\n\n' % (options.algo_descent))
            info.flag = values.fail_on_argument
            return info, options, values
    if isfield_(options, 'dxmin'):
        if (options.dxmin <= 0):
            if options.verbose:
                fprintf_(options.fout, '\n### ecdfo: options.dxmin = %g must be > 0\n\n' % (options.dxmin))
            info.flag = values.fail_on_argument
            return info, options, values
    else:
        options.dxmin = 1e-08
    if isfield_(options, 'inf'):
        if options.inf <= 0:
            if options.verbose:
                fprintf_('\n### ecdfo: incorrect value of options.inf %g (should be > 0)\n\n' % (options.inf))
            info.flag = values.fail_on_argument
            return info, options, values
    else:
        options.inf = inf
    if isfield_(options, 'miter'):
        if options.miter <= 0:
            if options.verbose:
                fprintf_('\n### ecdfo: incorrect value of options.miter %g (should be > 0)\n\n' % (options.miter))
            info.flag = values.fail_on_argument
            return info, options, values
    else:
        options.miter = 1000
    if isfield_(options,'msimul'):
        if options.msimul <= 0:
           if options.verbose:
              fprintf_('\n### ecdfo: incorrect value of options.msimul %g (should be > 0)\n\n',options.msimul)
           info.flag = values.fail_on_argument;
           return info, options, values
    else:
        options.msimul = 1000;
    if isfield_(options, 'tol'):
        if any_(options.tol <= 0):
            if options.verbose:
                fprintf_('\n### ecdfo: incorrect value of some options.tol (should be > 0)\n\n')
            info.flag = values.fail_on_argument
            return info, options, values
    else:
        options.tol = array([[1e-06], [1e-06], [1e-06]])
    if isfield_(options, 'hess_approx'):
        if 'bfgs' == options.hess_approx.lower().replace(" ", ""):
            options.hess_approx = values.bfgs
        elif 'model' == options.hess_approx.lower().replace(" ", ""):
            options.hess_approx = values.model
        else:
            if options.verbose:
                fprintf_(options.fout,
                         '\n### ecdfo: options.hess_approx "%s" not recognized\n\n' % (options.hess_approx))
            info.flag = values.fail_on_argument
            return info, options, values

    return info, options, values
