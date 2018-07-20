# -*- coding: utf-8 -*-

from runtime import *
from numpy import inf, array
from copy import copy
from ecdfo_values import ecdfoValues


def ecdfo_options_(info_=None,options_=None,*args,**kwargs):
    """
# [info,options,values] = ecdfo_options (info,options)
#
# Set the default options of the optimizer 'ecdfo'
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
            fprintf_('\n### ecdfo_options: options.fout = "%0i" is not a valid file identifier' % (options.fout))
            fprintf_('\n    options.fout is set to 1\n\n')
            options.fout = 1
    else:
        options.fout = 1  # default value

    if isfield_(options, 'verbose'):
        if (options.verbose < 0) or (options.verbose > 6):
            fprintf_(options.fout,
                     '\n### ecdfo_options: options.verbose = "%0i" and should be in [0,6], reset to 1\n\n' % (options.verbose))
            options.verbose = 1
    else:
        options.verbose = 1  # default value

    if isfield_(options, 'algo_descent'):
        algo_descent_str = options.algo_descent.lower().replace(" ", "")
        if (algo_descent_str == 'powell'):
            options.algo_descent = values.powell
        elif (algo_descent_str == 'wolfe'):
            options.algo_descent = values.wolfe
        else:
            if options.verbose:
                fprintf_(options.fout,
                         '\n### ecdfo_options: options.algo_descent "%s" not recognized\n\n' % (options.algo_descent))
            info.flag = values.fail_on_argument
            return info, options, values
    else:
        options.algo_descent = values.powell  # default value

    if isfield_(options, 'dxmin'):
        if (options.dxmin <= 0):
            if options.verbose:
                fprintf_(options.fout, '\n### ecdfo_options: options.dxmin = %g must be > 0\n\n' % (options.dxmin))
            info.flag = values.fail_on_argument
            return info, options, values
    else:
        options.dxmin = 1e-08  # default value

    if isfield_(options, 'inf'):
        if options.inf <= 0:
            if options.verbose:
                fprintf_('\n### ecdfo_options: incorrect value of options.inf %g (should be > 0)\n\n' % (options.inf))
            info.flag = values.fail_on_argument
            return info, options, values
    else:
        options.inf = inf  # default value

    if isfield_(options, 'miter'):
        if options.miter <= 0:
            if options.verbose:
                fprintf_('\n### ecdfo_options: incorrect value of options.miter %g (should be > 0)\n\n' % (options.miter))
            info.flag = values.fail_on_argument
            return info, options, values
    else:
        options.miter = 1000  # default value

    if isfield_(options,'msimul'):
        if options.msimul <= 0:
           if options.verbose:
              fprintf_('\n### ecdfo_options: incorrect value of options.msimul %g (should be > 0)\n\n',options.msimul)
           info.flag = values.fail_on_argument;
           return info, options, values
    else:
        options.msimul = 1000  # default value

    if isfield_(options, 'tol_grad'):
        if (options.tol_grad < 0):
            if options.verbose:
                fprintf_('\n### ecdfo_options: incorrect value of some options.tol_grad (should be > 0)\n\n')
            info.flag = values.fail_on_argument
            return info, options, values
        elif options.tol_grad==0:
                options.tol_grad=1e-5
    else:
        options.tol_grad = 1e-5  # default value

    if isfield_(options, 'tol_feas'):
        if (options.tol_feas < 0):
            if options.verbose:
                fprintf_('\n### ecdfo_options: incorrect value of some options.tol_feas (should be > 0)\n\n')
            info.flag = values.fail_on_argument
            return info, options, values
        elif options.tol_feas==0:
                options.tol_feas=1e-5
    else:
        options.tol_feas = 1e-5  # default value

    if isfield_(options, 'tol_bnds'):
        if (options.tol_bnds < 0):
            if options.verbose:
                fprintf_('\n### ecdfo_options: incorrect value of some options.tol_bnds (should be > 0)\n\n')
            info.flag = values.fail_on_argument
            return info, options, values
        elif options.tol_bnds==0:
                options.tol_bnds=1e-5
    else:
        options.tol_bnds = 1e-5  # default value

    if isfield_(options, 'hess_approx'):
        if 'bfgs' == options.hess_approx.lower().replace(" ", ""):
            options.hess_approx = values.bfgs
        elif 'model' == options.hess_approx.lower().replace(" ", ""):
            options.hess_approx = values.model
        else:
            if options.verbose:
                fprintf_(options.fout,
                         '\n### ecdfo_options: options.hess_approx "%s" not recognized' % (options.hess_approx))
                fprintf_('\n    options.hess_approx must be either "model" or "bfgs"\n\n')
            info.flag = values.fail_on_argument
            return info, options, values
    else:
        options.hess_approx = values.model  # default value

    if isfield_(options, 'whichmodel'):
        if 'subbasis' == options.whichmodel.lower().replace(" ", ""):
            options.whichmodel = values.subbasis
        elif 'frobnorm' == options.whichmodel.lower().replace(" ", ""):
            options.whichmodel = values.frobnorm
        elif 'l2norm' == options.whichmodel.lower().replace(" ", ""):
            options.whichmodel = values.l2norm
        elif 'regression' == options.whichmodel.lower().replace(" ", ""):
            options.whichmodel = values.regression
        else:
            if options.verbose:
                fprintf_(options.fout,
                         '\n### ecdfo_options: options.whichmodel "%s" not recognized' % (options.whichmodel))
                fprintf_('\n    options.whichmodel must be either "Subbasis", "Frobnorm", "L2norm" or "Regression"\n\n')
            info.flag = values.fail_on_argument
            return info, options, values
    else:
        options.whichmodel = 0  # default value

    if isfield_(options, 'final_degree'):
        if 'linear' == options.final_degree.lower().replace(" ", ""):
            options.final_degree = values.linear
        elif 'diagonal' == options.final_degree.lower().replace(" ", ""):
            options.final_degree = values.diagonal
        elif 'quadratic' == options.final_degree.lower().replace(" ", ""):
            options.final_degree = values.quadratic
        else:
            if options.verbose:
                fprintf_(options.fout,
                         '\n### ecdfo_options: options.final_degree "%s" not recognized' % (options.final_degree))
                fprintf_('\n    options.final_degree must be either "linear", "diagonal" or "quadratic"\n\n')
            info.flag = values.fail_on_argument
            return info, options, values
    else:
        options.final_degree = values.diagonal  # default value

    return info, options, values
