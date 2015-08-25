# -*- coding: utf-8 -*-
from __future__ import division
#try:
from runtime import *
from copy import copy
#except ImportError:
#    from smop.runtime import *
def sqplab_checkoptions_(nb=None,mi=None,me=None,ms=None,info_=None,options_=None,values=None,*args,**kwargs):
    """
% [info] = sqplab_checkoptions (nb,mi,me,ms,info,options,values);
%
% Check the structure 'options' to see whether the required options are
% compatible with the solver capabilities.

%-----------------------------------------------------------------------
%
% Author: Jean Charles Gilbert, INRIA.
%
% Copyright 2008, 2009, INRIA.
%
% SQPlab is distributed under the terms of the Q Public License version
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
    """
#    varargin = cellarray(args)
#    nargin = 7-[nb,mi,me,ms,info,options,values].count(None)+len(args)

    options=copy(options_)
    info=copy(info_)

    info.flag=values.success
    if options.algo_method == values.cheap_quasi_newton:
        if options.verbose:
            fprintf_(options.fout,'\n### sqplab_checkoptions: the cheap quasi-Newton method is not a valid\n')
            fprintf_(options.fout,'### approach when there is no state constraint\n\n')
        info.flag=values.fail_on_problem
        return info,options
    if values.unit_stepsize == options.algo_globalization:
        if options.algo_method != values.newton:
            if not isfield_(options,'algo_descent'):
                if options.verbose:
                    fprintf_(options.fout,'\n### sqplab_checkoptions: positive definiteness of the matrices is ensured\n')
                    fprintf_(options.fout,'### by Powell corrections\n\n')
                options.algo_descent=values.powell
            else:
                if options.algo_descent == values.wolfe:
                    if options.verbose:
                        fprintf_(options.fout,'\n### sqplab_checkoptions: positive definiteness of the matrices cannot be ensured\n')
                        fprintf_(options.fout,'### by the Wolfe linesearch when unit stepsize is required; Powell corrections\n')
                        fprintf_(options.fout,'### will be used instead\n\n')
                    options.algo_descent=values.powell
    else:
        if values.linesearch == options.algo_globalization:
            if options.algo_method == values.newton:
                if isfield_(options,'algo_descent'):
                    if options.verbose:
                        fprintf_(options.fout,"\n### sqplab_checkoptions: descent cannot be ensured for Newton's method\n")
                        fprintf_(options.fout,'### by using Powell corrections or Wolfe linesearch\n\n')
                    info.flag=values.fail_on_argument
                    return info,options
                else:
                    if options.verbose:
                        fprintf_(options.fout,"\n### sqplab_checkoptions: Armijo's linesearch can fail with Newton's method,\n")
                        fprintf_(options.fout,'###                      use unit step-size instead\n\n')
            else:
                if not isfield_(options,'algo_descent'):
                    if options.verbose:
                        fprintf_(options.fout,'\n### sqplab_checkoptions: descent ensured by Powell corrections\n')
                        if nb + mi + me + ms == 0:
                            fprintf_(options.fout,'### setting "options.algo_descent = \'Wolfe\'" should be better\n\n')
                    options.algo_descent=values.powell
                else:
                    if (options.algo_descent == values.wolfe) and (nb + mi + me + ms != 0):
                        if options.verbose:
                            fprintf_(options.fout,'\n### sqplab_checkoptions: positive definiteness of the matrices cannot be ensured\n')
                            fprintf_(options.fout,'### by the Wolfe linesearch when constraints are present; Powell corrections\n')
                            fprintf_(options.fout,'### will be used instead\n\n')
                        options.algo_descent=values.powell
    return info,options