# -*- coding: utf-8 -*-
from __future__ import division
from runtime import *
from copy import copy

def sqplab_bfgs_(M=None,y=None,s=None,first=None,info_=None,options=None,values=None,*args,**kwargs):
    """
    #
# [M,info,values] = sqplab_bfgs (M,y,s,first,info,options,values);
#
# This procedure computes the BFGS update of the matrix M, which is
# supposed to be positive definite approximation of some Hessian. If
# y'*s is not sufficiently positive and options.algo_descent is set to
# values.powell, Powell's correction is applied to y.
#
# On entry:
#   M: matrix to update
#   y: gradient variation
#   s: point variation
#   first: if true, the procedure will initialize the matrix M to a
#       multiple of the identity before the update
#
# On return:
#   M: updated matrix

#-----------------------------------------------------------------------
#
# Author: Jean Charles Gilbert, INRIA.
#
# Copyright 2008, 2009, INRIA.
#
# SQPlab is distributed under the terms of the Q Public License version
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

#-----------------------------------------------------------------------

# Parameter

  eta = 0.2;    # constant for Powell corrections
# eta = 0.6;    # constant for Powell corrections

#-----------------------------------------------------------------------
TEST INPUT VALUES: (problem 3 options.hess_approx = 'bfgs')
sqplab_bfgs
K>> M

M =

     1     0     0
     0     1     0
     0     0     1

K>> y

y =

  -0.009783923878659
  -0.000000000000000
                   0

K>> s

s =

  -0.503054026564966
  -1.000000000000000
  -0.184478531874161

K>> first

first =

     1

K>> info

info = 

         g: [3x1 double]
        ai: []
        ae: [2x3 double]
        hl: []
     niter: 1
      flag: 0
    nsimul: [0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        ci: []
        ce: [2x1 double]
         f: 0.099563123926544
      glag: [3x1 double]
     glagn: 0.666666666736111
     feasn: 3
     compl: 0

K>> info.g

ans =

  -0.009783923878659
   1.000000000000000
                   0

K>> info.ae

ans =

   1.000000000000000   0.999999999999999   1.000000000000000
   1.000000000000000   1.999999999999999   3.000000000000000

K>> info.ce

ans =

   0.312467441560872
  -0.056489622187450

K>> info.glag

ans =

  -0.333333333013890
   0.666666666736111
  -0.333333333513889

K>> options

options = 

           algo_method: 101
    algo_globalization: 112
           hess_approx: 130
          bfgs_restart: 0
          algo_descent: 120
                   tol: [1.000000000000000e-05 1.000000000000000e-05 1.000000000000000e-05]
                 dxmin: 1.000000000000000e-06
                 miter: 500
                msimul: 500
               verbose: 2
                  fout: 1
                   inf: Inf
                   df1: 0

K>> values

values = 

                       success: 0
              fail_on_argument: 1
               fail_on_problem: 2
                 fail_on_simul: 3
                 stop_on_simul: 4
              stop_on_max_iter: 5
             stop_on_max_simul: 6
                 stop_on_dxmin: 7
          fail_on_non_decrease: 8
            fail_on_ascent_dir: 9
           fail_on_max_ls_iter: 10
              fail_on_ill_cond: 11
    stop_on_small_trust_region: 15
             fail_on_null_step: 20
         fail_on_infeasible_QP: 21
          fail_on_unbounded_QP: 22
                  fail_strange: 99
                    nsimultype: 16
                max_null_steps: 1
                        newton: 100
                  quasi_newton: 101
            cheap_quasi_newton: 102
                 unit_stepsize: 110
                    linesearch: 111
                 trust_regions: 112
                        powell: 120
                         wolfe: 121
                          bfgs: 130
                         model: 131
                         dline: '--------------------------------------------------------------------------------------'
                         eline: '======================================================================================'
                         sline: '**************************************************************************************'

TEST OUTPUT VALUES:

 M

M =

   0.205243860649550   0.003553460424288   0.000655537162145
   0.003553460424288   0.195307727402361  -0.000901167227317
   0.000655537162145  -0.000901167227317   0.200026425015353

K>> pc

pc =

   0.803070936030031

K>> info

info = 

         g: [3x1 double]
        ai: []
        ae: [2x3 double]
        hl: []
     niter: 1
      flag: 0
    nsimul: [0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        ci: []
        ce: [2x1 double]
         f: 0.099563123926544
      glag: [3x1 double]
     glagn: 0.666666666736111
     feasn: 3
     compl: 0

K>> info.g

ans =

  -0.009783923878659
   1.000000000000000
                   0

K>> info.ae

ans =

   1.000000000000000   0.999999999999999   1.000000000000000
   1.000000000000000   1.999999999999999   3.000000000000000

K>> info.ce

ans =

   0.312467441560872
  -0.056489622187450

K>> info.glag

ans =

  -0.333333333013890
   0.666666666736111
  -0.333333333513889

K>> values

values = 

                       success: 0
              fail_on_argument: 1
               fail_on_problem: 2
                 fail_on_simul: 3
                 stop_on_simul: 4
              stop_on_max_iter: 5
             stop_on_max_simul: 6
                 stop_on_dxmin: 7
          fail_on_non_decrease: 8
            fail_on_ascent_dir: 9
           fail_on_max_ls_iter: 10
              fail_on_ill_cond: 11
    stop_on_small_trust_region: 15
             fail_on_null_step: 20
         fail_on_infeasible_QP: 21
          fail_on_unbounded_QP: 22
                  fail_strange: 99
                    nsimultype: 16
                max_null_steps: 1
                        newton: 100
                  quasi_newton: 101
            cheap_quasi_newton: 102
                 unit_stepsize: 110
                    linesearch: 111
                 trust_regions: 112
                        powell: 120
                         wolfe: 121
                          bfgs: 130
                         model: 131
                         dline: '--------------------------------------------------------------------------------------'
                         eline: '======================================================================================'
                         sline: '**************************************************************************************'
    """
#    varargin = cellarray(args)
#    nargin = 7-[M,y,s,first,info,options,values].count(None)+len(args)

    info=copy(info_)
 # Parameter
   # eta = 0.6;    # constant for Powell corrections
    eta=0.2 # constant for Powell corrections
# initalization
    n=length_(s)
    pc=1
    info.flag=values.success
# prepare the update of M

    if norm_(s) == 0:
        info.flag=values.fail_strange
        if options.verbose >= 3:
            fprintf_(options.fout,'\n### sqplab_bfgs: null step s\n\n')
        return M,pc,info,values
    ys=y.T .dot(s)
    if options.verbose >= 4:
        fprintf_(options.fout," y'*s/(s'*s) = %9.3e\n"%(ys / (s.T.dot(s))))
    Ms=M.dot(s)
    sMs=s.T.dot(Ms)
    if sMs <= 0:
        info.flag=values.fail_strange
        if options.verbose >= 3:
            fprintf_(options.fout,'\n### sqplab_bfgs: BFGS Hessian approximation is not positive definite:\n')
            fprintf_(options.fout,"            s'*M*s = %g <= 0\n\n"%(sMs))
        return M,pc,info,values
#  disp_(['sqplab_bfgs: curvature condition fails:', str(ys < eta*sMs)])

    if (options.algo_descent == values.powell) and (ys < eta * sMs): # Powell's correction
        pc=(1 - eta) * sMs / (sMs - ys)
        if options.verbose >= 4:
            fprintf_(options.fout,"  Powell's corrector = %7.1e\n"%(pc))
        y=pc * y + (1 - pc) * Ms
        ys=y.T.dot(s)  # update ys, since y has changed
        if options.verbose >= 4:
            fprintf_(options.fout," (new y'*s/(s'*s) = %7.1e\n)"%(ys / (s.T.dot(s))))
        if ys <= 0:
            info.flag=values.fail_strange
            if options.verbose >= 4:
                fprintf_(options.fout,"\n### sqplab_bfgs: y'*s = %9.3e not positive despite correction:\n\n"%(ys).T)
            return M,pc,info,values
    else:
        if ys <= 0:
            if options.verbose >= 4:
                fprintf_(options.fout,"\n### sqplab_bfgs: y'*s = %9.3e is nonpositive\n\n"%(ys).T)
            info.flag=values.fail_strange
            return M,pc,info,values
    # case when the matrix has to be initialized

    if first:
     #   ol = ys/(s.T.dot(s));
        ol=(y.T.dot( y)) / ys
        M=ol * eye_(n)
        if options.verbose >= 4:
            fprintf_(options.fout,'  OL coefficient = %g\n'%(ol))
        Ms=ol*s
        sMs=s.T .dot(Ms)
    #   if options.verbose >= 6:
    #     values.file_eigM = fopen_('sqplab_eigM.m','w');
    #     fprintf_(values.file_eigM,'  Powell''s corrector    Eigenvalues of the updated matrix M\nn=%0i\n',n);
# update the matrix M (TODO: update with the formula used in SQPpro, which is more stable)

    M=M - (Ms.dot(Ms.T)) / sMs + (y .dot(y.T)) / ys   # BFGS formula
    if options.verbose >= 6:
        eigM=sort_(eig_(M))
        fprintf_(options.fout,'  eig(M): min = %g, max = %g, cond = %g\n'%(min_(eigM),max_(eigM),max_(eigM) / min_(eigM)))
    #   fprintf_(values.file_eigM,'%22.15e',pc);
    #   for i in range(0,n):
    #     fprintf_(values.file_eigM,'  %22.15e',eigM(i));
    #   fprintf_(values.file_eigM,'\n');
    return M,pc,info,values