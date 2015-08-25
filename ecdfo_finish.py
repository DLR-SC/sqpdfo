# -*- coding: utf-8 -*-
from __future__ import division
#try:
from runtime import *
#except ImportError:
    #from smop.runtime import *
def ecdfo_finish_(nb=None,mi=None,me=None,info=None,options=None,values=None,*args,**kwargs):
    """
    %
% Prints output status
%
%-----------------------------------------------------------------------

% Authors: Jean Charles Gilbert, INRIA.
%      and Anke Troeltzsch, DLR.
%
% Copyright 2008, 2009, INRIA. 2013, DLR.
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

%===============================================================================
    """
#    varargin = cellarray(args)
#    nargin = 6-[nb,mi,me,info,options,values].count(None)+len(args)

    if options.verbose >= 2:
        fprintf_(options.fout,'%s\n'%(values.dline))
        fprintf_(options.fout,'  Exit code %i: '%(info.flag))
        if values.success == info.flag:
            fprintf_(options.fout,'converged')
        elif values.fail_on_argument == info.flag:
            fprintf_(options.fout,'wrong argument')
        elif values.fail_on_problem == info.flag:
            fprintf_(options.fout,'unaccepted problem structure')
        elif values.fail_on_simul == info.flag:
            fprintf_(options.fout,'error when calling the simulator')
        elif values.stop_on_simul == info.flag:
            fprintf_(options.fout,'the simulator wants to stop')
        elif values.stop_on_max_iter == info.flag:
            fprintf_(options.fout,'max iteration reached')
        elif values.stop_on_max_simul == info.flag:
            fprintf_(options.fout,'max simulation reached')
        elif values.stop_on_dxmin == info.flag:
            fprintf_(options.fout,'too small variation in x (dxmin active)')
        elif values.fail_on_non_decrease == info.flag:
            fprintf_(options.fout,'the merit function cannot be decreased')
        elif values.fail_on_ascent_dir == info.flag:
            fprintf_(options.fout,'ascent direction encountered in linesearch')
        elif values.fail_on_max_ls_iter == info.flag:
            fprintf_(options.fout,'too many stepsize trials in linesearch')
        elif values.fail_on_ill_cond == info.flag:
            fprintf_(options.fout,'ill conditioning')
        elif values.fail_on_null_step == info.flag:
            fprintf_(options.fout,'null step d is solution of %0i QPs'%(values.max_null_steps + 1))
        elif values.fail_on_infeasible_QP == info.flag:
            fprintf_(options.fout,'infeasible QP')
        elif values.fail_on_unbounded_QP == info.flag:
            fprintf_(options.fout,'unbounded QP')
        elif values.fail_strange == info.flag:
            fprintf_(options.fout,'strange failure, call a guru')
        elif values.stop_on_small_trust_region == info.flag:
            fprintf_(options.fout,'trust region radius small')
        fprintf_(options.fout,'\n')
        fprintf_(options.fout,'%s\n'%(values.dline))
        fprintf_(options.fout,'  Final function value                     %12.5e\n'%(info.f))
        fprintf_(options.fout,'  Optimality conditions:\n')
        if nb + mi + me > 0:
            fprintf_(options.fout,'  . gradient of the Lagrangian (inf norm)  %11.5e\n'%(info.glagn))
            fprintf_(options.fout,'  . feasibility                            %11.5e\n'%(info.feasn))
            if nb + mi:
                fprintf_(options.fout,'  . complementarity                        %11.5e\n'%(info.compl))
        else:
            fprintf_(options.fout,'  Gradient of the cost function (inf norm)  %11.5e\n'%(info.glagn))
        fprintf_(options.fout,'  Counters:\n')
        fprintf_(options.fout,'  . nbr of iterations                   %4i\n'%(info.niter))
        fprintf_(options.fout,'  . nbr of function evaluations         %4i\n'%(info.nsimul[1] + info.nsimul[3]))
        fprintf_(options.fout,'  . nbr of gradient evaluations         %4i\n'%(info.nsimul[2] + info.nsimul[3]))
    fprintf_(options.fout,'%s\n'%(values.sline))
    return
