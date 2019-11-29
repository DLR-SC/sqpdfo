# -*- coding: utf-8 -*-

from sqpdfo.runtime import *


def sqpdfo_finish_(nb=None,mi=None,me=None,info=None,options=None,values=None,*args,**kwargs):

###############################################################################
# Prints output status
###############################################################################

    if options.verbose >= 1:
        fprintf_(options.fout,'%s\n'%(values.dline))
        fprintf_(options.fout,'  Exit code %i: '%(info.flag))
        if values.success == info.flag:
            fprintf_(options.fout,'converged')
        elif values.fail_on_argument == info.flag:
            fprintf_(options.fout,'wrong argument')
        elif values.fail_on_problem == info.flag:
            fprintf_(options.fout,'unaccepted problem structure')
        elif values.fail_on_simul == info.flag:
            fprintf_(options.fout,'error during function evaluation')
        elif values.stop_on_simul == info.flag:
            fprintf_(options.fout,'function can not be evaluated')
        elif values.stop_on_max_iter == info.flag:
            fprintf_(options.fout,'maximum number of iterations')
        elif values.stop_on_max_simul == info.flag:
            fprintf_(options.fout,'maximum number of function evaluations')
        elif values.stop_on_dxmin == info.flag:
            fprintf_(options.fout,'too small variation in x (dxmin active)')
        elif values.fail_on_non_decrease == info.flag:
            fprintf_(options.fout,'the merit function cannot be decreased')
        elif values.fail_on_ill_cond == info.flag:
            fprintf_(options.fout,'ill conditioning')
        elif values.fail_on_null_step == info.flag:
            fprintf_(options.fout,'null step d is solution of %0i QPs'%(values.max_null_steps + 1))
        elif values.fail_on_infeasible_QP == info.flag:
            fprintf_(options.fout,'infeasible QP')
        elif values.fail_on_unbounded_QP == info.flag:
            fprintf_(options.fout,'unbounded QP')
        elif values.fail_unexpected == info.flag:
            fprintf_(options.fout,'unexpected (possible premature) termination')
        elif values.stop_on_small_trust_region == info.flag:
            fprintf_(options.fout,'trust region radius too small')
        fprintf_(options.fout,'\n')
        fprintf_(options.fout,'%s\n'%(values.dline))
        fprintf_(options.fout,'  Final function value                     %12.5e\n'%(info.f))
        fprintf_(options.fout,'  Optimality conditions:\n')
        if nb + me > 0:
            fprintf_(options.fout,'  . gradient of the Lagrangian (inf norm)  %11.5e\n'%(info.glagn))
            fprintf_(options.fout,'  . feasibility                            %11.5e\n'%(info.feasn))
        else:
            fprintf_(options.fout,'  Gradient of the objective function (inf norm)  %11.5e\n'%(info.glagn))
        fprintf_(options.fout,'  Counters:\n')
        fprintf_(options.fout,'  . nbr of iterations                   %4i\n'%(info.niter))
        fprintf_(options.fout,'  . nbr of function evaluations         %4i\n'%(info.nsimul[1] + info.nsimul[3]))
        fprintf_(options.fout,'%s\n'%(values.sline))
    return
