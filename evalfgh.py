# -*- coding: utf-8 -*-

from runtime import *
from ecdfo_func import *
from ecdfo_global_variables import *
from numpy import array


def evalfgh_(key=None,xy=None,lm=None,*args,**kwargs):
    """
# [msg,f,ci,ce] = evalfgh (key,xy,lm)
#
# Inputs:
#   key switches options, possible values are:
#      2: compute f, ci, and ce
#   xy: variables to optimize
#   lm: KKT (or Lagrange) multiplier associated with the constraints
#     lm(1:n) is associated with bound constraints on xy 
#     lm(n+1:n+mi) is associated with the inequality constraints 
#     lm(n+mi+1:n+mi+me) is associated with the equality constraints 
#
# Outputs:
#   msg describes the result of the simulation
#     0: the required computation has been done
#    -1: xy is out of an implicit domain
#    -2: stop the optimization
#   f: function value
#   ci: inequality constraint values
#   ce: equality constraint values
    """

    nargin = 3-[key is None,xy is None,lm is None].count(True)+len(args)

    # Global and persistent variables

    #If somewhere we have done set_fileoutput (in the tests for instance), then
    #we get this value, otherwise we take '1' by default.
    try:
        fileoutput = get_fileoutput()
    except:
        fileoutput=1  

    # On the output arguments

    msg=0
    outf=array([])
    outci=array([])
    outce=array([])

    # Set global and persistent variables (only once)

    simul_not_initialized = get_simul_not_initialized()
    if simul_not_initialized:

        # iteration counter
        set_iter(0)

        # initialisation has been done
        set_simul_not_initialized(0)

    # call to ecdfo_func() to compute functions : f, ci, ce

    if key == 2:
        if nargin < 2:
            fprintf_(fileoutput,'\nevalfgh.py : not enough input arguments (%0i < 2) with key = %0i\n\n'%(nargin,key))
            msg=- 2
            return msg,outf,outci,outce
        msg,outf,outci,outce=ecdfo_func_(xy,nargout=4)
    else:
        fprintf_(fileoutput,'\nevalfgh.py : unexpected value of key (=%i)\n\n'%(key))
        msg=- 2
        return msg,outf,outci,outce
    return msg,outf,outci,outce
