# -*- coding: utf-8 -*-

from runtime import *
from ecdfo_func import *
from ecdfo_global_variables import *
from numpy import array




def evalfgh_(key=None,xy=None,lm=None,*args,**kwargs):
    """
    #
# [msg] = evalfgh (key,xy);
# [msg,f,ci,ce] = evalfgh (key,xy)
# [msg,hl] = evalfgh (key,xy,lm)
#
# On entry:
#   key switches between three options, possible values are:
#      1: free 
#      2: compute f, ci, and ce
#      5: compute hl
#   xy: variables to optimize 
#   lm: KKT (or Lagrange) multiplier associated with the constraints
#     lm(1:n) is associated with bound constraints on xy 
#     lm(n+1:n+mi) is associated with the inequality constraints 
#     lm(n+mi+1:n+mi+me) is associated with the equality constraints 
#
# On return
#   msg describes the result of the simulation
#     0: the required computation has been done
#    -1: xy is out of an implicit domain
#    -2: stop the optimization please (something wrong)
#   f: cost-function value (here the potential energy of the chain)
#   ci: inequality constraint value (here the floor contraint values)
#   ce: equality constraint value (here the gap between the bar lengths and
#     there required values)
#   hl: Hessian of the Lagrangian
    """

    nargin = 3-[key is None,xy is None,lm is None].count(True)+len(args)

    # Global and persistent variables


    #If somewhere we have done set_fileoutput (in the tests for instance), then
    #we get this value, otherwise we take '1' by default.
    try:
        fileoutput = get_fileoutput()
    except:
        fileoutput=1  
    simul_not_initialized = get_simul_not_initialized()
    # On the output arguments

    msg=array([])
    out2=array([])
    out3=array([])
    out4=array([])
    # Set global and persistent variables (only once)

    if simul_not_initialized:
        # fid of the file with the xy coordinates of the hanging chain

        set_foutxy(fopen_('results.out','w'))
        # iteration counter

        set_iter(0)
        # initialisation has been done

        set_simul_not_initialized(0)
    # Call the appropriate function, depending on the value of key
    # free call

    if key == 1:
        if nargin < 1:
            fprintf_(fileoutput,'\n(simulopt) >>> not enough input arguments (%0i < 1) with key = %0i\n\n'%(nargin,key))
            msg=- 2
            return msg,out2,out3,out4
    # functions : f, ci, ce

    elif key == 2:
        if nargin < 2:
            fprintf_(fileoutput,'\n(simulopt) >>> not enough input arguments (%0i < 2) with key = %0i\n\n'%(nargin,key))
            msg=- 2
            return msg,out2,out3,out4
        msg,out2,out3,out4=ecdfo_func_(xy,nargout=4)
   # requiring hl

    elif key == 5:
        if nargin < 3:
            fprintf_(fileoutput,'\n(simulopt) >>> not enough input arguments (%0i < 3) with key = %0i\n\n'%(nargin,key))
            msg=- 2
            return msg,out2,out3,out4
        if nargout < 2:
            fprintf_(fileoutput,'\n(simulopt) >>> not enough output arguments (%0i < 2) with key = %0i\n\n'%(nargout,key))
            msg=- 2
            return msg,out2,out3,out4
        msg,out2=ecdfo_hessian_lagr_(xy,lm,nargout=2)
    else:
        fprintf_(fileoutput,'\n(simulopt) >>> unexpected value of key (=%i)\n\n'%(key))
        msg=- 2
        return msg,out2,out3,out4
    return msg,out2,out3,out4
