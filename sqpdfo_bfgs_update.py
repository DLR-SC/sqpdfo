# -*- coding: utf-8 -*-

from runtime import *
from copy import copy

def ecdfo_bfgs_update_(H=None,y=None,s=None,first=None,info_=None,options=None,values=None,*args,**kwargs):

###############################################################################
# ecdfo_bfgs_update_(H,y,s,first,info,options,values)
#
# This procedure computes the BFGS update of the matrix H, which is
# supposed to be a positive definite approximation of some Hessian.
# If y'*s is not sufficiently positive, Powell's correction is applied to y.
#
# Input:
#   H: matrix to update
#   y: change in gradient g
#   s: change in point x
#   first: if true, the procedure will initialize the matrix H
#
# Output:
#   H: updated matrix
###############################################################################

    info=copy(info_)
    eta=0.2
    n=length_(s)
    pc=1
    info.flag=values.success

    # prepare the update of H
    if norm_(s) == 0:
        info.flag=values.fail_unexpected
        if options.verbose >= 3:
            fprintf_(options.fout,'\n### ecdfo_bfgs_update: null step s\n\n')
        return H,pc,info,values

    ys=y.T .dot(s)
    Hs=H.dot(s)
    sHs=s.T.dot(Hs)
    if sHs <= 0:
        info.flag=values.fail_unexpected
        if options.verbose >= 3:
            fprintf_(options.fout,'\n### ecdfo_bfgs_update: BFGS Hessian approximation is not positive definite:\n')
            fprintf_(options.fout,"            s'*H*s = %g <= 0\n\n"%(sHs))
        return H,pc,info,values

    if (ys < eta * sHs):
        if options.verbose >= 3:
            fprintf_(options.fout,'\n### ecdfo_bfgs_update: curvature condition fails: ys < eta*sHs\n')
        # do Powell's correction
        pc=(1 - eta) * sHs / (sHs - ys)
        y=pc * y + (1 - pc) * Hs
        ys=y.T.dot(s)  # update ys, since y has changed
        if ys <= 0:
            info.flag=values.fail_unexpected
            if options.verbose >= 3:
                fprintf_(options.fout,"\n### ecdfo_bfgs_update: y'*s = %9.3e not positive despite correction:\n\n"%(ys).T)
            return H,pc,info,values
    else:
        if ys <= 0:
            info.flag=values.fail_unexpected
            if options.verbose >= 3:
                fprintf_(options.fout,"\n### ecdfo_bfgs_update: y'*s = %9.3e is nonpositive\n\n"%(ys).T)
            return H,pc,info,values

    # when matrix H has to be initialized
    if first:
        para=(y.T.dot(y)) / ys
        H=para * eye_(n)
        Hs=para*s
        sHs=s.T.dot(Hs)

    # update the matrix H by the BFGS formula
    H = H - (Hs.dot(Hs.T)) / sHs + (y.dot(y.T)) / ys

    return H,pc,info,values