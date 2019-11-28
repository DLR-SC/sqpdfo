# -*- coding: utf-8 -*-

from runtime import *
import sqpdfo_global_variables as glob
from bcdfo_computeP import bcdfo_computeP_
from bcdfo_hessP import bcdfo_hessP_
from sqpdfo_bfgs_update import sqpdfo_bfgs_update_
from copy import copy
from numpy import zeros, ones, concatenate

def sqpdfo_computeHessian_(simul=None,x=None,null_step=None,constrained_pbl=None,\
    lm=None,M=None,n=None,me=None,mi=None,s=None,gx=None,gci=None,gce=None,\
    info_=None,options=None,values=None,fcmodel=None,Y=None,fY=None,ciY=None,\
    ceY=None,sigma=None,scale=None,shift_Y=None,QZ=None,RZ=None,\
    whichmodel=None,ind_Y=None,i_xbest=None,m=None,*args,**kwargs):

###############################################################################
# approximation of the Hessian of the Lagrangian
###############################################################################

    info=copy(info_)
    nbr_slacks = glob.get_nbr_slacks()
    sl = glob.get_slacks()

    lagfY=zeros(size_(Y,2))
    norm_ceY = zeros(size_(Y,2))
    norm_lm = zeros(size_(Y,2))

    pc=1.0
    if not null_step:

        # which Hessian approximation : BFGS <---> DFO model Hessian
        if options.hess_approx == values.bfgs:
            if options.verbose >= 4:
                if constrained_pbl:
                    fprintf_(options.fout,'\nBFGS update:\n')
                else:
                    fprintf_(options.fout,'\nBFGS inverse update:\n')

            # update y with the data at the previous point
            y=- info.g
            if me:
                y=y - info.ae.T.dot(lm[n + mi:n + mi + me])

            # save gradient of function and constraints
            info.g=gx
            info.ai=gci
            info.ae=gce

            # update y with the data at the new point
            y=y + info.g
            if me:
                y=y + info.ae.T.dot(lm[n + mi:n + mi + me])

            # compute the BFGS approximation of the Hessian of the Lagrangian 
            # (with Powell corrections)
            if options.bfgs_restart > 0 \
            and mod_(info.nsimul[1],options.bfgs_restart) == 0:
                szM = size_(M)
                M=eye_(szM[0])
                pc=2.0
            else:
                first=0
                if info.niter == 1:
                    first=1
                M,pc,info,values=\
                sqpdfo_bfgs_update_(M,y,s,first,info,options,values,nargout=4)
                if info.flag == values.fail_unexpected:
                    M=eye_(size_(M)[0])
                    M,pc,info,values=\
                    sqpdfo_bfgs_update_(M,y,s,first,info,options,values,nargout=4)
                    if info.flag == values.fail_unexpected:
                        disp_('sqpdfo_computeHessian: fail_unexpected ')
                        return M,pc,info

        elif options.hess_approx == values.model:
            info.g=gx
            info.ai=gci
            info.ae=gce

            cur_degree=size_(Y,2)

            # constrained case
            
            if not isempty_(ceY):

                #  Compute Lagrange function values for all points

                if nbr_slacks > 0:
                    for i in range(0,cur_degree):
                        norm_ceY[i] = norm_(ceY[:,i] - \
                                   concatenate((zeros((me-nbr_slacks,1)),sl**2)).T[0])   
                        norm_lm[i] = lm[n:n+me].T.dot(ceY[:,i] - \
                                   concatenate((zeros((me-nbr_slacks,1)),sl**2)).T[0])
                else:
                    for i in range(0,cur_degree):
                        norm_ceY[i] = norm_(ceY[:,i])
                        norm_lm[i] = lm[n:n+me].T.dot(ceY[:,i])

                lagfY = fY + norm_lm

                # Compute model of the Lagrangian
                
                model=\
                bcdfo_computeP_(QZ,RZ,Y,lagfY.reshape(1,-1),\
                               whichmodel,fcmodel[[0],:],ind_Y,\
                               i_xbest,m,gx,scale,shift_Y)
                
                # Compute Hessian of the Lagrangian
                
                M=bcdfo_hessP_(model,x,x,scale,shift_Y)
                
            else:
                M=bcdfo_hessP_(fcmodel[[0],:],x,x,scale,shift_Y)

        else:
            if options.verbose:
                fprintf_(options.fout,'\n### sqpdfo_computeHessian: options.hess_approx',\
                ' not recognized\n\n')
            info.flag=values.fail_on_argument
            return M,pc,info
    return M,pc,info
