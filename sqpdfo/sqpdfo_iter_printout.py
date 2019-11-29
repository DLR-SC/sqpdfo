# -*- coding: utf-8 -*-

from sqpdfo.runtime import *


def sqpdfo_iter_printout_(info=None,old_delta=None,norms=None,pc=None,itype=None,values=None,nb=None,mi=None,options=None,constrained_pbl=None,merit=None,*args,**kwargs):
    """
#  iteration printout    
    """

    if options.verbose < 4 and options.verbose >= 1:
        fprintf_(options.fout,'%4i  %4i'%(info.niter,info.nsimul[1] + info.nsimul[3]))
        fprintf_(options.fout,'  %+14.8e  %+14.8e  '%(info.f,merit))
        fid = fopen_('convhist.m', 'a')
        fprintf_(fid, '%6d  %6d  %+14.8e  %+14.8e  ' %(info.niter,info.nsimul[1]+info.nsimul[3],info.f,merit))
        if constrained_pbl:
            fprintf_(options.fout,'%10.4e  %10.4e'%(info.glagn,info.feasn))
            fprintf_(fid, '%10.4e  %10.4e  '%(info.glagn,info.feasn))
        else:
            fprintf_(options.fout,'%10.4e'%(info.glagn))
            fprintf_(fid,'%10.4e'%(info.glagn))

        if info.niter > 1:
            fprintf_(options.fout,'  %8.2e'%(old_delta))
            fprintf_(options.fout,'  %8.2e'%(norms))
            fprintf_(options.fout,'  %4.2f'%(pc))
            fprintf_(options.fout,'  %5s\n'%(itype))
        else:
            fprintf_(options.fout,'  \n')
        fprintf_(fid,'  \n')
        fclose_(fid)
    if options.verbose >= 4:
        fprintf_(options.fout,'%s\n'%(values.dline))
        fprintf_(options.fout,'iter %i,'%(info.niter))
        fprintf_(options.fout,'  cost %12.5e'%(info.f))
        if constrained_pbl:
            fprintf_(options.fout,',  glagn %11.5e,  feas %11.5e'%(info.glagn,info.feasn))
            if (nb + mi > 0):
                fprintf_(options.fout,',  compl %11.5e'%(info.compl))
        else:
            fprintf_(options.fout,',  grad %11.5e'%(info.glagn))
        if info.niter > 1:
            fprintf_(options.fout,'  %4.2f\n'%(pc))
        else:
            fprintf_(options.fout,'  \n')
    return
