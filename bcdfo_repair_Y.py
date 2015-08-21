# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:58:06 2015

@author: lien_ol
"""


from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from bcdfo_find_new_yj import bcdfo_find_new_yj_
from bcdfo_replace_in_Y import bcdfo_replace_in_Y_
from copy import copy
    
def bcdfo_repair_Y_(QZ_=None,RZ_=None,Y_=None,Delta=None,farfact=None,farthr=None,closethr=None,eps_L=None,xbase_=None,lSolver=None,whichmodel=None,hardcons=None,xl=None,xu=None,indfree=None,stratLam=None,scale_=None,shift_Y=None,normgx=None,kappa_ill=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 20-[QZ,RZ,Y,Delta,farfact,farthr,closethr,eps_L,xbase,lSolver,whichmodel,hardcons,xl,xu,indfree,stratLam,scale,shift_Y,normgx,kappa_ill].count(None)+len(args)

    Y=copy(Y_)
    QZ=copy(QZ_)
    RZ=copy(RZ_)
    xbase=copy(xbase_)
    scale=copy(scale_)

    replaced=np.array([])
    verbose=0
    max_improve_loops=20
    if (verbose):
        disp_('--- enter bcdfo_repair_Y  for Delta = ',num2str_(Delta),' ---')
    n,p1=size_(Y,nargout=2)
    d=zeros_(p1,1)
    for j in range(1,p1): #it was initially in matlab range(2,p1)
        d[j]=norm_(Y[:,j] - Y[:,0])
#    dsorted,jsorted=sort_(d,'descend',nargout=2)
    dsorted=np.sort(d.reshape(-1))[::-1]
    jsorted=np.argsort(d.reshape(-1))[::-1]
    if (verbose):
      print 'd='
      print d
      print 'dsorted='
      print dsorted
      print 'jsorted='
      print jsorted
    for j in range(0,p1):
        if (dsorted[j] > farfact * (1 + eps_L) * Delta):
            jmax=jsorted[j]
            if (hardcons == 1):
                y,improvement,msgTR=bcdfo_find_new_yj_bc_(QZ,RZ,Y,jmax,Delta,eps_L,xbase,lSolver,whichmodel,xl,xu,indfree,stratLam,scale,shift_Y,nargout=2)
            else:
                y,improvement,msgTR=bcdfo_find_new_yj_(QZ,RZ,Y,jmax,Delta,eps_L,xbase,lSolver,whichmodel,scale,shift_Y,nargout=2)
            if (verbose):
                disp_('lambda = ',num2str_(improvement),' at j=',num2str_(jmax))
            QZ,RZ,Y,xbase,scale=bcdfo_replace_in_Y_(QZ,RZ,y,Y,jmax,xbase,whichmodel,scale,shift_Y,Delta,normgx,kappa_ill,nargout=5)
#            replaced=matlabarray([replaced,jmax])
            replaced  = np.append( replaced, jmax )
            d[jmax]=norm_(y - Y[:,[0]])
        else:
            Y_radius=dsorted[j]
            break
    if (verbose):
        print 'replaced='
        print replaced
        poised,Y_radius=bcdfo_poisedness_Y_(QZ,RZ,Y,eps_L,xbase,lSolver,whichmodel,hardcons,xl,xu,indfree,stratLam,scale,shift_Y,nargout=2)
        disp_('after distant replace: poisedness(Y) = ',num2str_(poised),' Y_radius  = ',num2str_(Y_radius))
    for k in range(1,max_improve_loops+1):
        maximprove=0
        for j in range(1,p1):
            if (hardcons == 1):
                y,improvement, msgTR=bcdfo_find_new_yj_bc_(QZ,RZ,Y,j,Delta,eps_L,xbase,lSolver,whichmodel,xl,xu,indfree,stratLam,scale,shift_Y,nargout=2)
            else:
                y,improvement, msgTR=bcdfo_find_new_yj_(QZ,RZ,Y,j,Delta,eps_L,xbase,lSolver,whichmodel,scale,shift_Y,nargout=2)
            if (verbose > 1):
                disp_(' ==> j = ',int2str_(j),' improve = ',num2str_(improvement))
                print 'y='
                print y
            if (improvement > maximprove):
                maximprove=copy(improvement)
                jmax=copy(j)
                ymax=copy(y)
        if (maximprove < closethr or jmax == 0):
            Y_radius=max_(d)
            if (verbose):
                print 'replaced='
                print replaced
                disp_('maximprove(small)= ',num2str_(maximprove),', jmax= ',int2str_(jmax))
                poised,Y_radius=bcdfo_poisedness_Y_(QZ,RZ,Y,eps_L,xbase,lSolver,whichmodel,hardcons,xl,xu,indfree,stratLam,scale,shift_Y,nargout=2)
                disp_('after everything: poisedness(Y) = ',num2str_(poised),' Y_radius  = ',num2str_(Y_radius))
                disp_('--- exit 1 bcdfo_repair_Y in round k=',num2str_(k),' ---')
            if (isempty_(replaced)):
                maximprove=0
            return QZ,RZ,Y,replaced,maximprove,Y_radius,xbase,scale
        if (verbose):
            disp_('maximprove= ',num2str_(maximprove),', jmax= ',int2str_(jmax))
        QZ,RZ,Y,xbase,scale=bcdfo_replace_in_Y_(QZ,RZ,ymax,Y,jmax,xbase,whichmodel,scale,shift_Y,Delta,normgx,kappa_ill,nargout=5)
        d[jmax]=norm_(ymax - Y[:,[0]])
        if (length_(find_(replaced == jmax)) == 0):
            replaced  = np.append( replaced, jmax )
    Y_radius=max_(d)
    if (verbose):
        print 'replaced='
        print replaced
        poised,Y_radius=bcdfo_poisedness_Y_(QZ,RZ,Y,eps_L,xbase,lSolver,whichmodel,hardcons,xl,xu,indfree,stratLam,scale,shift_Y,nargout=2)
        disp_('after everything: poisedness(Y) = ',num2str_(poised),' Y_radius  = ',num2str_(Y_radius))
        disp_('--- exit 2 bcdfo_repair_Y after round k=',num2str_(k),'  ---')
    return QZ,RZ,Y,replaced,maximprove,Y_radius,xbase,scale    