# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 15:15:59 2014

@author: jaco_da
"""
from runtime import isempty_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_

def ecdfo_swap_in_Y_(i,j,QZ,RZ,Y,ind_Y,fY,ciY,ceY,xbase,whichmodel,scale,shift_Y,Delta,normgx,kappa_ill, nargout=None):
	if (i > j):
		ii=j
		jj=i
	elif (i < j):
            ii=i
            jj=j
	else:
            return QZ,RZ,Y,ind_Y,fY,ciY,ceY,xbase,scale
		
	y=Y[:,ii]
	Y[:,ii]=Y[:,jj]
	Y[:,jj]=y												
	
	ind=ind_Y[ii]
	ind_Y[ii]=ind_Y[jj]
	ind_Y[jj]=ind
	
	f=fY[ii]
	fY[ii]=fY[jj]
	fY[jj]=f
	
	if not isempty_(ciY):
		ci=ciY[:,ii]
		#print "ciY", ciY								
		#print "type ciY", type(ciY)								
		#print "jj", jj								
		ciY[:,ii]=ciY[:,jj]
		ciY[:,jj]=ci
		
	if not isempty_(ceY):
		ce=ceY[:,ii]
		ceY[:,ii]=ceY[:,jj]
		ceY[:,jj]=ce
		

	QZ,RZ,xbase,scale=bcdfo_build_QR_of_Y_(Y,whichmodel,shift_Y,Delta,normgx,kappa_ill,nargout=4)
	return QZ,RZ,Y,ind_Y,fY,ciY,ceY,xbase,scale		