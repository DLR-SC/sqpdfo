# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 11:09:58 2015

@author: jaco_da
"""

from runtime import matlabarray, isempty_
import numpy as np

#A = np.array([[1,2,3],[4,5,6],[7,8,9]])
#print "A\n", A
#B = np.array([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
#print "B\n", B

#ilb = np.array([[True, False, True]]).T#,[False, False, False],[False, False, False]])
#ilb = np.array([[True, False, True], [False, False, False], [True, False, True]])
#print "ilb\n", ilb

#indfree = np.array([[0,1,2], [4,5,6], [7,8,9]])
#print "indfree\n", indfree

#print "A[ilb]\n", A[ilb]
#print "A[ilb]\n", A.T[ilb.T]

#B.T[ilb.T] = A.T[ilb.T]
#print "B[ilb] = A[ilb]\n",B
#print "A[indfree[ilb]]\n", A[indfree[ilb]]

#print "A[indfree[ilb]].T\n", A[indfree[ilb]].T
#B[ilb] = A[ilb]
#print "B after\n", B

#---------------------------------

#B2 = matlabarray([[10, 11, 12, 13, 14, 15, 16, 17, 18]]).T#[[10, 11, 12],[13, 14, 15], [16, 17, 18]])
#indfree2 = matlabarray([[0,1,2, 3,4,5,6, 7,8]]).T
#ilb2 = matlabarray([[True, False, True, False, False, False, True, False, True]]).T#[[True, False, True], [False, False, False], [True, False, True]])
#ilb2 = matlabarray([[False, False, False, False, False, False, False, False, False]]).T#[[True, False, True], [False, False, False], [True, False, True]])

B2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
indfree2 = matlabarray([[0,1,2], [3,4,5],[6, 7,8]])
#ilb2 = matlabarray([[True, False, True], [False, False, False], [True, False, True]])
ilb2 = matlabarray([[False, False, False], [False, False, False], [False, False, False]])


print "B2\n", B2
print "indfree2\n", indfree2
print "ilb2\n", ilb2


ret = indfree2[ilb2]
print "indfree2[ilb2]\n", ret



if isempty_(ret):
	print "Isempty"
else:
	print "np.asarray(ret.T)", np.asarray(ret.T)[0]
	np.asarray(B2.T).__setitem__(np.asarray(ilb2.T), np.asarray(ret.T)[0])
	#np.asarray(B2.T).__setitem__(np.asarray(ilb2), np.asarray(ret.T)[0])
	
print B2