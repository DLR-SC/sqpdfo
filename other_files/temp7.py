# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 11:43:34 2015

@author: jaco_da
"""

from runtime import matlabarray
import numpy as np

free = matlabarray([[1,2,3,4]])
atb = matlabarray([[2]])

print "free\n", free
print "atb\n", atb

#free = np.delete(np.asarray(free), np.asarray(atb) - 1)
free[atb]=[]

print "free after\n", free
