# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 10:12:56 2015

@author: jaco_da
"""

from runtime import *

sigma = matlabarray([3])

sigmab = 6.27 # 0.05
sigma = max_(sigmab, 1.5 * sigma )

print "sigma\n", sigma