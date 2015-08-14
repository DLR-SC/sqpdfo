# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 13:04:50 2015

@author: lien_ol
"""

from __future__ import division

try:
    from runtime import *
except ImportError:
    from smop.runtime import *


from bcdfo_evalZ import bcdfo_evalZ_

def bcdfo_evalP_(P=None,x=None,xbase=None,scale=None,shift_Y=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 5-[P,x,xbase,scale,shift_Y].count(None)+len(args)

    if (shift_Y):
        value=P.dot(  bcdfo_evalZ_((x - xbase) * scale[1],size_(P,2)))
    else:
        value=P.dot(  bcdfo_evalZ_(x,size_(P,2)))
    return value