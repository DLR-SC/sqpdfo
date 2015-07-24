# -*- coding: utf-8 -*-
"""
Created on Fri Nov 07 15:15:41 2014
%
% [info] = sqplab_badsimul (outdic,info,options,values)
%
% Print a message, modify info, and return

%-----------------------------------------------------------------------
%
% Author: Jean Charles Gilbert, INRIA.
%
% Copyright 2008, 2009, INRIA.
%
% SQPlab is distributed under the terms of the Q Public License version
% 1.0.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the Q Public
% License version 1.0 for more details.
%
% You should have received a copy of the Q Public License version 1.0
% along with this program.  If not, see
% <http://doc.trolltech.com/3.0/license.html>.
%
%-----------------------------------------------------------------------
@author: jaco_da
"""

# Autogenerated with SMOP version 
# c:\Users\jaco_da\AppData\Local\Continuum\Anaconda\Scripts\smop-script.py sqplab_badsimul.m

#from __future__ import division
#try:
from runtime import *
#except ImportError:
#    from smop.runtime import *


def sqplab_badsimul_(outdic=None,info_=None,options=None,values=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 4-[outdic,info,options,values].count(None)+len(args)
    
    info=copy_(info_)

    if outdic == 2:
        if options.verbose:
            fprintf_(options.fout,char('\\n\\n### sqplab: the simulator wants to stop\\n'))
        info.flag=values.stop_on_simul
    else:
        if outdic > 2:
            if options.verbose:
                fprintf_(options.fout,char('\\n\\n### sqplab: error with the simulator (outdic = %0i)\\n'),outdic)
            info.flag=values.fail_on_simul
    return info

#def sqplab_badsimul_(outdic=None,info=None,options=None,values=None,*args,**kwargs):
#    #varargin = cellarray(args)
#    #nargin = 4-[outdic,info,options,values].count(None)+len(args)
#
#    if outdic == 2:
#        if options.verbose:
#            fprintf_(options.fout,char('\\n\\n### sqplab: the simulator wants to stop\\n'))
#        info.flag=values.stop_on_simul
#    else:
#        if outdic > 2:
#            if options.verbose:
#                fprintf_(options.fout,char('\\n\\n### sqplab: error with the simulator (outdic = %0i)\\n'),outdic)
#            info.flag=values.fail_on_simul
#    return info
