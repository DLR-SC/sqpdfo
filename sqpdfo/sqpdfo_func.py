# -*- coding: utf-8 -*-

from sqpdfo.runtime import *
import sqpdfo.sqpdfo_global_variables as glob
from numpy import array, zeros, concatenate, zeros_like


def sqpdfo_func_(x=None,*args,**kwargs):
    """
    #-----------------------------------------------------------------------
    # Computation of f, ci, ce
    #-----------------------------------------------------------------------
    """

    # Initialization
    c=array([])
    dummy=array([])
    msgf = 0
    msgc = 0
    msg = 0

    ffunc = glob.get_filename_f()
    cfunc = glob.get_filename_cons()

    # GTlab problem
    #scalefacX = array([[1000., 0.1, 10., 1., 1.]]).T
    #x = copy(x) / scalefacX

    xc = copy(x.T)
    xc = list(xc[0])
    try:
        f = ffunc(xc)
    except:
        f = np.inf
        msgf = 1
    if cfunc != '':
        try:
            c = cfunc(xc)
            c = array(c).reshape(-1, 1)
        except:
            msgc = 1
    if msgf > 0 or msgc > 0:
        msg = 'Error: error during calculating f (and/or c) in '
        msg = msg + 'user-defined problem !'


    return msg,f,dummy,c
   
