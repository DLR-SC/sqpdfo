# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 07:57:18 2015

@author: lien_ol


This file is used to set and get the global variables used in the ECDFO algorithm.

List of variables :
    prob : the problem number
    threshold : a threshold, typically 1e-08
    fileoutput : the output (a file or screen) of many fprintf_ functions
    simul_not_initialized : a variable indicating if simul has been initialized (0) or not (1)
    foutxy :  an output (a file or screen) which is not used anywhere yet (07.08.2015)
    _iter : count iterations in sqplab_tcg
    
    
Variables are not initialized here so that we are sure that the user set them himself, and do not forget that these variables exist (this is important at least for the tests)
"""


def set_prob(value):
	global prob
	prob = value


def set_threshold(value):
	global threshold
	threshold = value

 
def set_fileoutput(value):
    global fileoutput
    fileoutput=value
    

def set_simul_not_initialized(value):
    global simul_not_initialized
    simul_not_initialized=value
    
def set_foutxy(value):
    global foutxy
    foutxy=value
    
def set_iter(value):
    global _iter
    _iter=value

    
def get_prob():
    return prob
 
def get_threshold():
    return threshold
    
def get_fileoutput():
    return fileoutput
    
def get_simul_not_initialized():
    return simul_not_initialized
    
def get_foutxy():
    return foutxy
    
def get_iter():
    return _iter