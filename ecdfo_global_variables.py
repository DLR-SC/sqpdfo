# -*- coding: utf-8 -*-
"""
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

def set_prob_cuter(prob_cuter):
    global cproblem
    #Warning : here the CUTEr interface from this website has to be installed in order to use CUTEr problems :
    #http://fides.fe.uni-lj.si/~arpadb/software-pycuter.html. Thanks to Prof. Dr. Árpád Bűrmen
    from pycutermgr import clearCache, prepareProblem, importProblem
    clearCache(prob)
    prepareProblem(prob)
    cproblem=importProblem(prob)

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
    
def set_check_condition(value):
    global check_condition
    check_condition=value

    
def get_prob():
    return prob

def get_prob_cuter():
    return cproblem
 
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
    
def get_check_condition():
    return check_condition