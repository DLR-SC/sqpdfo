# -*- coding: utf-8 -*-
import sys
import types

"""
This file is used to set and get the global variables used in the ECDFO algorithm.

List of variables :
    prob : the problem number
    threshold : a threshold, typically 1e-08
    fileoutput : the output (a file or screen) of many fprintf_ functions
    simul_not_initialized : a variable indicating if simul has been initialized (0) or not (1)
    foutxy :  an output (a file or screen) which is not used anywhere yet (07.08.2015)

Variables are not initialized here so that we are sure that the user set them himself, and do not forget that these variables exist (this is important at least for the tests)
"""


class ecdfo_global_variables:
    # default values, in case the global variable is not set by the user
    global prob
    prob = None
    global fileoutput
    fileoutput = 1
    global threshold
    threshold = 1e-8  # Threshold for violated bounds
    global check_condition
    check_condition = 1  # Check condition of the interpolation matrix regularly
    global simul_not_initialized
    simul_not_initialized = 1  # Initial value
    global filename_f
    filename_f = None
    global filename_cons
    filename_cons = None
    global nbr_slacks
    nbr_slacks = 0
    global slacks
    slacks = None

def set_prob(value):
	global prob
	prob = value

def set_prob_cuter(prob_cuter, params):
    global cproblem
    #Warning : here the pyCUTEst interface from this website has to be 
    #installed in order to use CUTEst problems :
    #https://jfowkes.github.io/pycutest/_build/html/index.html 
    #Thanks to Jaroslav Fowkes and Lindon Roberts
    import pycutest
    len_params = len(params)
    if len_params == 0:
        cproblem = pycutest.import_problem(prob_cuter,efirst=True)
    elif len_params == 2:
        cproblem = pycutest.import_problem(prob_cuter, sifParams={params[0]:float(params[1])},efirst=True)
    elif len_params == 4:
        cproblem = pycutest.import_problem(prob_cuter, sifParams={params[0]:float(params[1]),params[2]:float(params[3])},efirst=True)
    
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

def set_filename_f(value):
    global filename_f
    if value is None:
        sys.exit('Error: Definition of the objective function in func() is missing !')
    elif isinstance(value, types.FunctionType):
        filename_f = value
    else:
        sys.exit('Error: function handle for the objective is not of type function !')

def set_filename_cons(value):
    global filename_cons
    if isinstance(value, types.FunctionType):
        filename_cons = value
    elif value is '' or value is None:
        filename_cons = ''
    else:
        sys.exit('Error: function handle for the constraints is not of type function !')

def set_nbr_slacks(value):
    global nbr_slacks
    nbr_slacks = value
    
def set_slacks(value):
    global slacks
    slacks = value
    

def get_prob():
    if prob is None:
        sys.exit('Problem number is not set!\n'\
            'Please import ecdfo_global_variables and use set_prob(nbr)\n'\
            'where nbr is 1,...,5 for test examples or 100 for a user defined problem.')
    else:
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

def get_filename_f():
    return filename_f

def get_filename_cons():
    return filename_cons
    
def get_nbr_slacks():
    return nbr_slacks
    
def get_slacks():
    return slacks
