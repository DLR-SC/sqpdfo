# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 07:57:18 2015

@author: lien_ol


This file is used to set and get the global variables used in the ECDFO algorithm.

List of variables :
    prob : the problem number
    threshold : a threshold, typically 1e-08
    
    
Variables are not initialized here so that we are sure that the user set them himself, and do not forget that these variables exist (this is important at least for the tests)
"""


def set_prob(value):
	global prob
	prob = value


def get_prob():
	return prob
 
def set_threshold(value):
	global threshold
	threshold = value


def get_threshold():
    return threshold
 
