# -*- coding: utf-8 -*-
"""
Created on Fri Nov 07 15:17:32 2014

@author: jaco_da
"""
import sys
sys.path.append("../")

import unittest
from sqplab_badsimul import *
#import numpy as np
import helper

class dummyInfo():
    def __init__(self):
        self.flag = False

class Test_sqplab_badsimul(unittest.TestCase):
    """
          Reminder :
          This class is a test for sqplab_badsimul which prints a message, modifies info, and returns
    """
    def setUp(self):
        self.dummyOptions = helper.dummyOptions()
        self.dummyInfo = dummyInfo()
        self.dummyValues = helper.dummyValues()

    def test_sqplab_bad_simul1(self):
        """
               Test when oudic=2
        """
        self.dummyInfo = sqplab_badsimul_(2,self.dummyInfo,self.dummyOptions,self.dummyValues)
        
#        print self.dummyInfo.flag
        self.assertTrue(self.dummyInfo.flag==self.dummyValues.stop_on_simul)
    
    def test_sqplab_bad_simul2(self):
        """
               Test when oudic>2
        """
        self.dummyInfo = sqplab_badsimul_(3,self.dummyInfo,self.dummyOptions,self.dummyValues)
        
        #print self.dummyInfo.flag
        self.assertTrue(self.dummyInfo.flag==self.dummyValues.fail_on_simul)


if __name__ == '__main__':
    unittest.main()