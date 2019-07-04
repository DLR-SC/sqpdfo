import pycutest
import sys, shutil
sys.path.append("../")
import unittest
import helper
from runtime import *
from ecdfo_init_prob import ecdfo_init_prob_
from ecdfo_global_variables import *
from ecdfo import ecdfo_
import numpy as np

class Test_ecdfo_on_cuter_equalities(unittest.TestCase):

   def setUp(self):
      self.write_baseline_test_to_file = 0
      if self.write_baseline_test_to_file == 1:
          self.read_baseline_to_compare_with_new_test = 0
      else:
          self.read_baseline_to_compare_with_new_test = 1
      self.abs_tol=1e-5
      self.rel_tol=1e-8
      self.props_list = []
      self.results_list = []
      self.results_from_file = []
      
      # set options
      self.options = helper.dummyUnionStruct()
      self.options.hess_approx='model'
      self.options.bfgs_restart=0
      self.options.algo_descent='Powell'
      self.options.tol_grad=1e-04
      self.options.tol_feas=1e-04
      self.options.tol_bnds=1e-04
      self.options.miter=5000
      self.options.msimul=5000
      self.options.verbose=1
      self.options.whichmodel = 'subbasis'
      self.options.final_degree = 'quadratic'

      # open file with problem names and parameters
      self.f = open('cuterprobs_eq_only.dat','r')
      for line in self.f: 
          self.props_list.append(line)
      self.f.close()
      
      if self.read_baseline_to_compare_with_new_test:
          # open file with baseline results of CUTEst problems
          self.file = open('cutest_results_list_5000.dat','r')
          self.results_from_file = self.file.readlines()
          self.file.close()      

   def test_ecdfo_on_cuter_equalities(self):
      """
      Test ECDFO on CUTEst
      """
      #set parameters
      set_prob(1000)
      prob=get_prob()
      set_check_condition(0)     
      
      # check list for parameters
      #for i in range(0,len(self.props_list)):
      for i in range(0,5):
      
          print('**************************************************')
          print('problem number '+str(i))
          # unwrap problem names and parameters (in list from file)
          nl1 = self.props_list[i].split('\n')[0]
          nl2 = nl1.split('[')[1]
          nl3 = nl2.split(']')[0]
          nl4 = nl3.split(',')
          name = nl4[0].split('\'')[1]
          print(name)
          if (len(nl4) == 1):
              params = []
          elif (len(nl4) == 3):
              pname1 = nl4[1].split('\'')[1]
              pvalue1 = nl4[2]
              params = [pname1,pvalue1]
          elif (len(nl4) == 5):
              pname1 = nl4[1].split('\'')[1]
              pvalue1 = nl4[2]
              pname2 = nl4[3].split('\'')[1]
              pvalue2 = nl4[4]
              params = [pname1,pvalue1, pname2,pvalue2]
          
          # initialize problem in ECDFO
          set_prob_cuter(name,params)
          x,lb,ub,dxmin,li,ui,dcimin,infb,n,nb,mi,me,info=ecdfo_init_prob_(prob)
          print('CUTEst problem: ',name)
          print('n: ',n,', me: ',me)
          self.options.dxmin=dxmin
          lm=np.array([])
          
          # run problem with ECDFO
          x,lm,info=ecdfo_(x,lm,lb,ub,self.options)
          
          if self.write_baseline_test_to_file:
              # print results on the screen and write to string
              print('flag: ',info.flag)
              print('f: ',info.f)
              print('gnorm: ',info.glagn)
              print('cnorm: ',info.feasn)
              print('neval: ',info.nsimul[1])
              results_string = name+','+str(info.flag)+','+str(info.f)
              results_string = results_string+','+str(info.glagn)+','+str(info.feasn)
              results_string = results_string+','+str(info.nsimul[1])+'\n'
              self.results_list.append(results_string)
          
          if self.read_baseline_to_compare_with_new_test:
              # unwrap baseline results (from list results_from_file)
              reslist = self.results_from_file[i].split(',')
              resname = reslist[0]
              resflag = int(reslist[1])
              resf = float(reslist[2])
              resg = float(reslist[3])
              resc = float(reslist[4])
              resneval = int(float(reslist[5].split('\n')[0]))
          
              # comparison of actuasl test with baseline results from file
              self.assertEqual(info.flag,resflag)
              self.assertAlmostEqual(info.f,resf,places=5)
              self.assertAlmostEqual(info.glagn,resg,places=4)
              self.assertAlmostEqual(info.feasn,resc,places=4)
              self.assertEqual(info.nsimul[1],resneval)
          
          shutil.move('convhist.m', 'experimentsOnCutest/convhist_'+name+'.m')
          
      if self.write_baseline_test_to_file:
          # write results of CUTEst problems in file
          file = open('cutest_results_list.dat','w')
          for i in range(0,len(self.results_list)):
              file.write(self.results_list[i])
          file.close()           
      
      # check number of problems
      self.assertEqual(292,len(self.props_list))
   
if __name__ == '__main__':
    unittest.main()
       
