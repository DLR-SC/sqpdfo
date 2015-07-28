# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 13:31:30 2015

@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from runtime import *
import numpy as np
#import helper

class Test_runtime(unittest.TestCase):

    def setUp(self):
        pass

#    def test_runtime_boolean_indexing1(self):
#        B2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
#        indfree2 = matlabarray([[0,1,2], [3,4,5],[6, 7,8]])
#        #ilb2 = matlabarray([[True, False, True], [False, False, False], [True, False, True]])
#        ilb2 = matlabarray([[False, False, False], [False, False, False], [False, False, False]])
#
#
#        print "B2\n", B2
#        print "indfree2\n", indfree2
#        print "ilb2\n", ilb2
#
#
#        ret = indfree2[ilb2]
#        print "indfree2[ilb2]\n", ret
#        
#        B2[ilb2] = indfree2[ilb2]
#        
#        print "B2\n", B2
#        
#        
#    def test_runtime_boolean_indexing2(self):
#        B2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
#        indfree2 = matlabarray([[0,1,2], [3,4,5],[6, 7,8]])
#        #ilb2 = matlabarray([[True, False, True], [False, False, False], [True, False, True]])
#        ilb2 = matlabarray([[True, False, True], [False, False, False], [True, False, True]])
#
#
#        print "B2\n", B2
#        print "indfree2\n", indfree2
#        print "ilb2\n", ilb2
#
#
#        ret = indfree2[ilb2]
#        print "indfree2[ilb2]\n", ret
#        
#        B2[ilb2] = indfree2[ilb2]
#        
#        print "B2\n", B2
#        
#    def test_runtime_boolean_indexing3(self):
#        B2 = matlabarray([[10, 11, 12, 13, 14, 15, 16, 17, 18]]).T
#        indfree2 = matlabarray([[0,1,2, 3,4,5, 6, 7,8]]).T
#        #ilb2 = matlabarray([[True, False, True], [False, False, False], [True, False, True]])
#        ilb2 = matlabarray([[False, True, False, True, True, True, False, True, False]]).T
#
#
#        print "B2\n", B2
#        print "indfree2\n", indfree2
#        print "ilb2\n", ilb2
#
#
#        ret = indfree2[ilb2]
#        print "indfree2[ilb2]\n", ret
#        
#        B2[ilb2] = indfree2[ilb2]
#        
#        print "B2\n", B2
#
#    def test_runtime_boolean_indexing4(self):
#        lbounds = matlabarray([[-np.inf, -np.inf, -np.inf]]).T
#        ilb = matlabarray([[False, True, False]]).T
#        lb = matlabarray([[-0.5, 0.0, -np.inf]]).T
#        indfree = matlabarray([[1,2,3]]).T
#
#        print "lbounds\n", lbounds
#        print "ilb\n", ilb
#        print "lb\n", lb
#        print "indfree\n", indfree
#
#        lbounds[ilb]=lb[indfree[ilb]]
#
#        print "lboundsafter\n", lbounds

    def test_runtime_eig_(self):
           M = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
           print "eigenvalues of M : \n", eig_(M)
           print "eigenvalues and eigenvectors of M : \n", eig_(M,2)
           
    def test_runtime_isvector_(self):
          M = matlabarray([[10, 11, 12]])
          M2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
          M3 = matlabarray([[10],[15], [18]])
          self.assertTrue(isvector_(M))
          self.assertFalse(isvector_(M2))
          self.assertTrue(isvector_(M3))
          
    def test_runtime_diag_(self):
          M = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
          V=diag_(M)
          print "diag of a matrix M : \n", V
          M=diag_(V)
          print "diag of a vector V : \n", M
#          
#          
#          
          
          #EXAMPLES OF UTILISATION OF MATLABARRAYS IN PYTHON :
          
    def test010(self):
        """Two-dimensional assignment to []"""
        a = matlabarray()
        a[a.shape[0]+1,[1,2,3]] = [123,456,789]
        a[a.shape[0]+1,[1,2,3]] = [123,456,789]
        a[a.shape[0]+1,[1,2,3]] = [123,456,789]
        self.assertTrue(isequal_(a, [[123,456,789],
                                     [123,456,789],
                                     [123,456,789]]))
    def test020(self):
        """Two-dimensional assignment to []
           Expand, use list [1,2,3] for indexing"""
        a = matlabarray()
        a[a.shape[0]+1,[1,2,3]] = 123
        a[a.shape[0]+1,[1,2,3]] = 123
        a[a.shape[0]+1,[1,2,3]] = 123
        self.assertTrue(isequal_(a, [[123,123,123],
                                     [123,123,123],
                                     [123,123,123]]))

    def test030(self):
        """Two-dimensional assignment to []
           Expand, use slice 1:3 for indexing"""
        a = matlabarray()
        #import pdb; pdb.set_trace()
        a[a.shape[0]+1,1:3] = 123
        a[a.shape[0]+1,1:3] = 123
        a[a.shape[0]+1,1:3] = 123
        #print a.shape
        self.assertTrue(isequal_(a, [[123,123,123],
                                     [123,123,123],
                                     [123,123,123]]))
    #@unittest.skip("FIXME")
    def test040(self):
        a = matlabarray()
        with self.assertRaises(IndexError):
            a[a.shape[0]+1,:] = 123
            a[a.shape[0]+1,:] = 123
            a[a.shape[0]+1,:] = 123
#            self.assertTrue(isequal_(a, [[123],
#                                         [123],
#                                         [123]]))
    @unittest.skip("wonders of matlab")
    def test050(self):
        """
        Compare to test060
        octave> a=[]
        a = []
        octave> a(:,:)=99
        a =  99
        octave> a
        a =  99
        octave> size(a)
        ans =
        1   1
        """
        a = matlabarray()
        a[:,:] = 99
        self.assertTrue(isequal_(a.item(0), 99))

    def test060(self):
        """One-dimensional assignment to empty array
        octave> a=[]
        a = []
        octave> a(:)=99
        a = []
        octave> a
        a = []
        """
        a = matlabarray()
        with self.assertRaises(IndexError):
            a[:] = 99
            self.assertTrue(isempty_(a))

    #@unittest.skip("wonders of matlab")
    def test062(self):
        """One-dimensional assignment to empty array
        octave> a=[]
        a = []
        octave> a(:)=[1 2 3]
        error: A(I) = X: X must have the same size as I
        """
        a = matlabarray()
        with self.assertRaises(Exception):
            a[:] = [1,2,3]

    #@unittest.skip("wonders of matlab")
    def test064(self):
        """One-dimensional assignment to empty array
        octave> a=[]
        a = []
        octave> a(1:3)=[1 2 3]
        1 2 3
        """
        a = matlabarray()
        a[1:3] = [1,2,3]
        self.assertTrue(isequal_(a, [1,2,3]))

    def test070(self):
        """
        wonders of matlab
        octave> c=[]
        c = []
        octave> c(1:end)=9
        c = []
        """
        a = matlabarray()
        a[1:a.shape[0]] = 9
        self.assertTrue(isempty_(a))

    @unittest.skip("wonders of matlab")
    def test080(self):
        """
        octave> a=[]
        a = []
        octave> a(1:end,5) = 5
        a = [](0x5)        % notice 0x5
        """

    def test084(self):
        """
        octave> a=[]
        a = []
        octave> a(:,5) = 5
        a =
        0   0   0   0   5
        """

    def test090(self):
        a = matlabarray([[11,22,33]])
        a[4] = 44
        self.assertTrue(isequal_(a,[[11,22,33,44]]))

    def test092(self):
        a = matlabarray([[11,22,33,44]])
        a[5:7] = [55,66,77]
        self.assertTrue(isequal_(a,[[11,22,33,44,55,66,77]]))

    def test094(self):
        a = matlabarray([[11,22,33,44,55,66,77]])
        a[[8,9]] = [88,99]
        self.assertTrue(isequal_(a,[[11,22,33,44,55,66,77,88,99]]))

    def test100(self):
        a = matlabarray([[1,3],
                         [2,4]])
        #a[: , a.shape[1]+1] = [5,6]
        a[: , 3] = matlabarray([[5,6]])
        self.assertTrue(isequal_(a,[[1,3,5],
                                    [2,4,6]]))

    def test110(self):
        a = zeros_(4,4,dtype=int)
        a[2:3,2:3] = 1
        #print a
        self.assertTrue(isequal_(a,[[0,0,0,0],
                                    [0,1,1,0],
                                    [0,1,1,0],
                                    [0,0,0,0]]))
                                    
    def test_indexing(self):
        """Here are good examples of working indexing, for getting values or setting values"""
        M = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M3 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M4 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M5 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M6 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M7 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M8 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])

        ind = matlabarray([[0,1,2], [3,4,5],[6, 7,8]])
        ind2=matlabarray([[1,2,3]])
        self.assertTrue(isequal_(M[1,3],[[12]]))
        self.assertTrue(isequal_(M[:,2],[[11],[14],[17]]))
        self.assertTrue(isequal_(M[:,2:2],[[11],[14],[17]]))
        self.assertTrue(isequal_(M[1,:],[[10,11,12]]))
        self.assertTrue(isequal_(M[1:1,:],[[10,11,12]]))

        self.assertTrue(isequal_(M[ind2,1],[[10],[13],[16]]))
        self.assertTrue(isequal_(M[3,ind2],[[16,17,18]]))
        self.assertTrue(isequal_(M[ind2.T,1],[[10],[13],[16]]))
        self.assertTrue(isequal_(M[ind2.T,1:1],[[10],[13],[16]]))
        self.assertTrue(isequal_(M[3,ind2.T],[[16,17,18]]))
        self.assertTrue(isequal_(M[3:3,ind2.T],[[16,17,18]]))

        M[1,3]=12
        self.assertTrue(isequal_(M[1,3],[[12]]))
        
        #WARNING : we cannot give here a row vector, this would not work.
        M2[:,2]=matlabarray([[1,2,3]]).T
        M3[:,2:2]=matlabarray([[1,2,3]]).T
        self.assertTrue(isequal_(M2, matlabarray([[10, 1, 12],[13, 2, 15], [16, 3, 18]])))
        self.assertTrue(isequal_(M3, matlabarray([[10, 1, 12],[13, 2, 15], [16, 3, 18]])))
        
        #WARNING : we cannot give here a column vector, this would not work.
        M4[1,:]=matlabarray([[1,2,3]])
        M5[1:1,:]=matlabarray([[1,2,3]])
        self.assertTrue(isequal_(M4, matlabarray([[1, 2, 3],[13, 14, 15], [16, 17, 18]])))
        self.assertTrue(isequal_(M5, matlabarray([[1, 2, 3],[13, 14, 15], [16, 17, 18]])))

                
        
if __name__ == '__main__':
    unittest.main()