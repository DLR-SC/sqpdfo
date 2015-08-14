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

    def test_runtime_isvector_(self):
          M = matlabarray([[10, 11, 12]])
          M2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
          M3 = matlabarray([[10],[15], [18]])
          self.assertTrue(isvector_(M))
          self.assertFalse(isvector_(M2))
          self.assertTrue(isvector_(M3))
#        print "lboundsafter\n", lbounds

    #Thoses commented tests were just to see if the functions computed solutions like we want : eig_ and diag_ normally work
#like the equivalent matlab functions

#    def test_runtime_eig_(self):
#           M = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
#           print "eigenvalues of M : \n", eig_(M)
#           print "eigenvalues and eigenvectors of M : \n", eig_(M,2)
#           

#          
#    def test_runtime_diag_(self):
#          M = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
#          V=diag_(M)
#          print "diag of a matrix M : \n", V
#          M=diag_(V)
#          print "diag of a vector V : \n", M
#          
#          
#          
          
          #EXAMPLES OF UTILISATION OF MATLABARRAYS IN PYTHON :

    def test_runtime_boolean_indexing1(self):
        B2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        indfree2 = matlabarray([[0,1,2], [3,4,5],[6, 7,8]])
        #ilb2 = matlabarray([[True, False, True], [False, False, False], [True, False, True]])
        ilb2 = matlabarray([[False, False, False], [False, False, False], [False, False, False]])

#
#        print "B2\n", B2
#        print "indfree2\n", indfree2
#        print "ilb2\n", ilb2
#
#
        ret = indfree2[ilb2]
#        print "indfree2[ilb2]\n", ret
        self.assertTrue(isempty_(ret))

        B2[ilb2] = indfree2[ilb2]
        self.assertTrue(isequal_(B2,matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])))

#        print "B2\n", B2
#        
#        
    def test_runtime_boolean_indexing2(self):
        B2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        indfree2 = matlabarray([[0,1,2], [3,4,5],[6, 7,8]])
        #ilb2 = matlabarray([[True, False, True], [False, False, False], [True, False, True]])
        ilb2 = matlabarray([[True, False, True], [False, False, False], [True, False, True]])


#        print "B2\n", B2
#        print "indfree2\n", indfree2
#        print "ilb2\n", ilb2
#

        ret = indfree2[ilb2]
#        print "indfree2[ilb2]\n", ret
        self.assertTrue(isequal_(ret,matlabarray([[0,2,6,8]])))

        B2[ilb2] = indfree2[ilb2]
        self.assertTrue(isequal_(B2,matlabarray([[0, 11, 2],[13, 14, 15], [6, 17, 8]])))

#        print "B2\n", B2
#        
    def test_runtime_boolean_indexing3(self):
        B2 = matlabarray([[10, 11, 12, 13, 14, 15, 16, 17, 18]]).T
        indfree2 = matlabarray([[0,1,2, 3,4,5, 6, 7,8]]).T
        #ilb2 = matlabarray([[True, False, True], [False, False, False], [True, False, True]])
        ilb2 = matlabarray([[False, True, False, True, True, True, False, True, False]]).T


#        print "B2\n", B2
#        print "indfree2\n", indfree2
#        print "ilb2\n", ilb2


        ret = indfree2[ilb2]
        self.assertTrue(isequal_(ret,matlabarray([[1],[3],[4],[5],[7]])))

#        print "indfree2[ilb2]\n", ret
        
        B2[ilb2] = indfree2[ilb2]
        self.assertTrue(isequal_(B2,matlabarray([[10, 1, 12, 3, 4, 5, 16, 7, 18]]).T))

#        print "B2\n", B2
#
    def test_runtime_boolean_indexing4(self):
        lbounds = matlabarray([[-np.inf, -np.inf, -np.inf]]).T
        ilb = matlabarray([[False, True, False]]).T
        lb = matlabarray([[-0.5, 0.0, -np.inf]]).T
        indfree = matlabarray([[1,2,3]]).T

#        print "lbounds\n", lbounds
#        print "ilb\n", ilb
#        print "lb\n", lb
#        print "indfree\n", indfree

        lbounds[ilb]=lb[indfree[ilb]]
        self.assertTrue(isequal_(lbounds,[[-np.inf], [0.0], [-np.inf]]))


          
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
#    @unittest.skip("wonders of matlab")
#    def test050(self):
#        """
#        Compare to test060
#        octave> a=[]
#        a = []
#        octave> a(:,:)=99
#        a =  99
#        octave> a
#        a =  99
#        octave> size(a)
#        ans =
#        1   1
#        """
#        a = matlabarray()
#        a[:,:] = 99
#        self.assertTrue(isequal_(a.item(0), 99))

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

#    @unittest.skip("wonders of matlab")
#    def test080(self):
#        """
#        octave> a=[]
#        a = []
#        octave> a(1:end,5) = 5
#        a = [](0x5)        % notice 0x5
#        """

#    def test084(self):
#        """
#        octave> a=[]
#        a = []
#        octave> a(:,5) = 5
#        a =
#        0   0   0   0   5
#        """

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
        """Here are some more good examples of working indexing, for getting values or setting values
        Theses examples are certainly not exhaustive, but they helped me debug a lot and more may be
        coming in the future, if more bugs are to be found. Also, previous tests are also important"""



#EXAMPLES OF GETTING VALUES
        M = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])

        ind=matlabarray([[1,2,3]])
        self.assertTrue(isequal_(M[1,3],[[12]]))
        self.assertTrue(isequal_(M[:,2],[[11],[14],[17]]))
        self.assertTrue(isequal_(M[:,2:2],[[11],[14],[17]]))
        self.assertTrue(isequal_(M[1,:],[[10,11,12]]))
        self.assertTrue(isequal_(M[1:1,:],[[10,11,12]]))

        self.assertTrue(isequal_(M[ind,1],[[10],[13],[16]]))
        self.assertTrue(isequal_(M[3,ind],[[16,17,18]]))
        self.assertTrue(isequal_(M[ind.T,1],[[10],[13],[16]]))
        self.assertTrue(isequal_(M[ind.T,1:1],[[10],[13],[16]]))
        self.assertTrue(isequal_(M[3,ind.T],[[16,17,18]]))
        self.assertTrue(isequal_(M[3:3,ind.T],[[16,17,18]]))

        M = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        ind=matlabarray([1,2,3])
        self.assertTrue(isequal_(M[ind,1],[[10],[13],[16]]))
        self.assertTrue(isequal_(M[3,ind],[[16,17,18]]))


        M[1,3]=12
        self.assertTrue(isequal_(M[1,3],[[12]]))
#        
        
        #EXAMPLES WITH VECTORS
        V=matlabarray([[2,3,4,5]])
        self.assertTrue(isequal_(V[1:3],[[2,3,4]]))
        self.assertTrue(isequal_(V[:],[[2],[3],[4],[5]])) #same result in matlab
        self.assertTrue(isequal_(V[ind],[[2,3,4]]))
        self.assertTrue(isequal_(V[ind.T],[[2,3,4]]))

        
        V=matlabarray([[2],[3],[4],[5]])
        self.assertTrue(isequal_(V[1:3],[[2],[3],[4]]))
        self.assertTrue(isequal_(V[:],[[2],[3],[4],[5]]))
        self.assertTrue(isequal_(V[ind],[[2],[3],[4]]))
        self.assertTrue(isequal_(V[ind.T],[[2],[3],[4]]))


        ind=matlabarray([[True,False,True,False]])
        V=matlabarray([[2,3,4,5]])
        self.assertTrue(isequal_(V[ind],[[2,4]]))
        self.assertTrue(isequal_(V[ind.T],[[2,4]]))

        ind=matlabarray([[True,False,True,False]])
        V=matlabarray([[2,3,4,5]]).T
        self.assertTrue(isequal_(V[ind],[[2],[4]]))
        self.assertTrue(isequal_(V[ind.T],[[2],[4]]))

#EXAMPLES OF SETTING VALUES    


    #CHANGING COLUMNS OF THE MATRICES
        M2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M3 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M4 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M5 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M6 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
#        
        ind=matlabarray([[1,2,3]])
        M2[:,2]=matlabarray([[1,2,3]]).T       
        M3[:,2:2]=matlabarray([[1,2,3]]).T
        M4[ind,2]=matlabarray([[1,2,3]]).T
        M5[ind.T,2]=matlabarray([[1,2,3]]).T
##
        self.assertTrue(isequal_(M2, matlabarray([[10, 1, 12],[13, 2, 15], [16, 3, 18]])))
        self.assertTrue(isequal_(M3, matlabarray([[10, 1, 12],[13, 2, 15], [16, 3, 18]])))
        self.assertTrue(isequal_(M4, matlabarray([[10, 1, 12],[13, 2, 15], [16, 3, 18]])))
        self.assertTrue(isequal_(M5, matlabarray([[10, 1, 12],[13, 2, 15], [16, 3, 18]])))
#        
        ind=matlabarray([1,2,3])
        M6[ind,2]=matlabarray([[1,2,3]]).T
        self.assertTrue(isequal_(M6, matlabarray([[10, 1, 12],[13, 2, 15], [16, 3, 18]])))
#
#
        M2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M3 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M4 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M5 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M6 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        
        ind=matlabarray([[1,2,3]])
        M2[:,2]=matlabarray([[1,2,3]])       
        M3[:,2:2]=matlabarray([[1,2,3]])
        M4[ind,2]=matlabarray([[1,2,3]])
        M5[ind.T,2]=matlabarray([[1,2,3]])
#
        self.assertTrue(isequal_(M2, matlabarray([[10, 1, 12],[13, 2, 15], [16, 3, 18]])))
        self.assertTrue(isequal_(M3, matlabarray([[10, 1, 12],[13, 2, 15], [16, 3, 18]])))
        self.assertTrue(isequal_(M4, matlabarray([[10, 1, 12],[13, 2, 15], [16, 3, 18]])))
        self.assertTrue(isequal_(M5, matlabarray([[10, 1, 12],[13, 2, 15], [16, 3, 18]])))
#    
        ind=matlabarray([1,2,3])
        M6[ind,2]=matlabarray([[1,2,3]])
        self.assertTrue(isequal_(M6, matlabarray([[10, 1, 12],[13, 2, 15], [16, 3, 18]])))
        
#    ADDING COLUMNS
   
        M7 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        ind=matlabarray([[1,2,3]])
        M7[:,4]=ind
        self.assertTrue(isequal_(M7, matlabarray([[10, 11, 12,1],[13, 14, 15,2], [16, 17, 18, 3],])))
        
        M7 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M7[:,4]=ind.T
        self.assertTrue(isequal_(M7, matlabarray([[10, 11, 12,1],[13, 14, 15,2], [16, 17, 18, 3],])))
        
        M7 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M7[ind,4]=ind
        self.assertTrue(isequal_(M7, matlabarray([[10, 11, 12,1],[13, 14, 15,2], [16, 17, 18, 3],])))
        
        M7 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M7[ind,4]=ind.T
        self.assertTrue(isequal_(M7, matlabarray([[10, 11, 12,1],[13, 14, 15,2], [16, 17, 18, 3],])))
#        
     #CHANGING ROWS OF THE MATRICES
        ind=matlabarray([[1,2,3]])
        M2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M3 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M4 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M5 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M6 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        
        
        M2[1,:]=matlabarray([[1,2,3]])         
        M3[1:1,:]=matlabarray([[1,2,3]])
        M4[1,ind]=matlabarray([[1,2,3]])
        M5[1,ind.T]=matlabarray([[1,2,3]])
        self.assertTrue(isequal_(M2, matlabarray([[1, 2, 3],[13, 14, 15], [16, 17, 18]])))
        self.assertTrue(isequal_(M3, matlabarray([[1, 2, 3],[13, 14, 15], [16, 17, 18]])))        
        self.assertTrue(isequal_(M4, matlabarray([[1, 2, 3],[13, 14, 15], [16, 17, 18]])))
        self.assertTrue(isequal_(M5, matlabarray([[1, 2, 3],[13, 14, 15], [16, 17, 18]])))
        
        ind=matlabarray([1,2,3])
        M6[1,ind]=matlabarray([[1,2,3]])
        self.assertTrue(isequal_(M6, matlabarray([[1, 2, 3],[13, 14, 15], [16, 17, 18]])))
        
        ind=matlabarray([[1,2,3]])
        M2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M3 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M4 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M5 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        M6 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
        
        
        M2[1,:]=matlabarray([[1,2,3]]).T      
        M3[1:1,:]=matlabarray([[1,2,3]]).T
        M4[1,ind]=matlabarray([[1,2,3]]).T
        M5[1,ind.T]=matlabarray([[1,2,3]]).T
        self.assertTrue(isequal_(M2, matlabarray([[1, 2, 3],[13, 14, 15], [16, 17, 18]])))
        self.assertTrue(isequal_(M3, matlabarray([[1, 2, 3],[13, 14, 15], [16, 17, 18]])))        
        self.assertTrue(isequal_(M4, matlabarray([[1, 2, 3],[13, 14, 15], [16, 17, 18]])))
        self.assertTrue(isequal_(M5, matlabarray([[1, 2, 3],[13, 14, 15], [16, 17, 18]])))
        
        ind=matlabarray([1,2,3])
        M6[1,ind]=matlabarray([[1,2,3]]).T
        self.assertTrue(isequal_(M6, matlabarray([[1, 2, 3],[13, 14, 15], [16, 17, 18]])))
        
    
    #ADDING ROWS THIS WAY DOES NOT WORK DUE TO IMPLEMENTATION OF MATLAB ARRAYS AS F-CONTIGUOUS (i.e. Fortran contiguous, with a memory column by column)
   
#        M7 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
#        ind=matlabarray([[1,2,3]])
#        M7[4,:]=ind
#        self.assertTrue(isequal_(M7, matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18],[1,2,3]])))
#        
#        M7 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
#        M7[4,:]=ind.T
#        self.assertTrue(isequal_(M7, matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18],[1,2,3]])))
#        
#        M7 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
#        M7[4,:]=ind
#        self.assertTrue(isequal_(M7, matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18],[1,2,3]])))
#        
#        M7 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
#        M7[4,:]=ind.T
#        self.assertTrue(isequal_(M7, matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18],[1,2,3]])))
        
    #EXAMPLES WITH VECTORS
        
        #Row vector
        ind=matlabarray([[1,2,3]])
        V=matlabarray([[2,3,4,5]])
        V2=matlabarray([[2,3,4,5]])
        V3=matlabarray([[2,3,4,5]])
        V4=matlabarray([[2,3,4,5]])
        V[ind]=matlabarray([[8,9,10]])
        V2[1:3]=matlabarray([[8,9,10]])
        V3[:]=matlabarray([[8,9,10,11]])
        V4[1]=matlabarray([30])

        self.assertTrue(isequal_(V, [[8,9,10,5]]))
        self.assertTrue(isequal_(V2, [[8,9,10,5]]))
        self.assertTrue(isequal_(V3, [[8,9,10,11]]))
        self.assertTrue(isequal_(V4, [[30,3,4,5]]))

        ind=matlabarray([[1,2,3]])
        V=matlabarray([[2,3,4,5]])
        V2=matlabarray([[2,3,4,5]])
        V3=matlabarray([[2,3,4,5]])
        V4=matlabarray([[2,3,4,5]])
        V[ind]=matlabarray([[8,9,10]]).T
        V2[1:3]=matlabarray([[8,9,10]]).T
        V3[:]=matlabarray([[8,9,10,11]]).T
        V4[1]=matlabarray([30]).T

        self.assertTrue(isequal_(V, [[8,9,10,5]]))
        self.assertTrue(isequal_(V2, [[8,9,10,5]]))
        self.assertTrue(isequal_(V3, [[8,9,10,11]]))
        self.assertTrue(isequal_(V4, [[30,3,4,5]]))
        
        ind=matlabarray([[True,False,True,False]]).T
        V=matlabarray([[2,3,4,5]])
        V2=matlabarray([[2,3,4,5]])
        V3=matlabarray([[2,3,4,5]])
        V4=matlabarray([[2,3,4,5]])
        V[ind]=matlabarray([[8,9]])
        V2[ind]=matlabarray([[8,9]]).T
        V3[ind.T]=matlabarray([[8,9]])
        V4[ind.T]=matlabarray([[8,9]]).T
        self.assertTrue(isequal_(V, [[8,3, 9,5]]))
        self.assertTrue(isequal_(V, [[8,3, 9,5]]))
        self.assertTrue(isequal_(V, [[8,3, 9,5]]))
        self.assertTrue(isequal_(V, [[8,3, 9,5]]))
        
        #Column vector
        ind=matlabarray([[1,2,3]]).T
        V=matlabarray([[2,3,4,5]]).T
        V2=matlabarray([[2,3,4,5]]).T
        V3=matlabarray([[2,3,4,5]]).T
        V4=matlabarray([[2,3,4,5]]).T
        V[ind]=matlabarray([[8,9,10]])
        V2[1:3]=matlabarray([[8,9,10]])
        V3[:]=matlabarray([[8,9,10,11]])
        V4[1]=matlabarray([30])

        self.assertTrue(isequal_(V, [[8],[9],[10],[5]]))
        self.assertTrue(isequal_(V2, [[8],[9],[10],[5]]))
        self.assertTrue(isequal_(V3, [[8],[9],[10],[11]]))
        self.assertTrue(isequal_(V4, [[30],[3],[4],[5]]))
        
        ind=matlabarray([[1,2,3]]).T
        V=matlabarray([[2,3,4,5]]).T
        V2=matlabarray([[2,3,4,5]]).T
        V3=matlabarray([[2,3,4,5]]).T
        V4=matlabarray([[2,3,4,5]]).T
        V[ind]=matlabarray([[8,9,10]]).T
        V2[1:3]=matlabarray([[8,9,10]]).T
        V3[:]=matlabarray([[8,9,10,11]]).T
        V4[1]=matlabarray([30]).T

        self.assertTrue(isequal_(V, [[8],[9],[10],[5]]))
        self.assertTrue(isequal_(V2, [[8],[9],[10],[5]]))
        self.assertTrue(isequal_(V3, [[8],[9],[10],[11]]))
        self.assertTrue(isequal_(V4, [[30],[3],[4],[5]]))
        
        ind=matlabarray([[True,False,True,False]]).T
        V=matlabarray([[2,3,4,5]]).T
        V2=matlabarray([[2,3,4,5]]).T
        V3=matlabarray([[2,3,4,5]]).T
        V4=matlabarray([[2,3,4,5]]).T
        V[ind]=matlabarray([[8,9]])
        V2[ind]=matlabarray([[8,9]]).T
        V3[ind.T]=matlabarray([[8,9]])
        V4[ind.T]=matlabarray([[8,9]]).T
        self.assertTrue(isequal_(V, [[8],[3],[9],[5]]))
        self.assertTrue(isequal_(V2, [[8],[3],[9],[5]]))
        self.assertTrue(isequal_(V3, [[8],[3],[9],[5]]))
        self.assertTrue(isequal_(V4, [[8],[3],[9],[5]]))
        
    #ADDITION OF SCALARS TO A ROW VECTOR          
        V=matlabarray([[1]])
        V[:,2]=2
        self.assertTrue(isequal_(V, [[1,2]]))
        
        V=matlabarray([[1]])
        V[2]=2
        self.assertTrue(isequal_(V, [[1,2]]))        
        
        V=matlabarray([[2,3,4,5]])
        V[5]=6
        self.assertTrue(isequal_(V, [[2,3,4,5,6]]))
        V[6]=matlabarray([7])
        self.assertTrue(isequal_(V, [[2,3,4,5,6,7]]))
        V[7]=matlabarray([[8]])
        self.assertTrue(isequal_(V, [[2,3,4,5,6,7,8]]))
        V[8]=matlabarray([[9]]).T
        self.assertTrue(isequal_(V, [[2,3,4,5,6,7,8,9]]))
        V[1,9]=10
        self.assertTrue(isequal_(V, [[2,3,4,5,6,7,8,9,10]]))  
        
    #ADDITION OF SCALARS TO A COLUMN VECTOR DOES NOT WORK DUE TO IMPLEMENTATION OF MATLAB ARRAYS AS F-CONTIGUOUS (i.e. Fortran contiguous, with a memory column by column)
#        V=matlabarray([[2,3,4,5]]).T
#        V[5]=6
#        self.assertTrue(isequal_(V, matlabarray([[2,3,4,5,6]].T)))
#        V[6]=matlabarray([7])
#        self.assertTrue(isequal_(V, matlabarray([[2,3,4,5,6,7]].T)))
#        V[7]=matlabarray([[8]])
#        self.assertTrue(isequal_(V, matlabarray([[2,3,4,5,6,7,8]].T)))
#        V[8]=matlabarray([[9]]).T
#        self.assertTrue(isequal_(V, matlabarray([[2,3,4,5,6,7,8,9]].T)))
#        V[1,9]=10
#        self.assertTrue(isequal_(V,matlabarray([[2,3,4,5,6,7,8,9,10]].T))) 
#        
if __name__ == '__main__':
    unittest.main()