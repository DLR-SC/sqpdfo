# -*- coding: utf-8 -*-
"""
Created on Thu Nov 06 11:28:36 2014

@author: jaco_da
"""

# SMOP compiler runtime support library
# Copyright 2014 Victor Leikehman

# MIT license

"""
Main differences between Matlab matrices and numpy arrays
=========================================================

#. Array indices start with one, not with zero.  Accordingly, the last
   element of an N-element array is N, not N-1 as in C.

#. Matrix elements are ordered column-first, aka "Fortran" order.

#. Arrays auto-expand on out-of-bound lhs indexing.

#. In matlab, arrays can be updated before they are created::

      clear a
      a(17) = 42

   is legal in matlab, but not in numpy

#. Array data is not shared by copying or slice indexing. Instead there
   is copy-on-write.

#. There are no zero or one-dimensional arrays. Scalars are
   two-dimensional rather than zero-dimensional as in numpy.

#. Single subscript implies ravel.

#. Boadcasting rules are different

Coexistence of matlab matrices and numpy arrays
===============================================

#. Empty vector::

        []                  matlabarray()

#. Scalars are 1x1 matrices::

        17                  [[ 17 ]]

#. Rectangular char arrays::

        'hello world'       char('hello world')

#. Row vector::

        [1 2 3]             [[1, 2, 3]]

#. Column vector::
 
        [1;2;3]             [[1], [2], [3]]

#. Cell strings::

        cellstr('abc','hello',[97 98 99])       


(*) Such expressions _are_ legal in Octave.  TBD

"""
import warnings
from copy import copy
#warnings.simplefilter('error', RuntimeWarning)
import scipy
import numpy as np
import os,sys
import numpy #this is for the other files which import runtime and need numpy, for instance with numpy.linalg.solve
from numpy import inf, logical_or, logical_and #some other files which import runtine need inf
#import helper
try:
    from scipy.io import loadmat
except:
    pass

def isvector_or_scalar_(a):
    """
    one-dimensional arrays having shape [N],
    row and column matrices having shape [1 N] and
    [N 1] correspondingly, and their generalizations
    having shape [1 1 ... N ... 1 1 1].
    Scalars have shape [1 1 ... 1].
    Empty arrays dont count
    """
    try:
        return a.size and a.ndim-a.shape.count(1) <= 1
    except:
        return False
def isvector_(a):
    """
    one-dimensional arrays having shape [N],
    row and column matrices having shape [1 N] and
    [N 1] correspondingly, and their generalizations
    having shape [1 1 ... N ... 1 1 1]
    """
    try:
        return a.ndim-a.shape.count(1) == 1
    except:
        return False

def disp_(*args):
    out = ""
    for arg in args:
        try:
            out = out + str(arg)
        except Exception as e:
            print "Display Error:", e.message()
    print out#(args)

false = False

def false_(*args):
    if len(args) == 1:
        args += args
    return np.zeros(args,dtype=bool,order="F")

def find_(a,n=None,d=None,nargout=1):
    if d:
        raise NotImplementedError

    elif nargout == 1:
        i = np.flatnonzero(a).reshape(1,-1)
        if isempty_(i):
            return []
                                                
        if n is not None:
            i = i.take(range(n))
        return np.array(i)
    elif nargout == 2:
        i,j = np.nonzero(np.asarray(a))
        if n is not None:
            i = i.take(n)
            j = j.take(n)

        return np.array(i).reshape(1,-1), np.array(j).reshape(1,-1)
    
    else:
        raise NotImplementedError

def floor_(a):
    return np.floor_(np.asanyarray(a))

def fopen_(*args):
    try:
        fp = open(*args)
        assert fp != -1
        return fp
    except:
        return -1
                                
def fclose_(*args):
    if args[0] == -1:
        print "No file to close, continue..."
    else:
        args[0].close()

def fullfile_(*args):
    return os.path.join(*args)

def intersect_(a,b,nargout=1):
    if nargout == 1:
        c = sorted(set(a) & set(b))
        if isinstance(a,str):
            return "".join(c)
        elif isinstance(a,list):
            return c
        else:
            # FIXME: the result is a column vector if
            # both args are column vectors; otherwise row vector
            return np.array(c)
    raise NotImplementedError


def isempty_(a):
    if a is None or a==[] :
        return True
    try:
        return 0 in a.shape
    except AttributeError:
        return False

def isequal_(a,b):
    return np.array_equal(np.asanyarray(a),
                          np.asanyarray(b))

#def isvectorRow_(a):
#    if a.ndim==1:
#        return True
#    elif ndims_(a)==2:
#        if a.shape[0]==1 and a.shape[1]>=1:
#            return True
#            
#    return False
#    
#def isvectorColumn_(a):
#    if ndims_(a)==2:
#        if a.shape[0]>=1 and a.shape[1]==1:
#            return True
#            
#    return False
                          
def isscalar_(a):
    """np.isscalar returns True if a.__class__ is a scalar
    type (i.e., int, and also immutable containers str and
    tuple, but not list.) Our requirements are different"""
    try:
        return a.size == 1
    except AttributeError:
        return np.isscalar(a)

def length_(a):
    try:
        return max(np.asarray(a).shape)
    except ValueError:
        return 1

try:
    def load_(a):
        return loadmat(a) # FIXME
except:
    pass

def max_(a, d=None, nargout=None):
    """min_ and max_ function normally returns the same as matlab min and max in the following cases :
        min_(a,b) where a,b are integer, inf or nan
        min_(a) and min_(a, nargout=2) where a is an array containing numbers, inf or nan
        One case not implemented but not necessary yet is the following :
        min_(a,b) where both a and b are arrays containing NaNs. Python returns then NaNs in priority.
        (same for max_ obviously)
    """
    #The warnigns happens when NaNs are involved, but the function returns in any case what we want, so no need to print the warnings
#    warnings.simplefilter('ignore', RuntimeWarning)

    if isempty_(a):
        ret = np.array([])                    
    elif  d is not None:
        ret=np.fmax(a,d)
    else:
        ret=np.nanmax(a)

                                
    if nargout == 2:
        if isempty_(a):
            ret2 = np.array([])                                    
        else:
            ret2 = np.nanargmax(a)
#        warnings.simplefilter('default', RuntimeWarning)                                                                                              
        return ret, ret2
    else:
#        warnings.simplefilter('default', RuntimeWarning)                                              
        return ret    

def min_(a, d=None, nargout=None):#, nargout=0):
    #print "a", a
#    print "len(a)", len(a)     
    """min_ and max_ function normally returns the same as matlab min and max in the following cases :
        min_(a,b) where a,b are integer, inf or nan
        min_(a) and min_(a, nargout=2) where a is an array containing numbers, inf or nan
        One case not implemented but not necessary yet is the following :
        min_(a,b) where both a and b are arrays containing NaNs. Python returns then NaNs in priority.
        (same for max_ obviously)
    """
    
    
    #The warnigns happens when NaNs are involved, but the function returns in any case what we want, so no need to print the warnings
#    warnings.simplefilter('ignore', RuntimeWarning)
    if isempty_(a):
        ret = np.array([])                    
    elif d is not None:
        ret=np.fmin(a,d)
    else:
        ret = np.nanmin(a)   
    if nargout == 2:
        if isempty_(a):
            ret2 = np.array([])                                    
        else:
            ret2 = np.nanargmin(a)  
#        warnings.simplefilter('default', RuntimeWarning)                                              
        return ret, ret2
    else:
#        warnings.simplefilter('default', RuntimeWarning)
        return ret                                

    warnings.simplefilter()
def mod_(a,b):
    try:
        return a % b
    except ZeroDivisionError:
        return a
def ndims_(a):
    return a.ndim

def numel_(a):
    return a.size

def rand_(*args,**kwargs):
    if not args:
        return np.random.rand()
    if len(args) == 1:
        args += args
    try:
        return np.random.rand(np.prod(args)).reshape(args,order="F")
    except:
        pass

def ravel_(a):
    return np.asanyarray(a).reshape(-1,1)

def round_(a):
    return np.round(np.asanyarray(a))

def rows_(a):
    return np.asarray(a).shape[0]

def columns_(a):
    return np.asarray(a).shape[1]

def size_(a, b=0, nargout=1):
    """
    >>> size_(zeros_(3,3)) + 1
    matlabarray([[4, 4]])
    """
    s = a.shape
    #not sure what the 2 following lines are for
#    if s is ():
#        return 1 if b else (1,)*nargout
    # a is not a scalar
    try:
        if b:
            return s[b-1]
        else:
            return np.array(s) if nargout <= 1 else s
    except IndexError:
        return 1

def strread_(s, format="", nargout=1):
    if format == "":
        a = [float(x) for x in s.split()]
        return tuple(a) if nargout > 1 else np.asanyarray([a])
    raise NotImplementedError

def strrep_(a,b,c):
    if isinstance(a,str):
        return a.replace(b,c)
    raise NotImplementedError # cell arrays

def sum_(a, dim=None):
    if dim is None:
        return a.sum(0)
    else:
        return a.sum(dim-1)

true = True

def true_(*args):
    if len(args) == 1:
        args += args
    return np.ones(args,dtype=bool,order="F")


def zeros_(*args,**kwargs):

    if not args:
        return 0.0
    if len(args) == 1:
        if numel_(args[0])==1:
            args += args
        else:
            return np.zeros(args[0],**kwargs)
    return np.zeros(args,**kwargs)


def ones_(*args,**kwargs):

    if not args:
        return 1.0
    if len(args) == 1:
        if numel_(args[0])==1:
            args += args
        else:
            return np.ones(args[0],**kwargs)
    return np.ones(args,**kwargs)
                
#------------------------------------------------------------------------------
#                Added Functions Start here.
#------------------------------------------------------------------------------
# -> matrices given as arguments  are matlabarray. Otherwise these functions may not do their job. We are doing our best here so that
#    the linear algebra functions return the arguments like in matlab. Therefore we do modifications on the outputs of classical python functions.
# -> normally all the functions below return matlabarrays when given matlabarrays as inputs

def eig_(A, nargout=1):
    if nargout == 1:
        return np.linalg.eigvals(A)
    else:
        D,V = np.linalg.eig(A) #when A is a matlabarray, linalg.eig(A) returns D as a python array and V as a matlabarray
        return V,np.diag(D)
            
def cond_(A):
    return np.linalg.cond(A)
    
def svd_(A, full_matrices=1, *args,**kwargs):
    u,s,v= np.linalg.svd(A, full_matrices)
    return u, np.diag(s), v.T
    
def chol_(A,nargout=2):
    try: #there is in python a exception when the matrix is not positive definite, which does not occur on matlab
        R = np.linalg.cholesky(A).T
        p = 0
    except:
        R = np.array([])
        p = 1
    return R,p   
    
def qr_(A,nargout=2):
    return np.linalg.qr(A)
    
def inv_(A):
    return np.linalg.inv(A)
     
def norm_(A, order=2, axis=None):
    if isempty_(A):
        return 0
        
#for some reason, python does not give exactly (so precision differs
#after the 10nth variable or so, but still this makes differences...) the same 2-norm
#when it does it calcultes the norm of an array([a b c...]) and array([[a b c...]]).
#And it is closer to matlab when it is like array([[a b c..]]). Hence those 2 lines
    if isvector_or_scalar_(A):
        #A=A.reshape(-1,1) <- this also works
        A=A.reshape(1,-1)        
        
    if axis==None:
        return np.linalg.norm(A, order,axis)
    else:
        return np.linalg.norm(A, order,axis-1)

    
def pinv_(A):
    return np.linalg.pinv(A)
    
def solve_(A,b):
    return np.linalg.solve(A,b)

def fprintf_(fid,*args, **kwargs):
    """This is supposed to work approximately like the matlab fprintf function in the following case in the cases where
    fid=1 (returns output on the screen), fid=2 (returns output in standard error ouput) or fid is a file identifier"""
    out = ""
    for arg in args:
        try:
            out = out + str(arg)
        except Exception as e:
            print "Display Error:", e.message()
    if  fid==1:
        print out,
    elif fid==2:
        print >>sys.stderr, out,
    elif type(fid)==file:
        print >>fid, out,
    else:
        print fid+out,
        
def poly1d_(A, r=0):
    #Careful : A.r which gives then the roots of the polynom is an array instead of a matlab array
    return np.poly1d(A,r) 
    
def eye_(n):
    return np.eye(n)
    
def concatenate_(arrs, axis=0):
    copy_arrs=copy(arrs)
    for arr in arrs:
        if isempty_(arr):
            copy_arrs.remove(arr)
    if copy_arrs==[]:
        return np.array([])
    else:
        return np.concatenate(copy_arrs, axis)
    
def sign_(A):
    return np.sign(A)
    
def setdiff_(a,b):
    return np.setdiff1d(a,b)
    
def isfield_(obj, name):
    return obj.__dict__.has_key(str(name))
    
def lower___(strng):
    #print "Warning: Lower just returns the string as it is"
    return strng
    
def lower__(strng):
    return lower___(strng)
    
def regexprep___(string1, string2, string3):
    #print "Warning: regexprep__ just returns the string(1) as it is"
    return string1
    
def regexprep__(string1, string2, string3):
    return regexprep___(string1, string2, string3)

def strtrim___(strng):
    #print "Warning: strtrim___ just returns the string as it is"
    return strng

def strtrim__(strng):
    return strtrim___(strng)
    
def any_(A):
    return A.any()
    
def sqrt_(x):
    return np.sqrt(x)

def real_(x):
    ret = np.real(x)
    #print "type np real ret", type(ret)
    return ret

def imag_(x):
    ret = np.imag(x)
    return ret

def exp_(x):
    return np.exp(x)
    
def num2str_(num):
    return str(num)

def strcmp_(s1, s2):    
    #print "compare strings:"
    #print "s1", s1
    #print "s2", s2
    if s1 == s2:
        return 1
    else:
        return 0        

def null_(A, eps=1e-15):
    u, s, vh = scipy.linalg.svd(A)
    padding = max(0,np.shape(A)[1]-np.shape(s)[0])
    null_mask = np.concatenate(((s <= eps), np.ones((padding,),dtype=bool)),axis=0)
    null_space = scipy.compress(null_mask, vh, axis=0)
    return scipy.transpose(null_space)        
#def null_(A, eps=1e-15):
#    u, s, vh = scipy.linalg.svd(A)
#    null_mask = (s <= eps)
#    null_space = scipy.compress(null_mask, vh, axis=0)
#    return scipy.transpose(null_space)        
#def null_(A, eps=1e-3):
    #u,s,vh = np.linalg.svd(A,full_matrices=1,compute_uv=1)
    #print "u", u
    #print "s", s
    #print "vh", vh
    #null_space = np.compress(s <= eps, vh, axis=0)
    #return null_space.T

def int2str_(val):
    return str(val)
    
def full_(A):
    #print "full_ just returns A"
    return A

if __name__ == "__main__":
    import doctest
    doctest.testmod()

def tf_mapper(x):
    #print "x", x
    if x:
        return 1
    else:
        return 0
                
def logical_or_(a,b):
#    vtf_mapper = np.vectorize (tf_mapper)
#    ret = np.logical_or(a,b)                
#    return vtf_mapper(ret)
    return np.logical_or(a,b)
    
def logical_and_(a,b):
    return np.logical_and(a,b)
    
def strcat_(*args):
    ret = ""
    for arg in args:
        ret = ret + str(arg)
        
    #print "strcat_:", ret
    return char(ret)
    
    
def randn_(msg, number):
    print "Warning: randn does not do anything, because it's discouraged syntax"
    
def compare_array(x, y, abs_tol, rel_tol):
    """  Function which compares the matlabarray x and y and returns true if all the elements are
     approximately the same according to the absolute tolerance (abs_tol) and to the relative tolerance (rel_tol)
    """ 
    
    #Flatening of x and y
    x=x.flat[:];
    y=y.flat[:];     
     
    #Delete the indices wherever x and y are infinite with the same sign.
    indices=np.where(np.logical_or(np.logical_and(x==np.inf, y==np.inf) , np.logical_and(x==-np.inf, y==-np.inf) ))
    x=np.delete(x,indices)
    y=np.delete(y,indices)
     
    #compute error and relative error. sys.float_info.min=2.2250738585072014e-308  to avoid zero division
    error = x-y;
    try:
        rel_error=error/(x+np.ones(np.shape(x))*sys.float_info.min)  #This line may issue warnings, but the inf and nan are dealt with below on the code. This is not an issue.
    except:
        pass
    
    #Sets rel_error[i] = 0 if either |x[i]|<= 0 or |y[i]|<=abs_tol, since in this case it makes sense only  to look at the absolute error
    indices=np.where(np.logical_or(abs(x)<=abs_tol, abs(y)<=abs_tol) )
    rel_error[indices] = 0

    return (abs(error) < abs_tol).all() and (abs(rel_error) < rel_tol).all()