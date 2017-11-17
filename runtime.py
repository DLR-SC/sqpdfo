# -*- coding: utf-8 -*-

"""

This file contains functions which serve mostly 2 purposes :
    -simply providing a shortcut for matlab functions, for instance real(A) in matlab becomes real_(A) in python, and
    real_(A) actually calls np.real(A)
    -imitating at best the output of matlab functions. For instance, U,S,V=svd function in matlab and python does not give exactly
    outputs with the same shapes : S is a vector in python and a diagonal matrix in matlab, and V.transpose in python = V in matlab.

"""
#import warnings
from copy import copy
#warnings.simplefilter('error', RuntimeWarning)
import scipy
import numpy as np
import sys


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
    """
    function printing out arguments
    """
    out = ""
    for arg in args:
        try:
            out = out + str(arg)
        except Exception as e:
            print "Display Error:", e.message()
    print out#(args)
    
def false_(*args):
    """
    function creating an array of false
    """
    if len(args) == 1:
        args += args
    return np.zeros(args,dtype=bool,order="F")

def find_(a,n=None,d=None,nargout=1):
    """
    function which returns a column vector or a tuple of column vectors (depending on the value of
    nargout) containing indices where a[indices] are true in the vector/matrix 'a'
    Argument 'n' tells us how many of those indices we want
    """
    
    if d:
        raise NotImplementedError

    elif nargout == 1:
        i = np.flatnonzero(a).reshape(-1,1)
        if isempty_(i):
            return []
                                                
        if n is not None:
            i = i.take(range(n))
        return np.array(i)
    elif nargout == 2:
        i,j = np.nonzero(a)
        if n is not None:
            i = i.take(n)
            j = j.take(n)

        return np.array(i).reshape(-1,1), np.array(j).reshape(-1,1)
    
    else:
        raise NotImplementedError

def floor_(a):
    """
    function floor
    """
    return np.floor(a)

def fopen_(*args):
    """
    function which opens a file
    """
    try:
        fp = open(*args)
        assert fp != -1
        return fp
    except:
        return -1
                                
def fclose_(*args):
    """
    function which closes a file
    """
    if args[0] == -1:
        print "No file to close, continue..."
    else:
        args[0].close()

#Unused functions
#def fullfile_(*args):
#    return os.path.join(*args)
#
#def intersect_(a,b,nargout=1):
#    if nargout == 1:
#        c = sorted(set(a) & set(b))
#        if isinstance(a,str):
#            return "".join(c)
#        elif isinstance(a,list):
#            return c
#        else:
#            # FIXME: the result is a column vector if
#            # both args are column vectors; otherwise row vector
#            return np.array(c)
#    raise NotImplementedError


def isempty_(a):
    """
    function which returns true if the argument 'a' is an empty array (array([])) or an empty list []
    """
    if a==[] :
        return True
    try:
        return 0 in a.shape
    except AttributeError:
        return False


def isequal_(a,b):
    """
    function which tests if 2 arrays have the same values
    """
    return np.array_equal(a,
                          b)

                          
def isscalar_(a):
    """np.isscalar returns True if a.__class__ is a scalar
    type (i.e., int, and also immutable containers str and
    tuple, but not list.) Our requirements are different"""
    try:
        return a.size == 1
    except AttributeError:
        return np.isscalar(a)

def length_(a):
    """
    function which returns the length of an array according to matlab definition
    """
    return max(np.asarray(a).shape)

#Not sure what this it for
#try:
#    def load_(a):
#        return loadmat(a) # FIXME
#except:
#    pass

def max_(a, d=None, nargout=None):
    """min_ and max_ function normally returns the same as matlab min and max in the following cases :
        min_(a,b) where a,b are arrays containing integer, inf or nan
        min_(a) and min_(a, nargout=2) where a is an array containing numbers, inf or nan
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
        min_(a,b) where a,b are arrays containing integer, inf or nan
        min_(a) and min_(a, nargout=2) where a is an array containing numbers, inf or nan
        (same for max_ obviously)
    """
    
    
    #The warnigns happens when NaNs are involved, but the function returns in any case what we want, so no need to print the warnings or to worry about them.
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

def mod_(a,b):
    """
    function which returns a modulo b
    """
    try:
        return a % b
    except ZeroDivisionError:
        return a
        
def ndims_(a):
    """
    function which returns the number of dimension of a
    """
    return a.ndim

def numel_(a):
    """
    function which returns the number of elements of a
    """
    return a.size

def rand_(*args,**kwargs):
    """
    function which returns a matrix of random numbers, according to matlab definition
    """
    if not args:
        return np.random.rand()
    if len(args) == 1:
        args += args
    try:
        return np.random.rand(np.prod(args)).reshape(args,order="F")
    except:
        pass

#Unused functions
#def ravel_(a):
#    return a.reshape(-1,1)
#
#def round_(a):
#    return np.round(a)
#
#def rows_(a):
#    return a.shape[0]
#
#def columns_(a):
#    return a.shape[1]

def size_(a, b=0, nargout=1):
    """
    function which returns the size of the matrix a, according to matlab defintion
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

def sum_(a, dim=None):
    """
    function which sums the elements of a in the corresponding dimension
    """
    if dim is None:
        return a.sum(0)
    else:
        return a.sum(dim-1)

def true_(*args):
    """
    function which creates an array of True
    """
    if len(args) == 1:
        args += args
    return np.ones(args,dtype=bool,order="F")


def zeros_(*args,**kwargs):
    """
    function which creates an array of zeros
    """
    if not args:
        return 0.0
    if len(args) == 1:
        if numel_(args[0])==1:
            args += args
        else:
            return np.zeros(args[0],**kwargs)
    return np.zeros(args,**kwargs)


def ones_(*args,**kwargs):
    """
    function which creates an array of ones
    """
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
# The linear algebra functions return the arguments like in matlab. Therefore we do modifications on the outputs of classical python functions.


def eig_(A, nargout=1):
    """
    function returning the eigenvalues (and eigenvectors if nargout=2) of A
    """
    if nargout == 1:
        return np.linalg.eigvals(A)
    else:
        D,V = np.linalg.eig(A) 
        return V,np.diag(D)
            
def cond_(A):
    """
    function returning the condition number of A
    """
    return np.linalg.cond(A)
    
def svd_(A, full_matrices=1, *args,**kwargs):
    """
    function returning the SVD decomposition of A,
    such that U*S*V.T = A
    """
    u,s,v= np.linalg.svd(A, full_matrices)
    return u, np.diag(s), v.T
    
def chol_(A,nargout=2):
    """
    function returning the cholesky factorization of A
    """
    try: #there is in python a exception when the matrix is not positive definite, which does not occur on matlab
        R = np.linalg.cholesky(A).T
        p = 0
    except:
        R = np.array([])
        p = 1
    return R,p   
    
def qr_(A,nargout=2):
    """
    function returning the qr factorization of A
    """
    return np.linalg.qr(A)
    
def inv_(A):
    """
    function returning the inverse of A
    """
    return np.linalg.inv(A)
     
def norm_(A, order=2, axis=None):
    """
    function returning the norm of the matrix or vector A
    """
    if isempty_(A):
        return 0
        
#for some reason, python does not give exactly (precision differs
#after the 10nth variable or so, but still this makes differences...) the same 2-norm
#when it does it from the norm of an array such as array([a b c...]) or such as array([[a b c...]])/array([[a] [b] [c]...])
#And it is closer to matlab when it is like array([[a b c..]]).
#Also, python will compute a different 'inf' norm for a vector such as  array([[a] [b] [c]...])/array([a,b,c]) or such as array([[a b c...]])
#and it is array([[a] [b] [c]...])/array([a,b,c]) which corresponds to the normal vector-inf norm we are looking for generally.
#Hence the following modifications :
    if isvector_or_scalar_(A):
        if order==2:
            #A=A.reshape(-1,1)# <- this also works
            A=A.reshape(1,-1)  
        if order==np.inf:
            A=A.reshape(-1,1)# <- this also works
#            A=A.reshape(-1)
        
    if axis==None:
        return np.linalg.norm(A, order,axis)
    else:
        return np.linalg.norm(A, order,axis-1)

    
def pinv_(A):
    """
    function returning the pseudo inverse of A
    """
    return np.linalg.pinv(A)
    
def solve_(A,b):
    """
    Returns x in problem A.dot(x)=b
    """
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
    """
    A one dimension polynomial class
    """
    return np.poly1d(A,r) 
    
def eye_(n):
    """
    Returns the identity matrix
    """
    return np.eye(n)
    
def concatenate_(arrs, axis=0):
    """
    Concatenate the arrays horizontally if axis=0, vertically if axis=1
    """
    copy_arrs=copy(arrs)
    for arr in arrs:
        if isempty_(arr):
            copy_arrs.remove(arr)
    if copy_arrs==[]:
        return np.array([])
    else:
        if axis == 0:
            return np.vstack(copy_arrs)
        else:
            return np.hstack(copy_arrs)

def sign_(A):
    """
    Returns a matrix of the same shape of A, with values -1 if the sign is negative, 0 if null, 1 if positive
    """
    return np.sign(A)
    
def setdiff_(a,b):
    """
    Returns the sorted, unique values in a that are not in b
    """
    return np.setdiff1d(a,b)
    
def isfield_(obj, name):
    """
    Returns true if 'name' is a key for obj
    """
    return obj.__dict__.has_key(str(name))
    
def any_(A):
    """
    Returns true if there is any boolean True in A
    """
    return A.any()
    
def sqrt_(x):
    """
    Returns the square root of x
    """
    return np.sqrt(x)

def real_(x):
    """
    Returns the real part of x
    """
    ret = np.real(x)
    #print "type np real ret", type(ret)
    return ret

def imag_(x):
    """
    Returns the imaginary part of x
    """
    ret = np.imag(x)
    return ret

def exp_(x):
    """
    Returns the exponential of x
    """
    return np.exp(x)
    
def null_(A, eps=1e-15):
    """
    Returns the null space of A
    """
    u, s, vh = np.linalg.svd(A)
    padding = max(0,np.shape(A)[1]-np.shape(s)[0])
    null_mask = np.concatenate(((s <= eps), np.ones((padding,),dtype=bool)),axis=0)
    null_space = scipy.compress(null_mask, vh, axis=0)
    return scipy.transpose(null_space)        

#Unused function
#def tf_mapper(x):
#    #print "x", x
#    if x:
#        return 1
#    else:
#        return 0
                
def logical_or_(a,b):
    """
    Returns an array corresponding to element-wise operation 'or' between a and b
    """
#    vtf_mapper = np.vectorize (tf_mapper)
#    ret = np.logical_or(a,b)                
#    return vtf_mapper(ret)
    return np.logical_or(a,b)
    
def logical_and_(a,b):
    """
    Returns an array corresponding to element-wise operation 'and' between a and b
    """
    return np.logical_and(a,b)
    
def strcat_(*args):
    """
    Concatenates some arguments
    """
    ret = ""
    for arg in args:
        ret = ret + str(arg)
        
    #print "strcat_:", ret
    return ret
    
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
    
    
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()