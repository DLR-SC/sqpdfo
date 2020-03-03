# -*- coding: utf-8 -*-

"""

This file contains functions which serve mostly 2 purposes :
    -simply providing a shortcut for matlab functions, for instance real(A) in matlab becomes real_(A) in python, and
    real_(A) actually calls np.real(A)
    -imitating at best the output of matlab functions. For instance, U,S,V=svd function in matlab and python does not give exactly
    outputs with the same shapes : S is a vector in python and a diagonal matrix in matlab, and V.transpose in python = V in matlab.

"""
import warnings
from copy import copy
import scipy
import numpy as np
import sys
#warnings.simplefilter('error', RuntimeWarning)

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


def disp_(*args):
    print(strcat_(*args))
    

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
            i = i.take(list(range(n)))
        return np.array(i)
    elif nargout == 2:
        i,j = np.nonzero(a)
        if n is not None:
            i = i.take(n)
            j = j.take(n)

        return np.array(i).reshape(-1,1), np.array(j).reshape(-1,1)
    
    else:
        raise NotImplementedError


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
        print("No file to close, continue...")
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
    if np.array(a).size == 0 :
        return True
    try:
        return 0 in a.shape
    except AttributeError:
        return False


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
    warnings.simplefilter('ignore', RuntimeWarning)
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
        warnings.simplefilter('default', RuntimeWarning)
        return ret, ret2
    else:
        warnings.simplefilter('default', RuntimeWarning)
        return ret                                


def numel_(a):
    """
    function which returns the number of elements of a
    """
    return a.size



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
    assert (dim is None)

    return a.sum(0)



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
    # return np.linalg.cond(A)
    # we do not use the function cond() from numpy.linalg
    # because it has no safeguard in case the smallest eigenvalue is zero!

    s = np.linalg.svd(A, compute_uv=False)
    if s[...,-1] == 0.0:
        return np.infty
    else:
        return s[..., 0] / s[..., -1]
    
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
    

def fprintf_(fid,*args, **kwargs):
    """This is supposed to work approximately like the matlab fprintf function in the following case in the cases where
    fid=1 (returns output on the screen), fid=2 (returns output in standard error ouput) or fid is a file identifier"""
    out = ""
    for arg in args:
        try:
            out = out + str(arg)
        except Exception as e:
            print("Display Error:", e.message())

    if isinstance(fid,int):
        if  fid == 1:
            print(out, end=' ')
        elif fid == 2:
            print(out, end=' ', file=sys.stderr)
        else:
            print('### Warning from runtime.py: unknown integer file identifier (1 and 2 are known) ###')
    elif isinstance(fid,str):
        print(fid + out, end=' ')
    elif hasattr(fid, "read"):
        # writes to file
        print(out, end=' ', file=fid)
    else:
        print(fid + out, end=' ')
        print('### Warning from runtime.py: unknown file identifier type (int, str and file are known) ###')
        
    
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
    for i in range(0,len(copy_arrs)):
        j=0
        for arr in copy_arrs:
            if isempty_(arr):
                #copy_arrs.remove(arr)
                # If warnings are 'on', list.remove()
                # is stated as deprecated when running
                del copy_arrs[j]
                break
            j=j+1
    #print(copy_arrs)
    #print(isinstance(copy_arrs,np.ndarray))
    #print(np.array(copy_arrs))
    #print(0 in copy_arrs.shape)
    #if isinstance(copy_arrs,np.ndarray) and copy_arrs.size == 0:
    #    return np.array([])
    #elif not isinstance(copy_arrs,np.ndarray) and 0 in np.array(copy_arrs).shape:
    #    return np.array([])
    if copy_arrs == []:
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
    return str(name) in obj.__dict__
    
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
    null_space = np.compress(null_mask, vh, axis=0)
    return np.transpose(null_space)        

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
    print("Warning: randn does not do anything, because it's discouraged syntax")
    
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
