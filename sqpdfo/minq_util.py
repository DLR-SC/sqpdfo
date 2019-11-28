"""
Supplementary or Auxiliary Functions for Snobfit and Minq

Most of them are the functions of the Rough MATLAB-Numpy Equivalents.
"""

import numpy

# Some constants
eps   = numpy.finfo('d').eps
inf   = numpy.inf
nan   = numpy.nan

# Some internal functions of numpy
rand  = numpy.random.rand
norm  = numpy.linalg.norm
inv   = numpy.linalg.inv
isnan = numpy.isnan


# ----------------------------------------------------------
def within(low, x, high):
    return numpy.logical_and(low <= x, high >= x)


def find(*args,**kw):
    tup = numpy.where(*args,**kw)
    if len(tup)==0:
        return numpy.array([],'int')
    else:
        return tup[0]


def toCol( v ):
    """
    v is a vector
    """
    return numpy.matrix( v.reshape(v.size,1) )


def toRow( m ):
    """
    m is a  n x 1 matrix
    """
    n = m.shape[0]
    return numpy.array( numpy.asarray(m).reshape(1,n)[0] )


def mldiv(a,b):
    """ Backslash or left matrix divide.

    a\b is the matrix division of A into B, which is roughly the
    same as INV(A)*B
    """
    if a.shape[0] == a.shape[1]:
        return numpy.linalg.solve(a,b)
    else:
        return numpy.asarray( numpy.linalg.lstsq(a,b)[0] ).T[0]


def svd( a ):
    """
    a is a mat
    """
    (U,S,V) = numpy.linalg.svd(a,0)
    return (U,S,V.T)


def qr(a):
    return numpy.linalg.qr(a, mode='r')


def vector(n):
    """
    Return a vector of the given length.
    """
    return numpy.empty(n,'d')

def ivector(n=0):
    """
    Return a int vector of the given length.
    """
    return numpy.zeros(n,'int')


def sort(x):
    ind = numpy.argsort(x);
    return (x[ind], ind)
    

def std(x):
    """
    STD(X) normalizes by (N-1) where N is the sequence length.
    This makes STD(X).^2 the best unbiased estimate of the variance
    if X is a sample from a normal distribution.
    """
    return numpy.std(x,ddof=1)


def max_(x):
    """
    [Y,I] = MAX(X) returns the indices of the maximum values in vector I.
    If the values along the first non-singleton dimension contain more
    than one maximal element, the index of the first one is returned.
    """
    if len(x.shape) > 1:
        print( "It's a Matrix")

    idx = numpy.argmax(x)
    return ( x[idx], idx )


def min_(x):
    """
    [Y,I] = MIN(X) returns the index of the minimum value in vector I.
    If the values along the first non-singleton dimension contain more
    than one minimal element, the index of the first one is returned.
    """
    if len(x.shape) > 1:
        print( "It's a Matrix")
    idx = numpy.argmin(x)
    return ( x[idx], idx )



def isEmpty(x):
    if len(x)==0: return True
    if x.shape[0]==0: return True
    if x.shape[1]==0: return True
    return False


def removeByInd(a, ind):
    n = len(a) - len(ind)
    ret = numpy.zeros( n, 'int' )
    k = 0
    for i in xrange(len(a)):
        if (i not in ind):
            ret[k] = a[i]
            k += 1
    return ret


def dot3(a,b,c):
    # Calculate a*(b*c)
    return numpy.dot(a, numpy.dot(b,c) )


def crdot(m, v):
    n = len(v)
    _sum = 0.0
    for i in range(n):
        _sum += m[i,0]*v[i]
        
    return _sum

#---------------------------------------------------------------
def rsort(x,w=None):
    """
    Sort x in increasing order, remove multiple entries,
    and adapt weights w accordingly x and w must both be a row or a column
    default input weights are w=1
 
    If w==None, the weighted empirical cdf is computed at x
    dof = len(x) at input is the original number of degrees of freedom

    Warning: when you use this function, make sure x and w is row vector
    """
    if  w is None:
        w = numpy.ones( x.shape )
       
    ind = numpy.argsort(x)
    x   = x[ind]
    w   = w[ind]

    # Remove repetitions
    n = len(x)

    xnew = numpy.append(  x[1:n], inf )

    ind = find( xnew != x ) 
    nn  = len(ind)

    x=x[ind]

    ww = numpy.zeros(nn)
    ww[0] = sum( w[0:ind[0]+1 ] ) 
    for i in xrange( nn-1 ):
        ww[i+1] = sum( w[ ind[i]:ind[i+1]] )
    
    # Get cumulative sum of weights
    cdfx = numpy.cumsum(ww)

    # Adjust for jumps and normalize
    cdfx = (cdfx-0.5*ww)/cdfx[nn-1]
    dof  = n
    
    return (x,ww,cdfx,dof)


    
# ---------------------------------------------------------------
if __name__ == '__main__':
   print( rand(5) )
