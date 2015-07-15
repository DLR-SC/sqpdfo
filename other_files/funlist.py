#!/usr/local/bin/python

from numpy import *

def Rosenbrock2D(x):
    f = 100*(x[1]-x[0]**2)**2 + (1-x[0])**2
    return f
def beale(x):
    return (1.5-x[0]*(1-x[1]))**2+(2.25-x[0]*(1-x[1]**2))**2+(2.625-x[0]*(1-x[1]**3))**2
def rosenbrock(x):
    n = size(x,0)
    f = 0
    for i in range(0,n-1):
        f = f + (1-x[i])**2 + 100*(x[i+1]-x[i]**2)**2
    return f
def book1(x):
        return (6*x[0]-2)**2*sin(12*x[0]-4)
def book2(x):
        return (x[1] - (5.1*x[1]/(4*pi**2)) + 5*x[0]/pi - 6)**2 + 10* ((1-1/(8*pi))*cos(x[0])+1) + 5*x[0]
def quadratic(x):
        n = size(x,0)
        f = 0
        for i in range(0,n):
            f = f + x[i]**2
        return f
def cubic(x):
        n = size(x,0)
        f = 0
        for i in range(0,n):
            f = f + x[i]**3
        return f
def poly4(x):
        H = array([[4,-1],[-1,4]])
        b = array([[1],[1]])
        c = [1,1]
        xHx = dot(x.transpose(),dot(H,x))
        f = xHx[0]/2 + dot(b.transpose(),x)[0]
        f = f + x[0]**3
        return f
def Aquadratic(x):
        n = size(x,0)
        H = -ones((n,n), float)
        for i in range(n):
            H[i,i] = 4
        b = ones((n,1), float)
        xHx = dot(x.transpose(),dot(H,x))
        return xHx[0]/2 + dot(b.transpose(),x)[0]
def sin1(x):
        return sin(x[0]/3) + cos(x[1]/5)
def sin2(x):
        return 2*x[1]* sin(x[0]+x[1]) + x[0]* sin(x[0]+x[1]) + 3*x[1]**2 + x[0]**2
def sin4(x):
        f = x[0]**2 + x[1]**2
        if x[0] < 1.1 and x[0] > 0.9 and x[1] < 1.1 and x[1] > 0.9:
            f = f-e**abs(x[0]-1)
        if x[0] < 0.1 and x[0] > -0.1 and x[1] < -1 and x[1] > -0.7:
            f = f-e**abs(x[0])
        return  f
def sin3(x):
        return 2*(x[1]**3)* sin((x[0])*pi) + x[0]* sin((x[0]+x[1])*10) 
def plateau(x):
        fY = 0
        if x[0] < -3:
            if x[1] < -3:
                fY = -5
            elif x[1] < 0:
                fY = -3
            else:
                fY = 5
        elif x[0] < 0:
            if x[1] < 0:
                fY = -3
            else:
                fY = 5
        else:
            if x[1] < 0:
                fY = 5
            else:
                fY = 0
        return fY 
def plateau2(x):
        fY = 0
        if x[0] < -3:
            if x[1] < -3:
                fY = obj(x,'sin1', RS)
            elif x[1] < 0:
                fY = obj(x,'sin3', RS)
            else:
                fY = obj(x,'sin2', RS)
        elif x[0] < 0:
            if x[1] < 0:
                fY = obj(x,'sin3', RS)
            else:
                fY = obj(x,'sin2', RS)
        else:
            if x[1] < 0:
                fY = obj(x,'sin2', RS)
            else:
                fY = obj(x,'rosenbrock', RS)
        return fY
def linquad(x):
        n = size(x,0)
        sum = 0
        for i in range(1,n):
            sum = sum + x[i]
        return x[0]**2 + sum
def constant(x):
        return 1.
def svrQP(x):
        return RS.svrQP(x,0)
def obj(x, f_name):

    if ( f_name == 1 ) :                 # Rosenbrock                 x0 = [ -1.2 1]
        f  = 100 * ( x[1] - x[0]**2 ) **2 + (1-x[0])**2
        #   c  = [ 10 * ( x[1] - x[0]**2 )  (1-x[0]) ]
        #   f = norm( c )**2
        #   f = c'*c
    elif( f_name == '43bis' ) :       # Penalty I      x0 = [ 1 2 3 4 ]
        n = size(x,0)
        f = 0
        sum1 = 0
        for i in range(n):
            sum1 = sum1 + x[i]**2
            f = f + sqrt(1.e-5)*(x[i]-1)
        f = f + sum1-0.25

    elif ( f_name == 2 ) :             # Beale                      x0 = [ 1.0 1.0]
        f  = (1.5-x[0]*(1-x[1]))**2+(2.25-x[0]*(1-x[1]**2))**2+(2.625-x[0]*(1-x[1]**3))**2
        
    elif ( f_name == 3 ) :             # Powell                  x0 = [ 1.0 1.0 1.0 1.0]
        f1 = x[0]+10*(x[1])
        f2 = sqrt(5)*(x[2]-x[3])
        f3 = (x[1]-2*(x[2]))**2
        f4 = sqrt(10)*(x[0]-x[3])**2
        f  = f1**2 + f2**2 + f3**2 + f4**2
        
    elif ( f_name == 4 ) :             # Freudenstein and Roth      x0 = [ 0.5 -2 0.5 -1 ]
        f = ( -13 + x[0] + ( ( 5 - x[1] )*x[1]- 2)*x[1] )**2 + ( -29 + x[0] + (( x[1] + 1)*x[1] - 14)*x[1])**2
        
    elif ( f_name == 5 ) :             # Trigonometric (n=3)        x0 = [ 0.2 0.2 0.2 ]
        f1 = (3 - cos(x[0])-cos(x[1])-cos(x[2]) + 1 - cos(x[0]) - sin(x[0]))**2
        f2 = (3 - cos(x[0])-cos(x[1])-cos(x[2]) + 2*(1 - cos(x[1])) - sin(x[1]))**2
        f3 = (3 - cos(x[0])-cos(x[1])-cos(x[2]) + 3*(1 - cos(x[2])) - sin(x[2]))**2
        f  = f1**2 + f2**2 + f3**2
        
    elif ( f_name == 6 ) :             # White and Holst            x0 = [ -1.2 1]
        f  = 100 * ( x[1] - x[0]**3 ) **2 + (1-x[0])**2
        
    elif ( f_name == 7 ) :             # Penalty (n=5)              x0 = [ 1 2 3 4 5 ]
        f = (x[0]-1)**2+(x[1]-1)**2+(x[2]-1)**2+(x[3]-1)**2+(x[0]**2+x[1]**2+x[2]**2+x[3]**2+x[4]**2-0.25)**2
        
    elif ( f_name == 8 ) :             # Hager (n=5)                x0 = [ 1 2 3 4 5 ]
        f = exp(x[0])-x[0] + exp(x[1])-sqrt(2)*x[1] + exp(x[2])-sqrt(3)*x[2] + exp(x[3])-sqrt(4)*x[3] + exp(x[4])-sqrt(5)*x[4]
        
    elif ( f_name == 9 ) :             # Himmelblau (n=2)           x0 = [ 1 1]
        f = (x[0]**2-x[1]-11)**2 + (x[0]+x[1]**2-7)**2
        
    elif ( f_name == 10 ) :            # Diagonal 5 ( n = 6 )       x0 = [ 1:10 ]
        fv = log( exp( x ) + exp( -x ) )
        f  = dot(fv.transpose(),fv)
        f = f[0]
        
    elif ( f_name == 11 ) :            # Quadratic (n = 5)          x0 = [5 -5 7 -7 15]
        A = array([
            [5.9274,    3.2013,    3.9526,    5.1713,    5.4177],
            [3.2013,    3.5665,    4.1836,    4.9734,    3.7165],
            [3.9526,    4.1836,    5.8019,    6.5200,    4.3957],
            [5.1713,    4.9734,    6.5200,    8.2392,    5.7992],
            [5.4177,    3.7165,    4.3957,    5.7992,    5.3939] ])
        f = dot(x.transpose(),dot(A,x))
        f = f[0]
    elif ( f_name == 12 ) :           # Cosh ( n = 5 )              x0 = [ 1 1 1 1 1 ]
        A = array([
            [1.14029,   0.97433,   1.59707,   0.78051,   1.05379],
            [0.97433,   0.45448,   1.51828,   1.41074,   0.55527],
            [1.59707,   1.51828,   1.60261,   0.63912,   0.73487],
            [0.78051,   1.41074,   0.63912,   1.96847,   0.68060],
            [1.05379,   0.55527,   0.73487,   0.68060,   0.27444] ])
        f = cosh(dot(x.transpose(),dot(A,x)) )
        f = f[0]
        
    elif ( f_name == 13 ) :                 # Sisser
        f = 3*x[0]**4-2*x[0]**2*x[1]**2+3*x[1]**4
        
    elif ( f_name == 14 ) :                 # SQRTMX
        n = int(sqrt(size(x)))
        X = reshape(x,(n,n))
        f = []
        for i in range(n):
            for j in range(n):
                p = 0
                for k in range(n):
                    p = p + X[i,k]*X[k,j]
                f = concatenate(( f, [p-sin(i+1)**2-sin(j+1)**2]), 0)
        
    elif ( f_name == 15 ) :                 # NONMSQRT
        n = int(sqrt(size(x)))
        X = reshape(x,(n,n))
        f = []
        for i in range(n):
            for j in range(n):
                p = 0
                for k in range(n):
                    p = p + X[i,k]*X[i,j]
                f = concatenate(( f, [p-sin(i+1)**2-sin(j+1)**2 ]), 0)
        
    elif ( f_name == 16 ) :                 # Mancino
        n = size(x)
        f = []
        for i in range(n):
            sumhij = 0
            for  j in range(n):
                if ( i != j ) :
                    vij     = sqrt( x[j]**2 + (float(i+1)/float(j+1)))
                    sij     = sin(log(vij))
                    cij     = cos(log(vij))
                    sumhij  = sumhij + vij*(sij**5+cij**5)
            f = concatenate(( f, sumhij+14*n*x[i]+(i+1-(n/2.))**3 ), 0)
        
    elif ( f_name == 17 ) :                 # Engval1
        n = size(x)
        f = []
        for i in range(n-1):
            f = concatenate(( f, (x[i]**2+x[i+1]**2)**2-4*x[i]+3 ), 0)
        
    elif ( f_name == 18 ) :                 # Engval2
        f1 = x[0]**2+x[1]**2+x[2]**2-1
        f2 = x[0]**2+x[1]**2+(x[2]-2)**2-1
        f3 = x[0]+x[1]+x[2]-1
        f4 = x[0]+x[1]-x[2]+1
        f5 = x[0]**3+3*x[1]**2+(5*x[2]-x[0]+1)**2-36
        f = [ f1[0], f2[0], f3[0], f4[0], f5[0] ]
        
    elif ( f_name == 19 ) :                 # BRKMCC
        p = -0.25*x[0]**2-x[1]**2+1
        h = x[0]-2*x[1]+1
        f = (x[0]-2)**2 + (x[1]-1)**2 + 1/(25*p)+5*h**2
        
    elif ( f_name == 20 ) :                 # Booth
        f = (x[0]+2*x[1]-7)**2 + (2*x[0]+x[1]-5)**2
        
        
    elif ( f_name == 21 ) :                 # Rosenbrock                 x0 = [ -1.2 1]
        f1 = 10 * ( x[1] - x[0]**2 )
        f2 = (1-x[0])
        f  = [ f1[0], f2[0] ]
        
    elif ( f_name == 22 ) :                 # Beale                      x0 = [ 1.0 1.0]
        f1 = 1.5-x[0]*(1-x[1])
        f2 = 2.25-x[0]*(1-x[1]**2)
        f3 = 2.625-x[0]*(1-x[1]**3)
        f  = [ f1[0], f2[0], f3[0] ]
        
    elif ( f_name == 23 ) :                   # Powell                  x0 = [ 1.0 1.0 1.0 1.0]
        f1 = x[0]+10*(x[1])
        f2 = sqrt(5)*(x[2]-x[3])
        f3 = (x[1]-2*(x[2]))**2
        f4 = sqrt(10)*(x[0]-x[3])**2
        f  = [ f1[0],  f2[0], f3[0], f4[0]]
        
    elif ( f_name == 24 ) :             # Freudenstein and Roth      x0 = [ 0.5 -2 0.5 -1 ]
        f1 = -13 + x[0] + ( ( 5 - x[1] )*x[1]- 2)*x[1]
        f2 = -29 + x[0] + (( x[1] + 1)*x[1] - 14)*x[1]
        f = [ f1[0],  f2[0] ]
        
    elif ( f_name == 25 ) :             # Trigonometric (n=3)        x0 = [ 0.2 0.2 0.2 ]
        f1 = (3 - cos(x[0])-cos(x[1])-cos(x[2]) + 1 - cos(x[0]) - sin(x[0]))**2
        f2 = (3 - cos(x[0])-cos(x[1])-cos(x[2]) + 2*(1 - cos(x[1])) - sin(x[1]))**2
        f3 = (3 - cos(x[0])-cos(x[1])-cos(x[2]) + 3*(1 - cos(x[2])) - sin(x[2]))**2
        f  = [ f1[0], f2[0], f3[0] ]
        
    elif ( f_name == 26 ) :             # White and Holst            x0 = [ -1.2 1]
        f1 = 10 * ( x[1] - x[0]**3 )
        f2 = 1-x[0]
        f  = [ f1[0], f2[0] ]
        
    elif ( f_name == 27 ) :             # Penalty (n=5)              x0 = [ 1 2 3 4 5 ]
        f1 = x[0]-1
        f2 = x[1]-1
        f3 = x[2]-1
        f4 = x[3]-1
        f5 = x[0]**2+x[1]**2+x[2]**2+x[3]**2+x[4]**2-0.25
        f = [ f1[0], f2[0], f3[0], f4[0], f5[0] ]
        
    elif ( f_name == 28 ) :             # Powell badly scaled        x0 = [ 0 1 ]
        f1 = 10000.*x[0]*x[1]-1
        f2 = exp(-x[0])+exp(-x[1])-1.0001
        f = [ f1[0], f2[0] ]
        
    elif ( f_name == 29 ) :             # Jenrich Sampson            x0 = [ 0.3 0.4 ]
        m = 10
        f = []
        for i in range(m):
            f = concatenate(( f, 2+2*(i+1)-( exp((i+1)*x[0])+exp((i+1)*x[1])) ), 0)
        
    elif ( f_name == 30 ) :             # Helical valley             x0 = [-1 0 0] 
        if ( x[0] > 0 ):
            theta = math.atan( x[1]/x[0] ) / ( 2*pi)
        elif ( x[0] < 0 ):
            theta = math.atan( x[1]/x[0] ) / (2*pi) + 0.5
        else:
            theta = 1/(2*pi)
        f1 = 10*(x[2]-10*theta)
        f2 = 10*(sqrt(x[0]**2+x[1]**2)-1)
        f = [ f1[0], f2[0] ]
        
    elif ( f_name == 31 ) :             # Bard                       x0 = [1 1 1]
        m = 15
        y = [ 0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58, 0.73, 0.96, 1.34, 2.10, 4.39]
        f = []
        for i in range(m) :
            vi = 16-(i+1)
            wi = min( i+1, vi )
            f = concatenate(( f, y[i]-(x[0]+float(i+1)/(vi*x[1]+wi*x[2])) ), 0)
        
    elif ( f_name == 32 ) :             # Gaussian                   x0 = [0.4 1 0]
        m = 15
        y = [ 0.0009, 0.0044, 0.0175, 0.0540, 0.1295, 0.2420, 0.3521, 0.3989, 0.3521, 0.2420, 0.1295, 0.0540, 0.0175, 0.0044, 0.0009 ]
        f = []
        for i in range(m) :
            ti = 0.5*(8-(i+1))
            f = concatenate(( f, x[0]*exp(-0.5*x[1]*(ti-x[2])**2)-y[i] ), 0)
        
    elif( f_name == 33 ) :              # Meyer                     x0 = [0.02 4000 250]
        m = 16
        y = [ 34780, 28610, 23650, 19630, 16370, 13720, 11540, 9744, 8261, 7030, 6005, 5147, 4427, 3820, 3307, 2872 ]
        f = []
        for i in range(m) :
            ti = 45 + 5*(i+1)
            f = concatenate(( f, x[0]*exp(x[1]/(ti+x[2]))-y[i] ), 0)
        
    elif( f_name == 34 ) :              # Gulf                     x0 = [ 5 2.5 0.15 ]
        m = 99
        f = []
        for i in range(m) :
            ti = 0.01 * (i+1)
            yi = 25. + (-50.*log(ti))**(2./3.)
            f = concatenate(( f, exp(-(abs(yi-x[1])**(x[2]))/x[0])-ti ), 0)
        
    elif( f_name == 35 ) :              # Box                      x0 = [ 0 10 20 ]
        m = 10
        f = []
        for i in range(m) :
            ti = 0.1 * (i+1)
            f = concatenate(( f, exp(-ti*x[0])-exp(-ti*x[1])-x[2]*(exp(-ti)-exp(-10*ti)) ), 0)
   
    elif( f_name == 36 ) :              # Wood                     x0 = [ -3 -1 -3 -1 ]
        f1 = 10*(x[1]-x[0]**2)
        f2 = 1-x[0]
        f3 = sqrt(90)*(x[3]-x[2]**2)
        f4 = 1-x[2]
        f5 = sqrt(10)*(x[1]+x[3]-2)
        f6 = (x[1]-x[3])/sqrt(10)
        f = [ f1[0], f2[0], f3[0], f4[0], f5[0], f6[0]]
        
    elif( f_name == 37 ) :              # Kowalik Osborne          x0 = [ 0.25 0.39 0.415 0.39 ]
        m = 11
        y = [0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627, 0.0456, 0.0342, 0.0323, 0.0235, 0.0246]
        u = [4,      2,      1,      0.5,    0.25,   0.1670, 0.1250, 0.1,    0.0833, 0.0714, 0.0625]
        f = []
        for i in range(m) :
            f = concatenate(( f, y[i]-x[0]*(u[i]**2+u[i]*x[1])/(u[i]**2+u[i]*x[2]+x[3])), 0)
        
    elif( f_name == 38 ) :              # Brown Dennis          x0 = [ 25 5 -5 -1 ]
        m = 20
        f = []
        for i in range(m) :
            ti = 0.2*(i+1)
            f = concatenate(( f, (x[0]+ti*x[1]-exp(ti))**2+(x[2]+x[3]*sin(ti)-cos(ti))**2), 0)
        
    elif( f_name == 39 ) :              # Osborne1              x0 = [ 0.5 1.5 -1 0.01 0.02 ]
        m = 33
        y = [0.844, 0.908, 0.932, 0.936, 0.925, 0.908, 0.881, 0.850, 0.818, 0.784, 0.751, 0.718, 0.685, 0.658, 0.628, 0.603, 0.580, 0.558, 0.538, 0.522, 0.506, 0.490, 0.478, 0.467, 0.457, 0.448, 0.438, 0.431, 0.424, 0.420, 0.414, 0.411, 0.406]
        f = []
        for i in range(m) :
            ti = 10*(i)
            f = concatenate(( f, y[i]-(x[0]+x[1]*exp(-ti*x[3])+x[2]*exp(-ti*x[4])) ), 0)
        
    elif( f_name == 40 ) :              # Bigs EXP   6          x0 = [ 1 10 1 5 4 3 ]
        m = 13
        f = []
        for i in range(m) :
            ti = 0.1*(i+1)
            yi = exp(-ti)-5*exp(-10*ti)+3*exp(-4*ti)
            f = concatenate(( f, x[2]*exp(-ti*x[0])-x[3]*exp(-ti*x[1])+x[5]*exp(-ti*x[4])-yi ), 0)
        
    elif( f_name == 41 ) :              # Osborne 2             x0 = [ 1.3 0.65 0.65 0.7 0.6 3 5 7 2 4.5 5.5 ]
        m = 65
        y=[]
        y[1:15]  = [ 1.366, 1.191, 1.112, 1.013, 0.991, 0.885, 0.831, 0.847, 0.786, 0.725, 0.746, 0.679, 0.608, 0.655, 0.616 ]
        y[16:30] = [ 0.606, 0.602, 0.626, 0.651, 0.724, 0.649, 0.649, 0.694, 0.644, 0.624, 0.661, 0.612, 0.558, 0.633, 0.495 ]
        y[31:45] = [ 0.500, 0.423, 0.395, 0.375, 0.372, 0.391, 0.396, 0.405, 0.428, 0.429, 0.523, 0.562, 0.607, 0.653, 0.672 ]
        y[46:60] = [ 0.708, 0.633, 0.668, 0.645, 0.632, 0.591, 0.559, 0.597, 0.625, 0.739, 0.710, 0.729, 0.720, 0.636, 0.581 ]
        y[61:65] = [ 0.428, 0.292, 0.162, 0.098, 0.054 ]
        f = []
        for i in range(m) :
            ti = 0.1*(i)
            f = concatenate(( f, y[i]-(x[0]*exp(-ti*x[4])+x[1]*exp(-x[5]*(ti-x[8])**2)+x[2]*exp(-x[6]*(ti-x[9])**2)+x[3]*exp(-x[7]*(ti-x[10])**2)) ), 0)
        
    elif( f_name == 42 ) :              # Watson 6              x0 = [ 0 0 0 0 0 0 ]
        m = 31
        n = 6
        f = []
        for i in range(m-2) :
            ti = (i+1)/29
            sum1 = 0
            sum2 = 0
            for j in range(1,n) :
                sum1 = sum1 + (j) * x[j] * ti**(j-1)
                sum2 = sum2 +         x[j] * ti**(j-0)
            f = concatenate(( f, sum1-sum2**2-1), 0)
        f = concatenate(( f, x[0],  x[1]-x[0]**2-1), 0)
        
    elif( f_name == 43 ) :              # Penalty I             x0 = [ 1 2 3 4 ]
        n = 4
        m = n +1
        f = []
        sum1 = 0
        for i in range(n):
            sum1 = sum1 + x[i]**2
            f = concatenate(( f, sqrt(1.e-5)*(x[i]-1)), 0)
        f = concatenate(( f, sum1-0.25), 0)
        
    elif( f_name == 44 ) :              # Penalty II             x0 = [ 0.5 0.5 0.5 0.5 ]
        n = 4
        a = sqrt(1.e-5)
        f =  x[0]-0.2 
        sum1 = n*x[0]**2
        for i in range(1,n) :
            yi = exp((i+1)/10.) + exp((i)/10.)
            f = concatenate(( f, a*(exp(x[i]/10.)+exp(x[i-1]/10.)-yi)), 0)
            sum1 = sum1 + (n-i)*x[i]**2
        for i in range(n, 2*n-1) : #= n+1:2*n-1
            yi = exp((i+1)/10) + exp((i)/10)
            f = concatenate(( f, a*(exp(x[i-n+1]/10.)+exp(-1./10.)-yi)), 0)
        f = concatenate(( f, sum1-1 ), 0)
        
    elif( f_name == 45 ) :              # Vardim                 x0 = [ 0.8 0.6 0.4 0.2 0 ]
        n = 5
        m = n+2
        f = []
        sum = 0
        for i in range(n):
            f = concatenate(( f, x[i]-1), 0)
            sum = sum + (i+1)*(x[i]-1)
        f = concatenate(( f, sum,  sum**2), 0)
        
    elif( f_name == 46 ) :              # Trigonometric MGH      x0 = 0.125*ones(8,1)
        n = 8
        m = n
        sum = 0
        for i in range(n):
            sum = sum + cos(x[i])
        f = []
        for i in range(m) :
            f = concatenate(( f, n-sum+(i+1)*(1-cos(x[i]))-sin(x[i])), 0)
        
    elif( f_name == 47 ) :              # Brownal                x0 = 0.5*ones(5,1)
        n = 4
        m =n
        sum = 0
        prd = 1
        for i in range(n):
            sum = sum + x[i]
            prd = prd * x[i]
        f = []
        for i in range(m-1) :
            f = concatenate(( f, x[i]+sum-n-1), 0)
        f = concatenate(( f, prd-1), 0)
        
    elif( f_name == 48 ) : # Bdvalue  x0 = -[0.0826 0.1488 0.1983 0.2314 0.2479 0.2479 0.2314 0.1983 0.1488 0.0826]
        n = 10
        m = n
        h = 1/(n+1)
        f = [ 2*x[0]-x[1]+0.5*h**2*(x[0]+h+1)**3 ]
        for i in range(1,n-1) :
            f = concatenate(( f, 2*x[i]-x[i-1]-x[i+1]+0.5*h**2*(x[i]+i*h+1)**3), 0)
        f = concatenate(( f, 2*x[n-1]-x[n-2]+0.5*h**2*(x[n-1]+n*h+1)**3), 0)
        
    elif( f_name == 49 ) :  # Integral equation  x0 = -[0.0826 0.1488 0.1983 0.2314 0.2479 0.2479 0.2314 0.1983 0.1488 0.0826]
        n = 10
        m = n
        h = 1/(n+1)
        f = []
        for i in range(n):
            sum1= 0
            for j in range(i):
                sum1 = sum1 + (j+1)*h*(x[j]+(j+1)*h+1)**3
            sum2= 0
            for j in range(i,n):  #i+1:n
                sum2 = sum2 + (1-(j+1)*h)*(x[j]+(j+1)*h+1)**3
            f = concatenate(( f, x[i]+0.5*h*((1-(i+1)*h)*sum1+(i+1)*h*sum2)), 0)

    elif( f_name == 50) :              # Broyden tridiagonal  x0 = -ones( 10, 1) 
        n = 10
        m = 10
        f = [ (3-2*x[0])*x[0]-2*x[1]+1 ]
        for i in range(1, m-1): # = 2:m-1
            f = concatenate(( f, (3-2*x[i])*x[i]-x(i-1)-2*x(i+1)+1), 0)
        f = concatenate(( f, (3-2*x[n-1])*x[n-1]-x[n-2]+1), 0)
        
    elif( f_name == 51) :              # Broyden banded        x0 = -ones( 10, 1) 
        n = 10
        m = 10
        mi = 5
        mu = 1
        f = []
        for i in range(m) :
            sum = 0
            for j in range(max(1,i-mi),i-1):
                sum = sum + x[j]*(1+x[j])
            for j in range(i+1,min(i+mu,n)) :
                sum = sum + x[i+1]*(1+x[i+1])
            f = concatenate(( f, x[i]*(2+5*x[i]**2)+1-sum), 0)

    elif( f_name == 52) :              # Linear FR        x0 = ones( 8, 1) 
        n = 8
        m = 10
        sum = 0
        for j in range(n):
            sum = sum + x[j]
            f = []
            for i in range(n):
                f = concatenate(( f, x[i]-(2./float(m))*sum-1), 0)
            for i in range(n,m): 
                f = concatenate(( f, -(2./float(m))*sum-1), 0)
                
    elif( f_name == 53) :              # Linear R1-1        x0 = ones( 8, 1) 
        n = 8
        m = 10
        sum = 0
        for j in range(n):
            sum = sum + (j+1)*x[j]
            f = []
            for i in range(m) :
                f = concatenate(( f, (i+1)*sum-1), 0)
                
    elif( f_name == 54) :              # Linear R1-2        x0 = ones( 8, 1) 
        n = 8
        m = 10
        sum = 0
        for j in range(1,n-1):  #= 2:n-1
            sum = sum + (j+1)*x[j]
        f = [ -1]
        for i in range(1,m-1):  #= 2:m-1
            f = concatenate(( f, (i)*sum-1), 0)
        f = concatenate(( f, [-1]), 0)
        
    elif( f_name == 55 ) :            # Powell badly scaled  x0 = (0, 1 )
        n = 2
        f1 = 10000*x[0]*x[1]-1
        f2 = exp(-x[0])+exp(-x[1])-1.0001
        f = [ f1[0], f2[0] ]
        
    elif( f_name == 56 ) :            # Brown badly scaled  x0 = (1, 1 )
        n = 2
        f1 = x[0]-1000000
        f2 = x[1]-0.00002
        f3 = x[0]*x[1]-2
        f = [ f1[0], f2[0], f3[0] ]
        
    elif( f_name == 57 ) :            # Chebyquad x0 = ( j/n+1 )
        n  = 5
        f  = []
        T0 = 1
        for j in range(n): 
            T1 = (4*x[j]-2)
            T2 = 2*((1-2*x[j])*T1-T0)
            T3 = 2*((1-2*x[j])*T2-T1)
            T4 = 2*((1-2*x[j])*T3-T2)
            T5 = 2*((1-2*x[j])*T4-T3)
            fj = (T1+T2+T3+T4+T5)/5
            if ( mod(j+1,2) == 0 ):
                fj = fj + 1./((j+1)**2-1)
            f = concatenate(( f, fj), 0)
        
    elif( f_name == 58 ) :            # Cliff
        n  = 2
        f  = ((x[0]-3)/100)**2-x[0]+x[1]+exp(20*(x[0]-x[1]))
        
    elif( f_name == 59 ) :            # Cluster
        n = 2
        f1 = (x[0]-x[1]**2)*(x[0]-sin(x[1]))
        f2 = (cos(x[1])-x[0])*(x[1]-cos(x[0]))
        f = [ f1[0], f2[0] ]
        
    elif( f_name == 60 ) :            # Cragg -Levy
        n = 4
        f1 = (exp(x[0])-x[1])**2
        f2 = 10*(x[1]-x[2])**3
        f3 = (tan(x[2]-x[3]))**2
        f4 = x[0]**4
        f5 = (x[3]-1)
        f = [ f1[0], f2[0], f3[0], f4[0], f5[0] ]
        
    elif( f_name == 61 ) :            # Gottfr
        n = 2
        f1 = x[0]-0.1136*(x[0]+3*x[1])*(1-x[0])
        f2 = x[1]+7.5*(2*x[0]-x[1])*(1-x[1])
        f = [ f1[0], f2[0] ]
        
    elif( f_name == 62 ) :            # Hypcir
        n = 2
        f1 = x[0]*x[1]-1
        f2 = (x[0])**2+x[1]**2-4
        f = [ f1[0], f2[0] ]
        
    elif( f_name == 63 ) :            # Himln3
        n = 2
        f = x[0]**3+x[1]**2-3*x[0]-2*x[1]+2
        
    elif( f_name == 64 ) :            # Himm25
        n = 2
        f = 4*(x[0]-5)**2+(x[1]-6)**2
        
    elif( f_name == 65) :             # Himm27
        n = 2
        f = (x[0]*x[1]*(1-x[0])*(1-x[1]-x[0]*(1-x[0])**5))**2
        
    elif( f_name == 66) :             # Himm28
        n = 2
        f1 = x[0]**2+x[1]-11
        f2 = x[0]+x[1]**2-7
        f = [ f1[0], f2[0] ]
        
    elif( f_name == 67) :             # Himm29
        n = 2
        f1 = x[0]**2+12*x[1]-1
        f2 = 49*x[0]**2+49*x[1]**2+84*x[0]+2324*x[1]-681
        f = [ f1[0], f2[0] ]
        
    elif( f_name == 68) :             # Himm30
        n = 3
        f1 = 10*(x[2]-0.25*(x[0]+x[1])**2)
        f2 = 1-x[0]
        f3 = 1-x[1]
        f = [ f1[0], f2[0], f3[0] ]
        
    elif( f_name == 69) :             # Himm32
        n = 4
        a = [ 0,     0.000428, 0.001000, 0.001610, 0.002090, 0.003480, 0.005250 ]
        b = [ 7.391, 11.18,    16.44,    16.20,    22.20,    24.02,    31.32 ]
        f = []
        for i in range(7):
            u = x[0]**2+a[i]*x[1]**2+a[i]**2*x[2]**2
            v = b[i]*(1+a[i]*x[3]**2)
            f = concatenate(( f, 100*((u/v)-1)), 0)
        
    elif( f_name == 70) :             # Himm33
        n = 2
        f = exp(-x[0]-x[1])*(2*x[0]**2+3*x[1]**2)
        
    elif( f_name == 71 ) :            # Tridia
        n = 10
        f = x[0]-1
        for i in range(1,n):
            f = concatenate(( f, sqrt(i+1)*(2*x[i]-x[i-1])), 0)
        
    elif ( f_name == 101) :      # branin 
        a=1
        b=5.1/(4*pi*pi)
        c=5/pi
        d=6
        h=10
        ff=1/(8*pi)
        x1 = x[0]
        x2 = x[1]
        f=a*(x2-b*x1**2+c*x1-d)**2+h*(1-ff)*cos(x1)+h
        
        #fglob = 0.397887357729739
        #xglob = [9.42477796   -3.14159265  3.14159265 
        #         2.47499998   12.27500000  2.27500000]
        #nglob = 3  
        
    elif ( f_name == 102) :     # Six-hump camel
        x1 = x[0]
        x2 = x[1]
        f=(4-2.1*x1**2+x1**4/3)*x1**2+x1*x2+(-4+4*x2**2)*x2**2       
        #fglob = -1.0316284535
        #xglob = [ 0.08984201  -0.08984201
        #          -0.71265640   0.71265640]
        #nglob = 2
        
    elif ( f_name == 103) :     # Goldstein-Price
        x1 = x[0]
        x2 = x[1]
        
        f =(1+(x1+x2+1)**2*(19-14*x1+3*x1**2-14*x2+6*x1*x2+3*x2**2))*(30+(2*x1-3*x2)**2*(18-32*x1+12*x1**2+48*x2-36*x1*x2+27*x2**2))
        #fglob = 3
        #xglob = [0. -1.]
        #nglob = 1
        
    elif ( f_name == 104) :     # Well conditioned problem
        H  = diag(range(10))/10.
        g0 = 10.*ones((10,1), float)
        f0 = 2.
        f = f0+dot(g0.transpose(),x)+1./2.*dot(x.transpose(),dot(H,x))
        #fglob = -1.4624841e+03
        #xglob = -H\g0
        #nglob = 1
        f = f[0]
        
    elif ( f_name == 200 ) :   # Powell's problem
        n = size(x)
#-------------------------
        #rand('seed',44)
        random('seed',44)
#-------------------------
        xstar = random(n,1)
#   normdiffx=norm(xstar-x,inf)
        sigma = 1 + round( 9 * random(n,1) )
        S = reshape( round( 200 * random(2*n**2,1) - 100 ), 2*n, n )
        C = reshape( round( 200 * random(2*n**2,1) - 100 ), 2*n, n )
        f = 0
        
        for i in range(2*n):
            gi = 0
            si = 0
            for j in range(n):
                gi = gi + S[i,j] * sin(xstar[j]/sigma[j]) + C[i,j] * cos(xstar[j]/sigma[j])
                si = si + S[i,j] * sin(x[j]/sigma[j]) + C[i,j] * cos(x[j]/sigma[j])
            f = f + (gi - si)**2
            
#    elif ( f_name == 0 ) :
#        f = cuter_obj( x )
#        end
        
########################################################################
# produce noisy function values
#        
#    if ( size( varargin, 2 ) > 1 ):
#        noise = varargin{2}
#    else:
#        noise = 0
#
#    if ( noise > 0 ):        
#        f = f + noise * randn(size(f))

    if size(f) > 1:
        tmp = 0
        for i in range(size(f)):
            tmp += f[i]**2
        f = array([tmp])

    return f















#gradient of the predefined functions

def dfRosenbrock2D(x):
    dfunc = []
    der = -400*(X[1]-X[0]**2)*X[0] - 2*(1-X[0])
    dfunc.append(der)
    der = 200*(X[1]-X[0]**2)
    dfunc.append(der)
    return dfunc
def dfbeale(x):
        g1 = 2*(1.5-x[0]*(1-x[1]))*(x[1]-1)+2*(2.25-x[0]*(1-x[1]**2))*(x[1]**2-1)+2*(2.625-x[0]*(1-x[1]**3))*(x[1]**3-1)
        g2 = 2*(1.5-x[0]*(1-x[1]))*x[0] + 2*(2.25-x[0]*(1-x[1]**2))*2*x[0]*x[1] + 2*(2.625-x[0]*(1-x[1]**3))*3*x[0]*x[1]**2
        return array([g1,g2])
def dfrosenbrock(x):
        n = size(x,0)
        g = ones((n,1), float)
        for i in range(1,n-1):
            g[i] = 2*(x[i]-1) + 400*x[i]*(x[i]**2-x[i+1]) + 200*(x[i]-x[i-1]**2)
        g[0] = 2*(x[0]-1) + 400*x[0]*(x[0]**2-x[1])
        g[n-1] = 200*(x[n-1]-x[n-2]**2)
        return g
def dfbook1(x):
        g1 = 12*( (6*x[0]-2)*sin(12*x[0]-4) + (6*x[0]-2)**2*cos(12*x[0]-4) )
        g2 = [0]
        return array([g1,g2])
def dfbook2(x):
        g1 = 2*(x[1] - (5.1*x[1]/(4*pi**2)) + 5*x[0]/pi - 6)*5/pi - 10*(1-1/(8*pi))*sin(x[0]) + 5
        g2 = 2*(x[1] - (5.1*x[1]/(4*pi**2)) + 5*x[0]/pi - 6)*(1-(5.1/(4*pi**2)))
        return array([g1,g2])
def dfquadratic(x):
        n = size(x,0)
        g = ones((n,1), float)
        for i in range(0,n):
            g[i] = 2*x[i]
        return g
def dfcubic(x):
        n = size(x,0)
        g = ones((n,1), float)
        for i in range(0,n):
            g[i] = 3*x[i]**2
        return g
def dfpoly4(x):
        H = array([[4,-1],[-1,4]])
        g = dot(H,x) + array([[1],[1]])
        g[0] = g[0] + 3*x[0]**2
        g[1] = g[1] 
        return g
def dfAquadratic(x):
        n = size(x,0)
        H = -ones((n,n), float)
        for i in range(n):
            H[i,i] = 4
        b = ones((n,1), float)
        return dot(H,x) + b
def dfsin1(x):
        g1 = cos(x[0]/3)/3
        g2 = -sin(x[1]/5)/5
        return array([g1,g2])
def dfsin2(x):
        g1 = 2*x[1]*cos(x[0]+x[1]) + x[0]*cos(x[0]+x[1]) + sin(x[0]+x[1]) + 2*x[0]
        g2 = 2*x[1]*cos(x[0]+x[1]) + 2*sin(x[0]+x[1]) + x[0]*cos(x[0]+x[1]) + 6*x[1]
        return array([g1,g2]) 
def dfsin4(x):
        c = 3
        g1 = 2/c*x[1]* cos((x[0]+x[1])/c) 
        g2 = 2/c*x[1]* cos((x[0]+x[1])/c) + 2* sin((x[0]+x[1])/c)
        return array([g1,g2]) 
def dfsin3(x):
        g1 = 2*(x[1]**3)* cos((x[0])*pi)*pi + sin((x[0]+x[1])*10) + x[0]* cos((x[0]+x[1])*10)*10
        g2 = 6*(x[1]**2)* sin((x[0])*pi) + x[0]* cos((x[0]+x[1])*10) *10
        return array([g1,g2]) 
def dfplateau(x):
        return array([0.,0.]) 
def dfplateau2(x):
        if x[0] < -3:
            if x[1] < -3:
                fY = grad(x,'sin1', RS)
            elif x[1] < 0:
                fY = grad(x,'sin3', RS)
            else:
                fY = grad(x,'sin2', RS)
        elif x[0] < 0:
            if x[1] < 0:
                fY = grad(x,'sin3', RS)
            else:
                fY = grad(x,'sin2', RS)
        else:
            if x[1] < 0:
                fY = grad(x,'sin2', RS)
            else:
                fY = grad(x,'rosenbrock', RS)
        return fY
def dflinquad(x):
        n = size(x,0)
        g = ones((n,1), float)
        g[0] = 2*x[0]
        return g
def dfconstant(x):
        return array([0.,0.])
def dfsvrQP(x):
        return RS.svrQP(x,1)
        














#Hessian of the predefined functions


def ddfbeale(x):
        h11 = 2*(x[1]-1)**2+2*(x[1]**2-1)**2+2*(x[1]**3-1)**2
        h22 = 2*x[0]**2 + 2*(2*x[0]*x[1])**2 + 2*2*x[0]*(2.25-x[0]*(1-x[1]**2)) + 2*(3*x[0]*x[1]**2)**2 + 2*(2.625-x[0]*(1-x[1]**3))*6*x[0]*x[1]
        h12 = h21 = 2*(1.5-x[0]*(1-x[1]))+2*x[0]*(x[1]-1) + 2*(2.25-x[0]*(1-x[1]**2))*2*x[1]+2*2*x[0]*x[1]*(x[1]**2-1) + 2*(2.625-x[0]*(1-x[1]**3))*3*x[1]**2+2*3*x[0]*x[1]**2*(x[1]**3-1)
        return array([[h11[0],h12[0]],[h21[0],h22[0]]])
def ddfrosenbrock(x):
        n = size(x,0)
        H = zeros((n,n), float)
        H[0,0] = 2 - 400*x[1] + 1200*x[0]**2
        H[0,1] = - 400*x[0]
        H[n-1,n-1] = 200
        H[n-1,n-2] = - 400*x[n-2]
        for i in range(1,n-1):
            H[i,i-1] = -400*x[i-1]
            H[i,i] = 2 - 400*x[i+1] + 1200*x[i]**2 + 200
            H[i,i+1] = -400*x[i]
        return H
def ddfbook1(x):
        h11 = 12**2*( sin(12*x[0]-4)/2 + (6*x[0]-2)*2*cos(12*x[0]-4) - (6*x[0]-2)**2*sin(12*x[0]-4))
        h12=h21=h22 = [0]
        return array([[h11[0],h12[0]],[h21[0],h22[0]]])
def ddfbook2(x):
        h11 = 2*(5/pi)**2 - 10*(1-1/(8*pi))*cos(x[0])
        h22 = 2*(1-(5.1/(4*pi**2)))**2
        h12 = h21 = 2*(5/pi)*(1-(5.1/(4*pi**2)))
        return array([[h11[0],h12],[h21,h22]])
def ddfquadratic(x):
    n = size(x,0)
    H = zeros((n,n), float)
    for i in range(n):
        H[i,i] = 2
    return H
def ddfcubic(x):
    n = size(x,0)
    H = ones((n,n), float)
    for i in range(0,n):
        H[i,i] = 6*x[i]
    return H
def ddfpoly4(x):
        h11 = 4 + 6*x[0]
        h12 = [-1]
        h21 = [-1]
        h22 = [4]
        return array([[h11[0],h12[0]],[h21[0],h22[0]]])
def ddfAquadratic(x):
        n = size(x,0)
        H = -ones((n,n), float)
        for i in range(n):
            H[i,i] = 4
        return H
def ddfsin1(x):
        h11 = -sin(x[0]/3)/9
        h22 = -cos(x[1]/5)/25
        h12 = h21 = [0]
        return array([[h11[0],h12[0]],[h21[0],h22[0]]])
def ddfsin2(x):
        h11 =2*cos(x[0]+x[1]) - (2*x[1] + x[0]) * sin(x[0]+x[1]) + 2 
        h22 = 4*cos(x[0]+x[1]) - (2*x[1] + x[0]) * sin(x[0]+x[1]) + 6
        h12 = h21 = 3*cos(x[0]+x[1]) - (2*x[1] + x[0]) * sin(x[0]+x[1])
        return array([[h11[0],h12[0]],[h21[0],h22[0]]])
def ddfsin3(x):
        g1 = 2*(x[1]**3)* cos((x[0])*pi)*pi + sin((x[0]+x[1])*10) + x[0]* cos((x[0]+x[1])*10)*10
        g2 = 6*(x[1]**2)* sin((x[0])*pi) + x[0]* cos((x[0]+x[1])*10) *10
        h11 = - 2*(x[1]**3)* sin((x[0])*pi)*(pi**2) + 2*cos((x[0]+x[1])*10)*10 - x[0]* sin((x[0]+x[1])*10)*100
        h22 = 12*x[1]* sin((x[0])*pi) - x[0]* sin((x[0]+x[1])*10) *100
        h12 = h21 = 6*(x[1]**2)* cos((x[0])*pi)*pi + cos((x[0]+x[1])*10)*10 - x[0]* sin((x[0]+x[1])*10)*100
        return array([[h11[0],h12[0]],[h21[0],h22[0]]])
def ddflinquad(x):
    n = size(x,0)
    H = zeros((n,n), float)
    H[0,0] = 2
    return H
def ddfconstant(x):
        return array([[0,0],[0,0]])
def ddfsvrQP(x):
        return RS.svrQP(x,2)



def clinquad(x):
    n = size(x,0)
    nC = 2
    A = zeros((nC,n), float)
    for i in range(n):
        if i%2>0:
            A[0,i] = 1
    A[1,0] = 2
    return dot(A,x)

def cquadratic(x):
    n = size(x,0)
    nC = 2
    A = zeros((nC,n), float)
    for i in range(n):
        if i%2>0:
            A[0,i] = 1
    A[1,0] = 2
    return dot(A,x)

def crosenbrock(x):
    n = size(x,0)
    nC = 2
    A = zeros((nC,n), float)
    for i in range(n):
        if i%2>0:
            A[0,i] = 1
    A[1,0] = 2
    return dot(A,x)

def cbeale(x):
    n = size(x,0)
    nC = 2
    A = zeros((nC,n), float)
    for i in range(n):
        if i%2>0:
            A[0,i] = 1
    A[1,0] = 2
    return dot(A,x)
