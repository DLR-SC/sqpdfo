""" 
Minimizes an affine quadratic form subject to simple bounds, 
using coordinate searches and reduced subspace minimizations 
using LDL^T factorization updates 
    min    fval = gamma + c^T x + 0.5 x^T G x } 
    s.t.   x in [xu,xo] 
   
where G is a symmetric (n x n) matrix, not necessarily definite 
(if G is indefinite, only a local minimum is found) 
 
If G is sparse, it is assumed that the ordering is such that 
a sparse modified Cholesky factorization is feasible 
   
Usage: 
------ 
[x,fval] = minq( gamma, c, G, xu, xo, prt, xx) 
 
History: 
-------- 
1) Ziwen Fu. 10/22/2008 
   Directly Translate from minq.m ( v2.1 ) 
  
2) Ziwen Fu. 10/25/2008 
   Python Class version of Minq 
""" 
   
__docformat__ = "restructuredtext en" 
   
   
import numpy 
from minq_util import toRow, toCol, find, mldiv, crdot
eps   = numpy.finfo('d').eps 
inf   = numpy.inf 
   
   
# ============================================================= 
class Minq: 
    """ 
    Minimizes an affine quadratic form subject to simple bounds, 
    using coordinate searches and reduced subspace minimizations 
    using LDL^T factorization updates 
        min    fval = gamma + c^T x + 0.5 x^T G x } 
        s.t.   x in [xu,xo] 
   
    where G is a symmetric (n x n) matrix, not necessarily definite 
    (if G is indefinite, only a local minimum is found) 
 
    If G is sparse, it is assumed that the ordering is such that 
    a sparse modified Cholesky factorization is feasible 
    """ 
    convex = 0 
   
    # The maximal number of iterative refinement steps 
    nitrefmax   = 3 
   
    # The perturbation in last two digits 
    neps = 100*eps 
   
    # The number of subspace steps 
    nsub = 0 
   
    # The allow variables to be freed in csearch? 
    unfix = 1 
   
    # No iterative refinement steps so far 
    nitref = 0 
   
    # Improvement expected 
    improvement = 1 
 
    # Best function value 
    fval = inf 
   
    nfree     = 0 
    nfree_old = -1 
  
    def __init__(self, gamma, c, G, xu, xo, prt=0, x0=None): 
        """ 
        Inputs: 
        ------- 
        gamma  a constant. 
        c      a colomn vector. 
        G      a symmetric (n x n) matrix. 
        xu     lower bound 
        xo     upper bound 
 
        Optional Inputs: 
        ---------------- 
        prt    print level. 
        xx     initial guess (optional) 
        """ 
        self.gamma = gamma 
        self.c   = c 
        self.G   = numpy.mat(G) 
        self.xu  = xu 
        self.xo  = xo 
        self.prt = prt 
        self.xTry  = x0 
        self.n  = self.G.shape[0] 
   
        # Initialization 
        if  prt>0: 
            self.print_init() 
   
   
    def print_minq(self): 
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' )
        print('!!!!!           Minq          !!!!!' )
        print('!!!!! incomplete minimization !!!!!' )
        print('!!!!!   too many iterations   !!!!!' )
        print('!!!!!   increase maxit        !!!!!' )
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' )
   
   
    def print_init(self): 
        print('====================================================' )
        print('======================  Minq  ======================' )
        print('====================================================' )
   
   
    def setSearchInit(self): 
   
        # Initialize trial point xx, function value fval and gradient g 
        if  self.xTry is None: 
            # Cold start with absolutely smallest feasible point 
            self.xTry = numpy.empty( self.n, 'd' )
  
        # Force starting point into the box 
        self.xTry = numpy.maximum( self.xu, numpy.minimum(self.xTry,self.xo) )

        # Maximal number of iterations 
        # This limits the work to about 1+4*maxit/n matrix multiplies 
        # Usually at most 2*n iterations are needed for convergence 
        self.maxit = 3 * self.n 
   
        # Regularization for low rank problems, perturbation in last two digits 
        for i in range( self.n ):
            self.G[i,i] *= (1+self.neps) 
   
        # Initialize LDL^T factorization of G_KK 
        self.K = numpy.zeros( self.n )  # initially no rows in factorization 
        self.L = numpy.eye(self.n) 
        self.d = numpy.ones( self.n ) 
   
        # Dummy initialization of indicator of free variables 
        # will become correct after first coordinate search 
        self.free = [0] * self.n 
   
   
    def _quadMin(self, low, high, b, a): 
        """ 
        Get minimizer x in [low, high] for a univariate quadratic 
   
        q(x) = 0.5*a*x^2 + b*x 
  
        Output: 
        ------ 
        x:    minimizer x in [low, high] 
        lba:  lower bound active 
        uba:  upper bound active 
  
        ier:  0 (finite minimizer) 
              1 (unbounded minimum) 
        """ 
        lba = 0  # lower bound is active or not? 
        uba = 0  # upper bound is active or not? 
   
        # Determine unboundedness 
        ier=0 
        if  (low  == -inf) and ( a<0 or (a==0 and b>0) ): 
            ier=1; lba=1 
   
        if  (high ==  inf) and ( a<0 or (a==0 and b<0) ): 
            ier=1; uba=1 
   
        if  ier: 
            return (numpy.nan, lba, uba, ier )
   
        # Determine activity 
        if  a==0 and b==0: 
            x=0 
  
        elif a<=0: 
            # Concave case minimal at a bound 
            if   low  == -inf:  lba = 0 
            elif high ==  inf:  lba = 1 
            else:               lba = (2*b+(low+high)*a>0) 
            uba = not lba 
   
        else: 
            x   = -float(b)/float(a)  # Unconstrained optimal step 
            lba = (x <= low )         # Lower bound is active? 
            uba = (x >= high)         # Upper bound is active? 
  
        if lba: x = low 
        if uba: x = high 
   
        return (x, lba, uba, ier ) 
   
   
    def _Info(self): 
   
        self.info = 0 
        if not self.improvement: # Good termination 
            if  self.prt: 
                print('Terminate: no improvement in coordinate search' )
   
        elif self.nitref > self.nitrefmax:  # Good termination 
            if  self.prt:  print('Terminate print_minq: nitref>nitrefmax' )
   
        elif self.nitref > 0 and self.nfree_old == self.nfree and self.fvalnew >= self.fval: # Good termination
            if  self.prt:
                print('Terminate:nitref>0&nfree_old==nfree&fvalnew>=fval' )
   
        else: 
            self.info = -1 
   
        if self.nitref==0 and self.nsub == self.maxit: 
            if self.prt: self.print_minq() 
            else:   print('iteration limit exceeded' )
            self.info = 99 
   
   
    def ldlrk1(self, L, d, alpha, u):
        """ 
        Computes LDL^T factorization for LDL^T + alpha * u * u^T 
        if alpha>=0 or if new factorization is definite(both signalled by p=[]) 
        otherwise, the original L,d and 
        a direction p of null or negative curvature are returned 
   
        d contains diag(D) and is assumed positive 
   
        Warning: does not work for dimension 0 
        """ 
        L = numpy.asarray(L) 
        p = numpy.array([]) 
        if alpha==0: 
            return 
        n    = u.size 
   
        # Save old factorization 
        L0 = L.copy() 
        d0 = d.copy() 
   
        # Update 
        for k in find(u!=0): 
   
            delta = d[k] + alpha * u[k]**2 
            if alpha<0 and delta <= self.neps: 
                # Update not definite 
                p      = numpy.zeros(n, 'd') 
                p[k]   = 1 
                p[0:k+1] = toRow ( mldiv( numpy.asmatrix(L[0:k+1,0:k+1].T), 
                                        toCol( numpy.asarray( p[0:k+1] ) ) 
                                        ) ) 
   
                # Restore original factorization 
                L = L0 
                d = d0 
                return (L, d, p) 
   
            q    = d[k]/delta 
            d[k] = delta 
   
            # ---------------------------------------------- 
            ind = range(k,n) 
            c   = L[ind,k]*u[k] 
            L[ind,k] = L[ind,k]*q + (alpha*u[k]/delta)*u[ind] 
            u[ind] -= c 
   
            alpha *= q 
            if alpha==0: 
                break 
   
        return (L, d, p) 
   
   
    def ldlup(self, L, d, j, g): 
        """ 
        Updates LDL^T factorization when a unit j-th row and column 
        are replaced by column g 
   
        If the new matrix is definite (signalled by p=[]); 
        otherwise, the original L,d and 
        a direction p of null or negative curvature are returned 
   
        d contains diag(D) and is assumed positive 
        Note that g must have zeros in other unit rows!!! 
        """ 
        p = numpy.array([]) 
        n = d.size 
        if j==0: 
            delta = g[j] 
            if delta <= self.neps: 
                p    = numpy.zeros( n, 'd') 
                p[0] = 1.0 
                return (L,d,p) 
   
            d[j] = delta 
            return (L, d, p) 
   
        # Now j>1, K nonempty 
        LII  = L[0:j,0:j]
        tmat = g[0:j]
        u    = mldiv(LII, tmat )
        v    = u / numpy.array([d[0:j]]).T
        delta = g[j] - numpy.dot(u.T, v)
        if  delta <= self.neps: 
            pLII = mldiv(LII.T,v) 
            p    = numpy.zeros( n, 'd') 
            for i in range( len(pLII) ):
                p[i] = pLII[i] 
            p[j] = -1 
            return (L, d, p) 
   
        LKI   = L[j+1:n,0:j] 
        LKI_u = LKI*numpy.mat( toCol(u) ) 
   
        w = ( toCol(g[j+1:n]) - LKI_u )/delta 
   
        (LKK,d[j+1:n],q) = self.ldlrk1( L[j+1:n,j+1:n], 
                                        d[j+1:n], 
                                        -delta, 
                                        toRow(w) 
                                        ) 
   
        if len(q)==0: 
            # Work around expensive sparse L(K,K)=LKK 
            L1  = L[0:j,:]
            L2  = numpy.mat( numpy.zeros(n) )
            L3  = numpy.mat(L[j,j+1:n])
            L_1 = numpy.mat( 1 ) 
            Lw  = numpy.mat( w ) 
            Lv  = numpy.mat( v.T )
            L4  = numpy.mat( numpy.zeros( (n-j-1,1) ) )
            #print('L1='+str(L1))
            #print('Lv='+str(Lv))
            #print('L_1='+str(L_1))
            #print('L3='+str(L3))
            #print('LKI='+str(LKI))
            #print('Lw='+str(Lw))
            #print('LKK='+str(LKK))
            L   = numpy.bmat( 'L1; Lv,L_1,L3; LKI,Lw,LKK' ) 
            d[j]=delta 
   
        else: # Work around expensive sparse L(K,K)=LKK 
            L1 = L[0:j+1,:] 
            L2 = L[j+1:n,j] 
            L  = numpy.bmat( 'L1; LKI,L2,LKK' ) 
            pi = crdot( w, q ) 
            LKI_q = LKI.T * numpy.mat( toCol(q) ) 
            pi_v  = pi*toCol(v) 
   
            p1 = mldiv( LII.T, pi_v-LKI_q) 
   
            p = numpy.vstack( ( toCol(p1), numpy.array(-pi), toCol(q) ) ) 
   
        return  (L, d, p) 
   
   
   
    # --------------------------------------------------------------------- 
    def ldldown(self, L, d, j): 
        """ 
        Downdates LDL^T factorization 
        when j-th row and column are replaced by j-th unit vector 
   
        d contains diag(D) and is assumed positive 
        """ 
        L = numpy.asmatrix(L) 
        n = d.size 
   
        w =  toRow( L[j+1:n,j] ) 
   
        if  j<n: 
            [ LKK, d[j+1:n], q ] = self.ldlrk1( L[j+1:n,j+1:n], 
                                                d[j+1:n], 
                                                d[j], 
                                                w 
                                                ) 
            # Work around expensive sparse L[K,K]=LKK 
            L1 = L[0:j,] 
            L2 = numpy.mat( numpy.zeros(n) ) 
            L3 = L[j+1:n, 0:j] 
            L4 = numpy.mat( numpy.zeros( (n-j-1,1) ) ) 
            L  = numpy.bmat( 'L1; L2; L3,L4,LKK' ) 
            L[j,j] = 1 
   
        else: 
            for i in range(n-1):
                L[n,i] = 0.0 
   
        d[j]=1 
   
        return (L,d) 
   
   
    def coordinateSearch(self): 
        """ 
        Coordinate search 
        """ 
        # The number of consecutive free steps 
        count = 0 
  
        # Current coordinate searched 
        k = -1 
        while 1: 
   
            while count <= self.n: 
                # Find next free index (or next index if unfix) 
                count += 1 
                if k == self.n-1:  k=0 
                k += 1  #? 
   
                if  self.free[k] or self.unfix: 
                    break 
   
            if  count > self.n: 
                # Complete sweep performed without fixing a new active bound 
                break 
   
            q       = self.G[:,k] 
            alpha_u = self.xu[k] - self.x[k] 
            alpha_o = self.xo[k] - self.x[k]  # bounds on step 
   
            # Find step size 
            [alpha, lba, uba, self.info] = self._quadMin( alpha_u, alpha_o, self.gradient[k], q[k] )
   
            if  self.info: 
                x = numpy.zeros( self.n  ) 
                if lba: x[k] = -1 
                else:   x[k] = 1 
                return 
   
            # New x 
            xnew = self.x[k] + alpha 
            if  self.prt and self.nitref>0: 
                print(xnew,alpha)
   
            if  lba or xnew <= self.xu[k]:  # Lower bound active 
                if  alpha_u !=0 : 
                    self.x[k] = self.xu[k] 
                    self.gradient += alpha_u * q 
                    count=0 
                self.free[k] = 0 
   
            elif uba or xnew >= self.xo[k]: # Upper bound active 
                if  alpha_o != 0: 
                    self.x[k] = self.xo[k] 
                    self.gradient += alpha_o * q 
                    count = 0 
                self.free[k] = 0 
   
            else: # No bound active 
                if  alpha != 0: 
                    if self.prt>1 and not self.free[k]: 
                        unfixstep=[ self.x[k], alpha ]
  
                    self.x[k] = xnew 
                    self.gradient += alpha * q 
                    self.free[k] = 1 
        # ============== End of coordinate search ================ 
   
   
    def subspace(self): 
        """ 
        Take a subspace step 
        """ 
        self.nsub += 1 
   
        # Downdate factorization 
        for j in find( self.free < self.K ):    # list of newly active indices 
            [self.L, self.d] = self.ldldown(self.L, self.d, j) 
            self.K[j] = 0 
   
        # Update factorization or find indefinite search direction
        definite = 1 
        for j in find( self.free > self.K ):  # List of newly freed indices 
            p      = numpy.zeros(( self.n,1 ))
            ind    = find( self.K > 0 )
            p[ind] = self.G[ind,j]
            p[j]   = self.G[j,j]
            [self.L, self.d, p] = self.ldlup(self.L, self.d, j, p)
   
            if p.any():
                definite = 0
                if self.prt: print('indefinite or illconditioned step' )
                break 
   
            self.K[j]=1

        if definite:  # Find reduced Newton direction
            p = numpy.zeros((self.n,1))
            ind = find( self.K > 0 )
            p[ind] = self.gradient[ind]
            p = -mldiv( self.L.T, mldiv(self.L,p)/numpy.array([self.d]).T )

        # Set tiny entries to zero
        p = list(p.T[0])
        p = (self.x + p) - self.x
        ind = find( p != 0 ) 
        if  len(ind)==0: # Zero direction 
            if self.prt:  print('zero direction')
            self.unfix=1 
   
        # Find range of step sizes
        pp = p[ind]
        oo = (self.xo[ind] - self.x[ind]) / pp
        uu = (self.xu[ind] - self.x[ind]) / pp 
   
        if len(oo[pp<0])==0: max_oo = -inf 
        else:                max_oo = max(oo[pp<0]) 
        if len(uu[pp>0])==0: max_uu = -inf 
        else:                max_uu = max(uu[pp>0]) 
        alpha_u = numpy.maximum( max_oo, max_uu ) 
   
        if len(oo[pp>0])==0: min_oo = inf 
        else:                min_oo = min(oo[pp>0]) 
        if len(uu[pp<0])==0: min_uu = inf 
        else:                min_uu = min(uu[pp<0]) 
        alpha_o = numpy.minimum( min_oo, min_uu ) 
   
        if  alpha_o<=0 or alpha_u>=0: 
            print('programming error: no alpha')
   
        # Find step size 
        gTp  = float( numpy.dot(self.gradient.T, p) ) 
        agTp = float( numpy.dot(abs(self.gradient.T), abs(p)) ) 
   
        # Linear term consists of roundoff only 
        if abs(gTp) < self.neps * agTp: 
            gTp=0 
   
        pTGp = numpy.dot(p, self.G*toCol(p) ) 
   
        if self.convex:                 pTGp = numpy.max(0,pTGp) 
        if (not definite) and (pTGp>0): pTGp = 0 
   
        [alpha, lba, uba, info] = self._quadMin(alpha_u, alpha_o, gTp, pTGp) 
        if  info: 
        #    self.x = numpy.zeros(n,1)
            if  lba: self.x=-p 
            else:    self.x=p 
  
        # Allow variables to be freed in csearch? 
        self.unfix = not (lba or uba) 
   
        # Update of xx 
        for k in range(len(ind)):
            # Avoid roundoff for active bounds 
            ik = ind[k] 
            if  alpha == uu[k]: 
                self.xTry[ik] = self.xu[ik] 
                self.free[ik] = 0 
   
            elif alpha == oo[k]: 
                self.xTry[ik] = self.xo[ik] 
                self.free[ik] = 0 
   
            else: 
                self.xTry[ik] += alpha * p[ik] 
   
            if  abs( self.xTry[ik] ) == inf: 
                print('infinite xx in Minq' )
   
        self.nfree   = numpy.sum( self.free ) 
        self.subdone = 1 
   
   
   
    def subspaceSearch(self): 
        """ 
        Subspace search 
        """ 
        if  not self.improvement or self.nitref > self.nitrefmax: 
            # Optimal point found - nothing done 
            pass 
   
        elif self.nitref > self.nitrefmax: 
            # Enough refinement steps - nothing done 
            pass 
   
        elif self.nfree == 0: 
            # No free variables - no subspace step taken 
            if  self.prt>0: 
                print('no free variables-no subspace step taken')
            self.unfix=1 
 
        else: 
            self.subspace() 
  
   
    def calcGain(self): 
        """ Calculate gain """ 
        return self.fval - self.gamma - 0.5* numpy.dot( self.x, self.c + toRow(self.gradient) )
   
   
    def calc_g(self): 
        """ Calculate gradient """ 
        return self.G * toCol(self.xTry) + toCol(self.c) 
   
 
    def calc_newfval(self): 
        return self.gamma + 0.5 * numpy.dot( self.c + toRow(self.gradient), self.xTry)
   
   
    def pr01(self, name,x): 
        """ 
        Prints a (0,1) profile of x and returns the number of nonzeros 
   
        x: a numpy array 
        """ 
        # Check the size of the vector 
        n     = len(x)
        print('n='+str(n))
        text  = name + ': ' 
        summe = 0 
        for k in range(n):
            if x[k]: 
                text += '1' 
                summe += 1 
            else: 
                text += '0' 
 
        print(text + '   ', summe, ' nonzeros')
        return  summe 
   
   
    def search(self): 
        """ 
        Main loop: alternating coordinate and subspace searches 
   
        Outputs: 
        -------- 
        x      minimizer (but unbounded direction if info=1) 
        fval   optimal function value 
        nsub   the number of subspace steps 
        info   0  (local minimizer found) 
               1  (unbounded below) 
               99 (maxiter exceeded) 
        """ 
        self.setSearchInit() 
        while 1: 
   
            self.gradient = self.calc_g() 
            self.fvalnew  = self.calc_newfval() 
 
            if  self.nitref==0: 
                self.fval = numpy.minimum(self.fval, self.fvalnew) 
            else:  # More accurate g and hence f if nitref>0 
                self.fval = self.fvalnew 
 
            self.x = self.xTry 
            self._Info() 
            if self.info == 0 or self.info == 99: break 
 
            # ------------------------------------------------------------ 
            # Coordinate search 
            # ------------------------------------------------------------ 
            self.coordinateSearch() 
            if  self.info: 
                return  (self.x, self.fval, self.info, self.nsub) 
 
            # ============================================================ 
            self.nfree = numpy.sum( self.free ) 
            if  self.unfix and self.nfree_old == self.nfree: 
                # In exact arithmetic, we are already optimal 
                # recompute gradient for iterative refinement 
                self.gradient = self.calc_g() 
                self.nitref += 1 
                if self.prt>0: print('Optimum found.iterative refinement tried' )
            else: 
                self.nitref = 0 
 
            self.nfree_old = self.nfree 
 
            # Calculate gain 
            self.gain_cs = self.calcGain() 
 
            # Is improvement? 
            self.improvement = (self.gain_cs>0) or (not self.unfix) 
 
            if  self.prt: 
                # print (0,1) profile of free and return the number of nonzeros 
                self.nfree = self.pr01('csrch ',self.free) 
                print(self.gain_cs)
 
            # ------------------------------------------------------------- 
            # Subspace search 
            # ------------------------------------------------------------- 
            self.subspaceSearch() 
            if  self.info: 
                return (self.x, self.fval, self.info, self.nsub) 
 
            # Show it 
            if  self.prt: 
                # Print (0,1) profile of free and return the number of nonzeros 
                self.nfree = self.pr01('ssrch ', self.free)
                if self.unfix and self.nfree < self.n:
                    print('Bounds may be freed in next csearch')
 
        return self.x, self.fval, self.info, self.nsub
 

# ========================================================================== 
def minq( c, G, xlow, xhigh, x0=None, gamma=0, prt=0 ): 
    """ 
    Minimizes an affine quadratic form subject to simple bounds, 
    using coordinate searches and reduced subspace minimizations 
    using LDL^T factorization updates 
         min    fval = gamma + c^T x + 0.5 x^T G x } 
         s.t.   x in [xu,xo] 
 
    Inputs: 
    ------- 
    c      a colomn vector. 
    G      a symmetric (n x n) matrix. 
    xlow   lower bound 
    xhigh  upper bound 
 
    Optional Inputs: 
    ---------------- 
    gamma  a constant. 
    prt    print level. 
    x0     initial guess. 
 
    Output: 
    ------- 
    x      minimizer (but unbounded direction if info=1) 
    info   0  (local minimizer found) 
           1  (unbounded below) 
           99 (maxit exceeded) 
    """ 
    m = Minq(gamma, c, G, xlow, xhigh, prt=prt, x0=x0) 
    [x, fval, info, nsub] = m.search()  # fval   optimal function value 
  
    return (x, info) 
  
  
  
  
# ---------------------------------------------------------------------------- 
if __name__ == '__main__': 
    cc = numpy.array([21.0370,23.9598,23.0826,10.5767,18.5893,19.8627,16.4453]) 
  
    G = numpy.matrix([[ 8.6627, 7.5820, 9.2232, 2.5919, 7.0716, 8.8219, 5.7268], 
                      [ 7.5820, 12.0645,8.8799, 5.3508, 6.9379, 6.0264, 7.3194], 
                      [ 9.2232, 8.8799, 11.5346,4.0942, 7.5342,10.6629, 6.4869], 
                      [ 2.5919, 5.3508, 4.0942, 4.8251, 3.4571, 2.7991, 3.1634], 
                      [7.0716,  6.9379, 7.5342, 3.4571, 6.7366, 7.0169, 4.4661], 
                      [8.8219, 6.0264, 10.6629, 2.7991, 7.0169, 11.2640,5.0162], 
                     [5.7268, 7.3194, 6.4869, 3.1634, 4.4661,  5.0162,6.2784]]) 
  
    yu = numpy.array([0,0,0,0,0,0,0]) - 100000000.0 
    yo = numpy.array([1000]*7)  + 10000000.0
    prt = 0
  
    m = Minq(0, cc, G, yu, yo, prt)
  
    [y,fval,info, nsub] = m.search() 
  
    print("y=", y)
    print("fval=", fval)
    print("info=", info)
    print("nsub=",nsub)
