# -*- coding: utf-8 -*-

from runtime import *
from sqpdfo_global_variables import *
from numpy import array


class tcgInfo():
    def __init__(self):
         self.flag = None
         self.iter = None
         self.prec = None
         self.curv = None
         self.plevel = 0

def sqpdfo_truncated_cg_(A=None,b=None,delta=None,max_iter=None,tol=None,*args,**kwargs):
    """
    # [x,info] = sqpdfo_truncated_cg (A,b,delta,max_iter,tol);
    #
    # This function solves the system A x = b for x, by Steihaug's conjugate gradient method.
    # A is a symmetric (possibly indefinite) matrix.
    # Calculation starts from x = 0 and the returned x will be the computed step.
    #
    # Inputs:
    # A : symmetric matrix
    # b : right-hand-side
    # delta : trust-region radius
    # max_iter : maximum number of CG-iterations
    # tol : accuracy tolerance
    #
    # Outputs:
    # x : computed solution of the linear system
    # info :
    #     info.flag is the return code
    #       = 0: convergence is obtained (up to the given tolerance tol)
    #       =-1: stop on max_iter,
    #       = 1: the boundary of the ball of radius delta is encountered:
    #            the final x is on the boundary of the ball
    #       = 2: a negative curvature direction has been encountered: the
    #            final x is on the boundary of the ball.
    #     info.iter is the number of CG iterations
    #     info.prec is the final precision, the l2-norm of the residual A*x-b
    """

    info = tcgInfo()

    # Initialization

    x = zeros_(size_(b));  # initial x = 0
    cost = 0;  # initial cost
    g = -b;  # initial gradient (initial iterate is 0)
    g2 = g.T.dot(g);  # initial |g|^2

    tol2 = tol * tol;  # square of the tolerance
    delta2 = delta * delta;  # square of the trust region radius

    dAd = array([]);  # for the test isempty below
    if info.plevel:
        print('    TCG solver; required tolerance %8.2e\n' % (tol))
        print('    iter       cost        |res|   curvature  stepsize   |step|\n')

# -----------------------------------------------------------------------
# Main loop.
# -----------------------------------------------------------------------
# At the beginning of the loop the following variables are supposed
# known:
# - g (gradient at the current iterate)
# - g2 = g'*g
# - g2_ = g_'*g_ (g_ is the previous gradient, if iter > 1)
# - d (previous direction, if iter > 1)
# -----------------------------------------------------------------------

    set_iter(0)

    while 1:

        set_iter(get_iter() + 1)
        if info.plevel:
            print('    %4i  %14.7e  %7.1e'%(get_iter(),cost,sqrt_(g2)))

        # stop on convergence

        if g2 <= tol2:
            info.flag=0
            break

        # stop on iteration

        if get_iter() > max_iter:
            set_iter(copy(max_iter))
            info.flag=- 1
            break

        # new conjugate direction

        if get_iter() == 1:
            d=- g
        else:
            d=- g + (g2 / g2_) * d

        # matrix*vector product and curvature verification

        Ad=A.dot(d)
        dAd=d.T.dot(Ad)
        if info.plevel:
            print('  %9.2e'%(dAd / (d.T.dot(d))))
        if dAd <= 0: # negative curvature direction
            x,alpha=comp_dir_in_delta(x,x + d,delta,nargout=2)
            info.flag=2
            if info.plevel:
                print('  %8.2e  %8.2e\n'%(alpha,norm_(x)))
                cost=0.5 * (x.T.dot(A.dot(x))) - b.T.dot(x)
                print('    %4i  %14.7e\n'%(get_iter() + 1,cost))
            break
        #else:
        #   fprintf_(fout,'\n');

        # new iterate

        alpha=- (g.T.dot(d)) / dAd
        xx=x + alpha * d
        if info.plevel:
            print('  %8.2e'%(alpha))

        # intersection with the sphere

        if xx.T.dot(xx) > delta2: 	# new iterate is outside the trust region
            x,alpha=comp_dir_in_delta(x,xx,delta,nargout=2)
            info.flag=1
            if info.plevel:
                print('  %8.2e  %8.2e\n'%(alpha,norm_(x)))
                cost=0.5 * (x.T.dot(A.dot(x))) - b.T.dot(x)
                print('    %4i  %14.7ee\n'%(get_iter()+ 1,cost))
            break
        else:
            x=copy(xx)

        if info.plevel:
            print('  %8.2e\n'%(norm_(x)))

        # new gradient and cost

        g=g + alpha * Ad
        g2_=copy(g2)
        g2=g.T.dot(g)
        if info.plevel:
            cost=0.5 * (x.T.dot((g - b)))

    # save final information

    info.iter=get_iter()
    info.prec=sqrt_(g2)
    if not isempty_(dAd):
        info.curv=dAd / (d.T.dot(d))
    return x,info

###############################################################################

def comp_dir_in_delta(dc=None,dn=None,delta=None,*args,**kwargs):
    """
    # [dir,t] = comp_dir_in_delta (dc, dn, delta)
    #
    # Find the step dir at the intersection of the sphere of radius delta
    # (>0) and the half-line dc -> dn. It is assumed that norm(dc) < min
    # (delta,norm(dn)).
    """

    dir=dn - dc
    dir_inner=dir.T.dot(dir)
    if dir_inner == 0:
        dir=copy(dc)
        t=0
        return dir,t
    bb=dc.T.dot(dir)
    dc_inner=dc.T.dot(dc) - delta ** 2
    if dc_inner >= 0:
        dir=copy(dc)
        t=0
        return dir,t
    t=(sqrt_(bb ** 2 - dir_inner * dc_inner) - bb) / dir_inner
    dir=dc + t*(dir)
    return dir,t
