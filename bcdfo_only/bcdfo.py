#!/usr/local/bin/python
import os
import sys
from numpy import *
#from numpy.linalg import lstsq
from math import *
import bcdfo_c
import bcdfo_build_QR_of_Y
from bcdfo_poisedness_Y import *
from bcdfo_augmX_evalf import *
from bcdfo_find_smallf import *
from bcdfo_computeP import *
from bcdfo_gradP import *
from bcdfo_projgrad import *
from bcdfo_main import *

def bcdfo( f, df, x0, *args ):
   """
 
   A model-based trust-region algorithm for derivative-free bound-constrained
   minimization. Its purpose is to find a local mimimizer of the problem
 
                      min      t( f( x ) )     s.t. xl < x < xu
                    x in R^n
 
   where f is a smooth (preferably twice continuously differentiable)
   function from R^n into R^m and t(.) is the identity for m = 1 and the
   Euclidean norm if m > 1. The derivatives of f are not needed by the algorithm
   whereas there exist several options to incorporate first derivative
   information when these are available and supplied by the user.
   A starting point x0 and the bounds xl and xu must be provided by the user.
 
   The method builds a local (at most fully quadratic) interpolation model
   from known function values, and then minimizes this model within a trust
   region whose radius is updated as necessary. Geometry checks, a common
   feature of such methods, is reduced to a minimum for improved performance.
 
   The algorithm is the trust-region method suggested by Scheinberg and
   Toint in the reference below, and is based on using the Lagrange fundamental
   polynomials for interpolating known function values.  These polynomials are
   not explicitely recurred, but a QR factorization of the interpolation
   system's matrix is maintained instead (the Lagrange polynomials are then
   given by the rows of the inverse of that matrix), possibly after shifting
   and scaling to improve numerical conditioning.
 
   The objective function value is computed by a user-supplied function
   whose call is given by
       fx = f( x )
   where x is the vector at which evaluation is required.   For example,
   minimizing the hyperbolic cosine over x (starting from 1)  can be achieved
   by the call
       [ x, fx, gx, msg, nit, neval ] = bcdfo( @cosh, 1 )
   Objective functions  dep ing on parameters as in
       fx = f( x, parameter1, parameter2 )
   may also be used. For instance, if the function myobj is specified by
       function fx = myobj( x, p )
       fx = x' * p + 0.5 * x' * x;
   then minizing this function over x can be achieved by the statements
       x0 = zeros( 5, 1 );
       p  = ones( size( x0 ) );
       [ x, fx, gx, msg, nit, neval ] = bcdfo( @(x) myobj( x, p ), x0 );
   See the description of the first input argument below for more details
   on how to specify the "function_handle" for the objective function.
 
   Termination occurs as soon as the Euclidean norm of the model's gradient
   is less or equal to epsilon AND the estimated error on this gradient
   (computed by estimating the poisedness of the interpolation set) is
   below factor_CV * epsilon.
 
   The algorithmic parameters used may be specifed by the user (overriding
   their default values) via optional arguments.  Each parameter is specified
   by two consecutive such arguments, the first being a string containing the
   keyword associated with the parameter and the second its value (a number
   or a string).  For example, minimizing the objective function myobj from x0
   with a termination accuracy of 0.00001, a maximum number of evaluations of
   200 and an initial linear model is specified by the call
 
    [ x, fx, gx, msg, nit, neval ] = bcdfo( @myobj, x0, 'epsilon', 0.00001, ...
                                    'maxeval', 200, 'degree_m0', 'linear' )
 
   The implemented algorithm allows for scaling of the variables (if parameter
   scaleX = 1 and a vector of scaling factors scalefacX is provided by the user)
   where         x(i) * scalefacX(i) = approx. 1 for all i
   at the solution.
   Please do not confuse this user-defined variables scaling with the internal
   shifting and scaling of the interpolation set which is applied if parameter
   shift_Y is set to 1 and which shifts the current iterate to the origin and
   scales the interpolation set to be contained in a ball of radius one.
 
   The algorithm contains the option (if parameter shrink_Delta is set to 1) to
   shrink the trust-region radius at each unsuccessful iteration. This has been
   shown to work well for non-noisy problems and to accelerate convergence of
   the algorithm in average. The option has not been shown to perform well for
   noisy problems and shouldn't been applied therefore if the function is known
   to be noisy.
 
   See below for the full list of algorithmic parameters and associated
   keywords.
 
   The package features a "restart" option.  If this option is activated, the
   algorithmic parameters and the current interpolation set are saved on a
   file ('bcdfo.restart' in the current directory), allowing the optimization
   to be restarted from the saved state.  If such a restart is  requested, the
   algorithm is restarted with default algorithmic parameter values equal to
   the saved values.  These can then be modified by the user in the calling
   sequence (just as for ordinary defaults). The package may therefore be used
   in four different modes:
   1) a mode ('none') where restart information is neither used or saved,
   2) a mode ('save') where restart information is saved before exit,
   3) a mode ('use') where saved information is used to restart the algorithm
      ( possibly overriding the specification of x0),
   4) a mode ('use_and_save') which both uses saved information to start and
      saves information before exit.
 
 
   INPUT (mandatory):
 
 
   f           : a function handle pointing to the function to be minimized.
                 If one wishes to minimize the function whose calling sequence
                 is "fx = func(x)", then this argument can be either "@func"
                 "@(x) func(x)".  If the objective function dep s on
                 parameters, these may be specified. For instance, if one
                 wishes to minimize (over x) the function whose calling
                 sequence is "fx=func(x,parameter1,parameter2)", then the
                 argument becomes  "@(x) func(x,parameter1,parameter2)", and
                 the parameters are then automatically passed to func when
                 it is called within bcdfo  (provided they have been properly
                 assigned values in the program calling bcdfo).
                 NOTE: On restart, it is the responsability of the user to
                       provide an objective function handle which is coherent
                       with that used for the call where the restart
                       information was saved.
   df          : a function handle pointing to the gradient
   x0          : the starting point for the minimization
 
 
   INPUT (optional):
 
 
   Delta0       : the initial trust-region radius (if <=0, Delta0 = 1 is used)
                  Default: 1
   epsilon      : the threshold such that convergence is declared as soon as the
                  norm of the model's gradient and its estimated error both fall
                  below it (>0)
                  Default: 0.00001
   maxeval      : the maximum number of objective function evaluations (>0)
                  Default: 200*length(x0)
   maxit        : the maximum number of iterations (>0)
                  Default: = maxeval
   degree_m0    : the size of the initial interpolation polynomial:
                  'linear' starts with a fully linear model
                  'diagonal' starts with a linear + diagonal quadratic model
                  'quadratic' starts with a fully quadratic model
                  Default: 'linear'
   degree_mR    : the minimum size of the interpolation polynomial after its
                  recomputation beyond the first iteration:
                  'linear' reconstructs a fully linear model
                  'diagonal' recostructs a linear + diagonal quadratic model
                  'quadratic' reconstructs a fully quadratic model
                  Default: 'diagonal'
   geometry_Y0  : strategy to choose initial interpolation set:
                  'random' selects interpolation by a random procedure
                          whose geometry is optimized
                  'simplex' selects the points as the vertices of the unit
                          simplex and the mid-points of its edges.  This
                          option is only available for 'quadratic' initial
                          geometry.
                  Default: 'simplex'
   verbosity    : the volume of output:
                  'silent' = no output,
                  'low'    = a one-line summary per iteration,
                  'medium' = a one-line summary per iteration + a sumary of x
                  'high  ' = more detail
                  Default: 'low'
   show_errg    : 1 = compute gradient error estimate at every iteration
                  and print it out in the iteration summary.  Note that
                  this is computationally costly, as it implies computing
                  the poisedness of the current interpolation set.
                  A value of 0 avoids this computation. Note that this option
                  is inactive if verbosity is 'silent'.
                  Default: 0
   restart      : a string defining if a restart (from the information contained
                  the the bcdfo.restart file) is desired:
                  'none' means that restart is not requested (the file bcdfo.restart
                         is then ignored)
                  'use'  means that the algorithm must be restarted from the
                         information saved in bcdfo.restart (an error is generated
                         if the file does not exist)
                  NOTE : Restarted calls MUST specify the same objective function
                         that used in the call at which restart information was
                         saved.
                  Default: 'none'
   save-freq    : an integer defining the frequency at which information is
                  saved (in the file bcdfo.restart) for a possible restart.
                  save-freq = 0 : no restart information is ever saved
                  save-freq > 0 : restart information is saved at termination
                                  of the algorithm
                  Default : 0
   localSolver  : a string defining the local solver to minimize quadratic model
                  'MS' = More-Sorensen (2-norm), 'CG' = truncated CG (inf-norm)
                  Default: 'MS'
   stratLam     : an integer defining the strategy to adjust lambda when using
                  option localSolver = 1 (bound-constrained More-Sorensen)
                  1 = Newtons method, 0 = bisection method
                  Default: 1
   whichmodel   : an integer defining the approach to build the models/Lagr.pols.:
                  0 = Subbasis model (as long as not quadratic)
                  1 = Frobenius-norm model (as long as not quadratic)
                  2 = minimum L2-norm (as long as not quadratic)
                  3 = regression (use option 2 until quadratic, then regression)
                  9 = regression model (use option 2 until quadratic) taking the
                      user-provided gradient information into account
                  Default: 0
   hardcons     : option to consider the given bounds when maximizing the Lagr. pols
                  0 = no bounds on the Lagrange polynomials in maximization
                  1 = take care of bounds in the maximization of the Lagr. pols
                  Default: 0
   noisy        : a boolean to define whether the function is supposed to be noisy
                  0 = no noise expected, 1 = noisy function expected
                  Default: 0
   shrink_Delta : an integer to define whether the trust-region is to shrink in each
                  unsuccessful iteration
                  Default: 0
   Deltamin     : stopping criterion: the minimal reasonable trust-region radius
                  Default: 1e-10
   scaleX       : an integer defining whether a scaling of the variables (using the
                  factors from the vector scalefacX) shall be applied during the
                  optimization process.
                  Default: 0
   scalefacX    : a n-dimensional vector specifying the user-defined scaling factors
                  to be used to scale the variables (they should be chosen such that
                  x(i) * scalefacX(i) = 1 at the solution )
   method       : a string to define the method to use;
                  'BCDFO'   - computes a set of points and builds interpolation
                              models without using derivatives,
                  'BFGS'    - uses a quasi-newton trust-region method after the
                              first function and gradient evaluation,
                  'BFGSDFO' - uses a quasi-newton trust-region method to start
                              with and switches to the derivative-free BCDFO if
                              needed (troubles with noisy gradient,...)
                  --> Please note that in the first and third case, when BCDFO is
                  involved, the parameter "whichmodel" has a relevance and the
                  according model is applied!!
                  Default: 'BCDFO'
 
 
   INPUT (optional and dangerous: don't change this unless you truly know
           what you are doing!):
 
 
   TR_reduction : the mechanism to reduce the TR radius when necessary:
                  'interpolation' uses interpolation
                  'fixed' uses a simple reduction factor (gamma2)
                  Default: 'interpolation'

   mu0          : the initial reduction in gradient norm before recomputing a
                  poised interpolation set
                  Default: 0.1
   mu           : the subsequent reduction in gradient norm before recomputing
                  a poised interpolation set
                  Default: 0
   factor_Dmax  : the multiple of the initial trust-region radius defining the
                  maximal radius for all subsequent iterations
   factor_fmax  : the multiple of the initial function value defining the maximum
                  objective function value for all subsequent iteration (this value
                  is truncated to max( 1.e25, factor_fmax abs(t(f(x0)))))
   kappa_ill    : threshold to declare a system matrix as ill-conditioned
                  Default: 1e+15
   kappa_th     : threshold for a safely nondegenerate set of points
                  Default: 200
   eps_bnd      : distance in which a bound is defined as nearly-active:
                  |x - bnd|<eps_bnd
                  Default: 0.1*epsilon
 
 
   OUTPUT:
 
 
   x            : the best approximation found to a local minimizer,
   fx           : the value of the objective function at x,
   gx           : the estimated gradient of the objective function at x,
   msg          : a short diagnostic message,
   nit          : the number of iterations required by the algorithm,
   neval        : the number of objective function evaluations required by
                  the algorithm.
 
   SOURCES: K. Scheinberg and Ph. L. Toint,
            "Self-correcting geometry in model-based algorithms
            for derivative-free unconstrained optimization ",
            SIAM J. on Opt., 20(6):3512-3532,2010.
 
            S. Gratton, Ph. L. Toint, and A. Troeltzsch,
            "An active-set trust-region method for derivative-free
            nonlinear bound-constrained optimization",
            Opt. Methods and Software, 26(4-5):875-896, 2011.
 
   PROGRAMMING: A. Troeltzsch, Ph. L. Toint, and S. Gratton, 2009-2011.
                (This version 15 XI 2011)
 
   DEPENDENCIES: bcdfo_gradP, bcdfo_repair_Y, bcdfo_save, bcdfo_poisedness,
                 bcdfo_restart, bcdfo_print_summary_vector, bcdfo_print_vector,
                 bcdfo_build_QR_of_Y, bcdfo_augmX_evalf, bcdfo_find_smallf,
                 bcdfo_projgrad, bcdfo_main, bcdfo_computeP, bcdfo_choose_lin
 
   TEST:
   Minimize the ever-famous Rosenbrock "banana valley" function with
   [ x, fx, gx, msg, nit, neval ] = bcdfo( @banana, [-1.2 1] )
   where
      function fx = banana( x )
      fx  = 100 * ( x(2) - x(1)^2 ) ^2 + (1-x(1))^2;
 
   CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
 
   """


###############################################################################
##########################  Initializations ###################################
###############################################################################

   #  Set some constants from class c() in file bcdfo_c.py

   c = bcdfo_c.c()

   #  Initialize to zero
  
   x           = array([[]])
   fx          = 0
   gx          = array([[]])
   nit         = 0
   nitold      = 0
   neval       = 0
   norms       = 0
   oldDelta    = 0
   rho         = 0
   errg        = Inf
   X           = array([[]], dtype='float64')
   Y           = array([[]], dtype='float64')
   fX          = array([], dtype='float64')
   fY          = array([], dtype='float64')
   xstatus     = array([], dtype='int32')
   sstatus     = array([], dtype='int32')
   dstatus     = array([], dtype='int32')
   sspace_save = []
   xspace_save = []
   ndummyY     = 0
   gmat        = array([[]], dtype='float64')
   hessian     = array([[]], dtype='float64')
   indfree     = array([], dtype='int32')
   indfix      = array([], dtype='int32')

   #  Miscellaneous initializations
   
   shapex_1,shapex_2 = shape(x0)
   if ( shapex_1 == 1 and shapex_2 > 1 ):        # make sure the starting point is
      x0 = x0.reshape(shapex_2,1)                      # a column vector

   n       = len( x0 );                        # the dimension of the space
   pquad   = ( ( n + 1 ) * ( n + 2 ) ) / 2;    # the size of a fully quadratic model
   pdiag   = 2 * n + 1;                        # the size of a diagonal model
   plin    = n + 1;                            # the size of a linear model
   xl      = -1e+20 * ones((1,n))              # default lower bounds
   xu      = 1e+20 * ones((1,n))               # default upper bounds
   msg              = 'Unexpected exit';       # no meaningful message at this point
   poisedness_known = 0;                       # the poisedness is not known yet
   eps_rho          = 1.0e-14;                 # stabilization constant for rho
   save_filename    = 'bcdfo.restart';         # the name of the restart file
   savef            = 1;                       # default saving of restart information
   stallfact        = finfo(float).eps         # termination (stalled) when
                                               # ||s|| <= stallfact * norm( x )
   factor_Dmax      = 1.e5;                    # the ratio between Deltamax and Delta0
   factor_fmax      = 1.e20;                   # the ratio between the upper bound on the objective function
                                               # value and the absolute value of t(f(x0))
   CNTsin           = 0;                       # variable to produce a quasi-random vector

   #  Check the parity of the variable argument list.

   noptargs = len( args )                      # the number of optional arguments

   if ( mod( noptargs, 2 ) > 0 ):
      msg   = ' BC-DFO error: the number of variable arguments must be even! ';
      msg   = [ msg, 'Default parameters used.' ];
      print msg;
      noptargs = 0;


   ########################################
   ######## Restart, if requested #########
   ########################################

   #   First look for a restart argument in the calling sequence, while verifying
   #   that every even optional argument is a string.

   restart = 0;                     # default
   lst_noptargs = range(0,noptargs,2)
   for i in lst_noptargs:
      if ( isinstance( args[i], str ) ):
         if (  args[i] == 'restart' ):
            if ( isinstance( args[i+1], str ) ):
               rtype = args[i+1];
               if ( rtype == 'none' ):
                  restart = 0;
               elif ( rtype == 'use' ):
                  restart = 1;
               else:
                  msg = ' BC-DFO error: unknown restart mode! Ignoring.';
                  print msg
            else:
               msg = ' BC-DFO error: unknown restart mode! Ignoring.';
               print msg
      else:
         msg = [' BC-DFO error: argument ', str(i), ' should be a keyword!' ];
         print msg
   
   if ( restart ):
   
      # Read values from file
   
      [ epsilon, eps_current, maxeval, neval, maxit, nitold, verbose, show_errg,   
            rep_degree, eta1, eta2, gamma1, gamma2, gamma3, interpol_TR, factor_CV,
            factor_FPU, factor_FPR, Lambda_XN, Lambda_CP, Lambda_FP, criterion_S, 
            criterion_FP, criterion_CP, mu0, mu, theta, eps_TR, eps_L, shift_Y,  
            Deltamax, fxmax, X, fX, Delta0, lSolver, stratLam, whichmodel, kappa_ill, 
            kappa_th, eps_bnd, i_xbest, hardcons, noisy, vstatus, dstatus, 
            sspace_save, xspace_save, scaleX, scalefacX, shrink_Delta, gmat, 
            useMethod, current_method, Deltamin, restok ] = bcdfo_restart( save_filename )
   
      if ( ~restok ):
         msg = [' BC-DFO error: could not open ', save_filename,'. Starting afresh from supplied x0!!'];
         print msg
         restart = 0;
      else:
         poised_model = 0;
         nit = nitold;
      
   if (  not restart ):
   
      #  Set defaults
   
      Delta0       = 1;                # the initial TR radius
      Deltamin     = 1.0e-10;          # termination if minimum Delta reached
      cur_degree   = plin;             # the degree for the initial model
      rep_degree   = plin;             # the minimum degree for the model after repair
      epsilon      = 1.0e-5;           # gradient termination accuracy
      maxeval      = 200 * n;          # maximum number of evaluations
      maxit        = maxeval;          # maximum number of iterations
      verbose      = 1;                # printout quantity
      show_errg    = 0;                # display of the gradient error estimate
      initial_Y    = 'simplx';         # geometry of the initial interpolation set
      eta1         = 0.0001;           # min rho value for successful iterations
      eta2         = 0.9;              # min rho value for very successful iterations
      gamma1       = 0.01;             # lower bound on Delta interval (unsuccessful its)
      gamma2       = 0.5 ;             # upper bound on Delta interval (unsuccessful its)
      gamma3       = 2.0;              # increase in Delta (successful its)
      interpol_TR  = 1;                # use interpolation for the radius of the trust region
      factor_CV    = 10;               # constant for termination (see above)
      Lambda_XN    = 1.0e-10;          # poisedness for new iterates
      Lambda_CP    = 1.2;              # poisedness for close points
      Lambda_FP    = 1.0e-10;          # poisedness for far points
      factor_FPU   = 1;                # multiple of TR radius defining far points
                                       # (for unsuccessful iterations)
      factor_FPR   = 10;               # multiple of TR radius defining far points
                                       # (for reconstruction of a poised interpolation set
                                       # for large model gradient)
      criterion_S  = 'distance';       # selection of outgoing point at successful iterations:
      criterion_FP = 'distance';       # the same, but for far points at unsuccessful iterations
      criterion_CP = 'standard';       # the same, but for close points at unsuccessful iterations
      mu0          = 0;                # initial reduction in gradient norm before repair
      mu           = 0;                # subsequent reduction in gradient norm before repair
      theta        = 1;                # ratio between gradient norm and radius after repair
      eps_TR       = 0.0001;           # rel. accuracy on the trust-region constraint for steps
      eps_L        = 0.001;            # rel. accuracy on the trust-region constraint for L max
      shift_Y      = 1;                # shifting and scaling of the set Y is used
      lSolver      = 1;                # local solver for minimizing the model
      stratLam     = 1;                # strategy to adjust lambda when solving the bc MS problem,
      kappa_ill    = 1e+15;            # threshold to declare a system matrix as ill-conditioned
      kappa_th     = 2000;             # threshold for a safely nondegenerate set of points
      eps_bnd      = epsilon/10;       # epsilon to define a bound as nearly-active: |x - bnd|<eps_bnd
      whichmodel   = 0;                # approach to build the local models
      hardcons     = 0;                # apply hard bounds when maximizing the Lagrange polynomials
      noisy        = 0;                # function supposed to be noisy
      scaleX       = 0;                # scaling of variables is applied
      scalefacX    = ones((1,n),float) # scaling factors initialized to one
      shrink_Delta = 0;                # shrink trust-region radius in every unsuccessful iteration
      useMethod    = 0;                # defines the method to use (BCDFO, BFGS)    
   
   ########################################
   ## Process the variable argument list ##
   ########################################
   
   for i in lst_noptargs:
   
      #  The initial trust-region radius
   
      if ( args[i] == 'Delta0' ):
         if ( isinstance( args[i+1], (int,float) ) ):
            var = args[i+1];
            if ( var > 0 ):
               Delta0 = var;
            else:
               msg = ' BC-DFO error: negative value for parameter Delta0! Default used.';
               print msg
            
         else:
            msg = ' BC-DFO error: wrong type of input for parameter Delta0! Default used.';
            print msg
         
      #  The minimal trust-region radius
   
      elif ( args[i] == 'Deltamin' ):
         if ( isinstance( args[i+1], (int,float) ) ):
            var = args[i+1];
            if ( var > 0 ):
               Deltamin = var;
            else:
               msg = ' BC-DFO error: negative value for parameter Deltamin! Default used.';
               print msg
             
         else:
            msg = ' BC-DFO error: wrong type of input for parameter Deltamin! Default used.';
            print msg
          
      #  lower bounds on the variables
         
      elif ( args[i] == 'lbounds' ):
         if ( len(args[i+1][0]) == n ):
            xl = args[i+1]
            if (size(xl,1)==1):
               xl = transpose(xl)
         else:
            msg = ' BC-DFO error: lower bounds empty or not of length n! No bounds used.';
            print msg
          
      #  upper bounds on the variables
         
      elif ( args[i] == 'ubounds' ):
         if ( len(args[i+1][0]) == n ):
            xu = args[i+1]
            if (size(xu,1)==1):
               xu = transpose(xu)
         else:
            msg = ' BC-DFO error: upper bounds empty or not of length n! No bounds used.';
            print msg

      #  The type of initial interpolation model

      elif ( args[i] == 'degree_m0' ):
         if ( isinstance( args[i+1], str ) ):
            stype = args[i+1];
            if ( stype == 'linear' ):
               cur_degree = plin;
            elif ( stype == 'diagonal' ):
               cur_degree = pdiag;
            elif ( stype == 'quadratic' ):
               cur_degree = pquad;
            else:  
               msg = ' BC-DFO error: unknown type of initial interpolation model! Default used.';
               print msg
         else:
            msg = ' BC-DFO error: wrong type of input for initial interpolation model! Default used.';
            print msg

      #  The type of interpolation model on model reconstruction

      elif ( args[i] == 'degree_mR' ):
         if ( isinstance( args[i+1], str ) ):
            stype = args[i+1];
            if ( stype == 'linear' ):
               rep_degree = n + 1;
            elif ( stype == 'diagonal' ):
               rep_degree = 2 * n + 1;
            elif ( stype == 'quadratic' ):
               rep_degree = pquad;
            else:
               msg = ' BC-DFO error: unknown type of interpolation model on recomputation! Default used.';
               print msg
         else:
            msg = ' BC-DFO error: unknown type of interpolation model on recomputation! Default used.';
            print msg

      #  The termination accuracy

      elif ( args[i] == 'epsilon' ):
         if ( isinstance( args[i+1], (float,int) ) ):
            if ( args[i+1] > 0 ):
               epsilon = args[i+1];
            else:
               msg = ' BC-DFO error: negative value for parameter epsilon! Default used.';
               print msg
         else:
            msg = ' BC-DFO error: wrong type of input for parameter epsilon! Default used.';
            print msg

      #  The maximum number of objective function's evaluations
   
      elif ( args[i] == 'maxeval' ):
         if ( isinstance( args[i+1], int ) ):
            maxeval = abs( args[i+1]);
         else:
            msg = ' BC-DFO error: wrong type of input for parameter maxeval! Default used.';
            print msg

      #  The maximum number of iterations
   
      elif ( args[i] == 'maxit' ):
         if ( isinstance( args[i+1], int ) ):
            maxit = abs( args[i+1]);
         else:
            msg = ' BC-DFO error: wrong type of input for parameter maxit! Default used.';
            print msg

      #  The printout ouput

      elif ( args[i] == 'verbosity' ):
         if ( isinstance( args[i+1], str ) ):
            vtype = args[i+1];
            if ( vtype == 'silent' ):
               verbose = 0;
            elif ( vtype == 'low' ):
               verbose = 1;
            elif ( vtype == 'medium' ):
               verbose = 2;
            elif ( vtype == 'high' ):
               verbose = 3;
            elif ( vtype == 'debug' ):
               verbose = 10;
            else:
               msg = ' BC-DFO error: unknown verbosity level! Default used.';
               print msg
         else:
            msg = ' BC-DFO error:  unknown verbosity level! Default used.';
            print msg

      #  The gradient error printout

      elif ( args[i] == 'show_errg' ):
         if ( isinstance( args[i+1], int ) ):
            show_errg = max( min( 1, args[i+1] ), 0 );
         else:
            msg = ' BC-DFO error: wrong type of input for parameter show_errg! Default used.';
            print msg

      #  The type of initial interpolation set

      elif ( args[i] == 'geometry_Y0' ):
         if ( isinstance( args[i+1], str ) ):
            mtype = args[i+1];
            if ( mtype == 'random'  ):
               initial_Y = 'random';
            elif ( mtype == 'simplex'  ):
               initial_Y = 'simplx';
            else:
               msg = ' BC-DFO error: unknown type of initial interpolation set! Default used.';
               print msg
         else:
            msg = ' BC-DFO error: unknown type of initial interpolation set! Default used.';
            print msg

      #  The mechanism for radius reduction

      elif ( args[i] == 'TR_reduction' ):
         if ( isinstance( args[i+1], str ) ):
            Dtype = args[i+1];
            if ( Dtype == 'interpolation' ):
               interpol_TR = 1;
            elif ( Dtype == 'fixed' ):
               interpol_TR = 0;
            else:
               msg = ' BC-DFO error: unknown mechanism for radius reduction! Default used.';
               print msg
         else:
            msg = ' BC-DFO error: unknown mechanism for radius reduction! Default used.';
            print msg

      # The initial gradient reduction before repair

      elif ( args[i] == 'mu0' ):
         if ( isinstance( args[i+1], float ) ):
            mu0 = abs( args[i+1]);
         else:
            msg = ' BC-DFO error: wrong type of input for parameter mu0! Default used.';
            print msg
       
      # The subsequent gradient reduction before repair

      elif ( args[i] == 'mu' ):
         if ( isinstance( args[i+1], float ) ):
            mu = abs( args[i+1]);
         else:
            msg = ' BC-DFO error: wrong type of input for parameter mu! Default used.';
            print msg

      # The ratio between the initial TR radius and the maximal one

      elif ( args[i] == 'factor_Dmax' ):
         if ( isinstance( args[i+1], (float, int) ) ):
            factor_Dmax = abs( args[i+1]);
         else:
            msg = ' BC-DFO error: wrong type of input for parameter factor_Dmax! Default used.';
            print msg
       
      # The ratio between the initial objective function value and the maximal one

      elif ( args[i] == 'factor_fmax' ):
         if ( isinstance( args[i+1], (float, int) ) ):
            factor_fmax = abs( args[i+1]);
         else:
            msg = ' BC-DFO error: wrong type of input for parameter factor_fmax! Default used.';
            print msg
    
      # The local solver to minimize the model

      elif ( args[i] == 'localSolver' ):
         if ( isinstance( args[i+1], str ) ):
            ctype = args[i+1];
            if ( ctype == 'MS' ):
               lSolver = 1;
            elif ( ctype == 'CG' ):
               lSolver = 2;
            else:
               msg = ' BC-DFO error: unknown value for parameter localSolver! Default used.';
               print msg
         else:
            msg = ' BC-DFO error: wrong type of input for parameter localSolver! Default used.';
            print msg
      
      # The strategy to find a new point in the bound-constrained trust-region solver

      elif ( args[i] == 'stratLam' ):
         if ( isinstance( args[i+1], int ) ):
            if ( args[i+1] > 1 or args[i+1] < 0 ):
               msg = ' BC-DFO error: unknown value for parameter stratLam! Default used.';
               print msg
            else:
               stratLam = args[i+1];
         else:
            msg = ' BC-DFO error: wrong type of input for parameter stratLam! Default used.';
            print msg

      # The threshold for an ill-conditioned system matrix

      elif ( args[i] == 'kappa_ill' ):
         if ( isinstance( args[i+1], (float, int) ) ):
            if ( args[i+1] < 10 ):
               msg = ' BC-DFO error: meaningless value for parameter kappa_ill! Default used.';
               print msg
            else:
               kappa_ill = args[i+1];
         else:
            msg = ' BC-DFO error: wrong type of input for parameter kappa_ill! Default used.';
            print msg

      # The threshold for a safely nondegenerate set of points

      elif ( args[i] == 'kappa_th' ):
         if ( isinstance( args[i+1], (float, int) ) ):
            if ( args[i+1] < 10 ):
               msg = ' BC-DFO error: meaningless value for parameter kappa_th! Default used.';
               print msg
            else:
               kappa_th = args[i+1];
         else:
            msg = ' BC-DFO error: wrong type of input for parameter kappa_th! Default used.';
            print msg
      
      # The epsilon to define a bound as nearly-active: |x - bnd| < eps_bnd
   
      elif ( args[i] == 'eps_bnd' ):
         if ( isinstance( args[i+1], (int,float) ) ):
            eps_bnd = abs( args[i+1]);
         else:
            msg = ' BC-DFO error: wrong type of input for parameter eps_bnd! Default used.';
            print msg

      # The model building method

      elif ( args[i] == 'whichmodel' ):
         if ( isinstance( args[i+1], int ) ):
            whichmodel = abs( args[i+1] );
         else:
            msg = ' BC-DFO error: wrong type of input for parameter whichmodel! Default used.';
            print msg

      # The information saving frequency

      elif ( args[i] == 'save-freq' ):
         if ( isinstance( args[i+1], int ) ):
            savef = max( min( 1, args[i+1]), 0 );
         else:
            msg = ' BC-DFO error: wrong type of input for parameter save-freq! Default used.';
            print msg

      # hard constraints

      elif ( args[i] == 'hardcons' ):
         if ( isinstance( args[i+1], int ) ):
            hardcons = min( max( args[i+1], 0 ), 1 );
         else:
            msg = ' BC-DFO error: wrong type of input for parameter hardcons! Default used.';
            print msg

      # noisy function expected

      elif ( args[i] == 'noisy' ):
         if ( isinstance( args[i+1], int ) ):
            noisy = min( max( args[i+1], 0 ), 1 );
         else:
            msg = ' BC-DFO error: wrong type of input for parameter noisy! Default used.';
            print msg
   
      # shrink trust-region at each unsuccessful iteration
   
      elif ( args[i] == 'shrink_Delta' ):
         if ( isinstance( args[i+1], int ) ):
            shrink_Delta = min( max( args[i+1], 0 ), 1 );
         else:
            msg = ' BC-DFO error: wrong type of input for parameter shrink_Delta! Default used.';
            print msg
   
      # define the method to use
   
      elif ( args[i] == 'method' ):
         if ( isinstance( args[i+1], str ) ):
            ctype = args[i+1];
            if ( ctype == 'BCDFO' ):
               useMethod = 0;
            elif ( ctype == 'BFGS' ):
               useMethod = 1;
            elif ( ctype == 'BFGSDFO' ):
               useMethod = 2;
            else:
               msg = ' BC-DFO error: unknown value for parameter "method"! Default used.';
               print msg
         else:
            msg = ' BC-DFO error: wrong type of input for parameter "method"! Default used.';
            print msg

      # apply variables scaling

      elif ( args[i] == 'scaleX' ):
         if ( isinstance( args[i+1], int ) ):
            scaleX = min( max( args[i+1], 0 ), 1 );
         else:
            msg = ' BC-DFO error: wrong type of input for parameter scaleX! Default used.';
            print msg
    
      # variables scaling factors

      elif ( args[i] == 'scalefacX' ):
            if ( len(args[i+1]) == n ):
               scalefacX = args[i+1];
               sizefac1, sizefac2 = shape( scalefacX )
               if ( sizefac1 == 1 and sizefac2 > 1 ):    # make sure the scaling array is
                  scalefacX = transpose( scalefacX );                    # a column vector
            else:
               msg = ' BC-DFO error: the vector scalefacX is not of length n! No scaling used.';
               print msg;
   
      # The restart mode (treated already)

      elif ( args[i] == 'restart' ):
         """ nothing """
      
      # Unidentified keyword

      else:
         msg = [ ' BC-DFO warning: undefined keyword ', args[i], '! Ignoring.' ];
         print msg

   #  Compute the maximal TR radius.

   if ( not restart ):
      Deltamax = factor_Dmax * Delta0;

   ###############################################################################
   #########################  Start the algorithm proper #########################
   ###############################################################################

   #  If the algorithm is not restarted from a saved configuration, the initial
   #  interpolation set must be computed.

   if ( not restart ):

      #  Checking the bounds and correct Delta0 if there is insufficient space
      #  between the bounds. Modification of x0 if construction of first
      #  interpolation model would interfere with given bounds.
   
      zero    = 0.0;
      nfix    = 0;
      xfix    = zeros((n,1),float);
      vstatus = zeros((n),int)
      
      #  Define the current method to use for the starting sequence
   
      if ( useMethod == 0 ):
         current_method = 0;       # BCDFO
      else:
         current_method = 1;       # BFGS or BFGSDFO

      #  Safeguard bounds and the starting point

      for j in range(0,n):

         #  Check lower and upper bounds.

         if (xl[0,j] > xu[0,j]):
            msg  = ['Error: Lower bound of component ' +
               str( j ) + ' exceeds upper bound !!'];
            print msg
            return x, fx, gx, msg, nit, neval
       
         #  Check difference between bounds.
       
         temp = xu[0,j] - xl[0,j];

         if (temp < Delta0 + Delta0):
            if (temp == zero):
               nfix       = nfix + 1
               indfix = append(indfix,j)
               vstatus[j] = c.alwaysfixed()
               xfix[0,j]  = xl[0,j]
            else:
               Delta0 = 0.5 * temp
               print ' Diff. between lower and upper bound of component '+str( j )+' is less than 2*Delta0 !! New Delta0='+str(Delta0)  

         #  Move the starting point inside the bounds if necessary

         templ = xl[0,j] - x0[j,0]
         tempu = xu[0,j] - x0[j,0]

         if (templ >= -Delta0):
           x0[j,0] = xl[0,j] + Delta0;
         elif (tempu <= Delta0):
           x0[j,0] = xu[0,j] - Delta0;

      #  Scale x0 and bounds if user-defined.

      if ( scaleX ):
         for i in range(0,n):
            if ( scalefacX[i] > 0 ):
               x0[0,i] = x0[0,i] * scalefacX[i];
               xl[0,i] = xl[0,i] * scalefacX[i];
               xu[0,i] = xu[0,i] * scalefacX[i];
            else:
               scalefacX[i] = 1;      

      #  Reset constants if fixed some variables.

      if (nfix > 0):
         nfree   = n - nfix
         if ( nfree <= 0 ):
            msg = 'No free variables. Please, enlarge search space!'
            print msg
            return x, fx, gx, msg, nit, neval
            
         indfull = list(range(0,n))
         indfree = list(set(indfull)-set(indfix))

         x0      = x0[indfree,:]

         if (cur_degree == plin):
            cur_degree = nfree + 1;
         elif (cur_degree == pdiag):
            cur_degree = 2*nfree+1; 
         elif (cur_degree == pquad):
            cur_degree = ((nfree + 1) * (nfree + 2)) / 2; 
       
         if (rep_degree == plin):
            rep_degree = nfree+1; 
         elif (rep_degree == pdiag):
            rep_degree = 2*nfree+1; 
         elif (rep_degree == pquad):
            rep_degree = ((nfree + 1) * (nfree + 2)) / 2; 
       
         plin    = nfree + 1;
         pdiag   = 2 * nfree + 1;
         pquad   = ((nfree + 1) * (nfree + 2)) / 2;
         n       = nfree;
      else:
         indfree = list(range(0, n))
         
      #  Include x0 as the first point in the interpolation set.

      Y = x0

      #  Differ between BCDFO and BFGS.

      if ( current_method == 0 ):        # BCDFO

         getfY = 1;
         while ( getfY ):
            if ( verbose >= 2 ):
               print ' Degree of the initial  model = '+str( cur_degree )

            #  Compute an initial poised interpolation set around the starting point.
 
            if ( initial_Y == 'random' ):
            
               #  Loop in case of an accidentally ill-conditioned initial system
         
               ill_init = 1;
               while ( ill_init ):

                  Y = append(Y, - ones((n,cur_degree-1),float) + 2 * random.rand(n, cur_degree-1), axis=1)

                  lst_cur_degree = range(1,cur_degree)
                  for  j in lst_cur_degree:
                     Y[:,j] = Y[:,0] + Y[:,j] * ( Delta0 / linalg.norm( Y[:,j] ) );

                  #  Build an initial factorization.

                  [ QZ, RZ, x, scale ] = bcdfo_build_QR_of_Y( Y, whichmodel, 
                     shift_Y, Delta0, 1, kappa_ill );

                  # check the condition of the interpolation matrix

                  ill_init = bcdfo_checkZ( RZ, kappa_ill )
             
               #  Make the set poised.

               [ QZ, RZ, Y, replaced, poised, Y_radius, x, scale ] = bcdfo_repair_Y(
                  QZ, RZ, Y, Delta0, factor_FPR,
                  Lambda_FP, Lambda_CP, eps_L, x, lSolver, whichmodel, hardcons, 
                  xl, xu, indfree, stratLam, scale, shift_Y, 1, kappa_ill );
               poisedness_known = 1;

            elif ( initial_Y == 'simplx' ):

               #  Compute the initial interpolation set (simplex plus midpoints).

               I = eye( n );

               # lin
               for j in range(0,n):
                  step1 = -Delta0;
                  Y = append(Y, x0 + step1 * I[:,j].reshape((n,1)), axis=1);
               
               # diag
               if ( cur_degree >= pdiag ) :
                  for j in range(0,n):
                     step2 =  Delta0;
                     Y = append(Y, x0 + step2 * I[:,j].reshape((n,1)), axis=1);

               # quad
               if ( cur_degree == pquad ):
                  for j in range(0,n):
                     for jj in range(j+1,n):
                        Y = append(Y, 0.5*(Y[:,j+1].reshape((n,1))
                        + Y[:,jj+1].reshape((n,1))), axis=1)

               #  Build the initial factorization.

               QZ, RZ, x, scale = bcdfo_build_QR_of_Y( Y, whichmodel, shift_Y, Delta0, 1.0, kappa_ill )

               x = copy(Y[:,0])

               #  Compute poisedness of initial interpolation set

               [ poised, Y_radius ] = bcdfo_poisedness_Y( QZ, RZ, Y, eps_L, x,
                  1, whichmodel, hardcons, xl, xu, indfree, stratLam, scale,
                  shift_Y )
 
               poisedness_known     = 1

            # The initial interpolation set is poised.
            
            poised_model = 1                   

            #  Compute the associated function values, possibly reducing Delta0 to
            #  ensure that the objective function remains finite at all interpolation
            #  points.

            for i in range(0,cur_degree):

               #  Storing the new points in X and evaluate f

               if ( whichmodel >= 8 ):
                  compute_g = 1;
               else:
                  compute_g = 0;

               [X, fX, neval, xstatus, sstatus, dstatus, gmat, msg] = bcdfo_augmX_evalf(
                  f, Y[:,i], i, X, fX, nfix, xfix, indfix,
                  indfree, 1.e25, neval, maxeval, verbose, show_errg, xstatus, 
                  c.inY(), sstatus, dstatus, scaleX, scalefacX, compute_g, gmat );

               fY = append( fY, fX[i] )
               
               if ( msg[0:5] == 'Error' ):
                  print msg;
                  #  including fixed variables at return
                  if ( nfix > 0 ):
                     I  = eye( n + nfix )
                     x  = I[:,indfix] * xl[indfix] + I[:,indfree] * x
                
                  return x, fx, gx, msg, nit, neval             

               #  If the computed function value is infinite, restart with a smaller Delta0.

               if ( abs( fY[i] ) > 1.0e25 ):
                  break             

               #  All functions values at points of Y are finite.  No need for
               #  another pass.

               if ( i+1 == cur_degree ):
                  getfY = 0;             
        
            #  Another pass is needed with a smaller Delta0 (at least one function
            #  value is infinite).

            if ( getfY ):

               #  Decrease Delta0.

               Delta0 = gamma1 * Delta0;

               #  Terminate with an error message if Delta0 becomes negligible wrt x0.
   
               if ( Delta0 < stallfact * norm( x0 ) ):
                  msg = 'Error: cannot find enough finite objective function values in the neighbourhood of the starting point! Terminating.'
                  print msg

                  #  including fixed variables at return
                  if (nfix>0):
                     I  = eye( n + nfix );
                     x  = I[:,indfix] * xl[indfix] + I[:,indfree] * x
                     gx = I[:,indfix] * zeros( nfix, 1 ) + I[:,indfree] * gx
                
                  return x, fx, gx, msg, nit, neval

             # while-loop

         fx0 = fY[0]
         m   = copy(cur_degree)-1
         ind_Y = arange(0,cur_degree)
   
         #  Move to the best point in the interpolation set, if different from x0.

         [x, fx, QZ, RZ, Y, fY, ind_Y, i_xbest, scale] = bcdfo_find_smallf(c, QZ, RZ, Y,
            fY, ind_Y, 0, cur_degree, indfree, x, xl, xu, fx0,
            dstatus, whichmodel, scale, shift_Y, Delta0, 1, kappa_ill);
         x = copy(Y[:,0])

      else:     # ( useMethod = BFGS or BFGSDFO, thus start with BFGS-method )

         compute_g = 1;
         [X, fX, neval, xstatus, sstatus, dstatus, gmat, msg] = bcdfo_augmX_evalf(
            f, x0, 0, X, fX, nfix, xfix, indfix,
            indfree, 1.e25, neval, maxeval, verbose, show_errg, xstatus, 
            c.inY, sstatus, dstatus, scaleX, scalefacX, compute_g, gmat );
         fY = append( fY, fX[0])

         if ( strcmp( msg[0:5], 'Error' ) ):
            print msg;
            x = x0;
            #  including fixed variables at return
            if ( nfix > 0):
               I  = eye( n + nfix )
               x  = I[:,indfix] * xl[indfix] + I[:,indfree] * x
             
            return x, fx, gx, msg, nit, neval       

         #  If the computed function value is infinite, stop.

         if ( abs( fY[1] ) > 1.0e25 ):
            msg = 'Error: absolute value of function at initial point > 10^25! Setup a new bound.'
            print msg
            return x, fx, gx, msg, nit, neval
       
         fx0     = copy(fY[0]);
         Y       = copy(X[indfree,0]);
         ind_Y   = [0];
         i_xbest = 0;
         x       = x0[indfree,0];
         fx      = fx0;
         m       = 0;
         QZ = [];
         RZ = [];
         model = [];
         poised_model = 0;
         poised = 0;
         Y_radius = 0;
         scale = [];
         hessian = eye(n);

         # end of useMethod

      #  Compute the maximal objective function value.

      fxmax = min( 1.e25, factor_fmax * abs( fx ) );

   else:

      #  If the algorithm is restarted from a saved configuration, then the interpolation
      #  set is already read from this configuration, but its factorization must still be
      #  recomputed.

      m       = size(X,2)-1
      sstatus = transpose( ones(m+1,1) )
      maxeval = maxeval + neval
      maxit   = maxit + nitold

      #  Look for fixed variables in the problem and reset the parameters appropriately.

      xfix               = zeros(n,1)
      indfix             = where( vstatus == 2 )

      if ( isempty( indfix ) ):
         nfix            = 0;
         indfree         = range(0,n)
      elif ( n == len( indfix ) ):
         msg = 'No free variables. Please, enlarge search space!';
         print msg
         return x, fx, gx, msg, nit, neval
      else:
         nfix            = len( indfix );
         xfix[indfix]    = xl[indfix];
         nfree           = n - nfix;
         indfree         = setdiff( range(0,n), indfix );
         if ( rep_degree == plin ):
            rep_degree   = nfree + 1;
         elif ( rep_degree == pdiag ):
            rep_degree   = 2 * nfree + 1;
         elif ( rep_degree == pquad ):
            rep_degree   = ((nfree + 1) * (nfree + 2)) / 2;
       
         plin            = nfree + 1;
         pdiag           = 2 * nfree + 1;
         pquad           = ((nfree + 1) * (nfree + 2)) / 2;
         n               = nfree;
    
      if ( current_method == 0 ):   # BCDFO

         #  Select an interpolation set Y from the set of points X
         #  (because we don't know if configuration was savedin the full space
         #  or in the subspace)

         [ i_xbest, m, X, fX, cur_degree, ind_Y, Y, neval, xstatus, sstatus, dstatus,
            ndummyY, CNTsin, scale, gmat, msg ] = bcdfo_choose_lin(
            f, i_xbest, xl, xu, X, fX, fxmax, Delta0, kappa_th, 
            eps_L, rep_degree, lSolver, vstatus, sstatus, dstatus, xfix, neval, 
            maxeval, ones(1,n+1), shift_Y, c, 'fullspace', plin, pdiag, pquad, whichmodel, 
            stratLam, hardcons, show_errg, 1, kappa_ill, CNTsin, scaleX, scalefacX, gmat );
         fY = fX[ind_Y];
         fx = fX[i_xbest];

         [ QZ, RZ, x, scale ] = bcdfo_build_QR_of_Y( Y, whichmodel, shift_Y, 
            Delta0, 1, kappa_ill );
         [ poised, Y_radius ] = bcdfo_poisedness_Y( QZ, RZ, Y, eps_L, x, 
            lSolver, whichmodel, hardcons, xl, xu, indfree, stratLam, scale, shift_Y );
         poisedness_known     = 1;

      else:    # BFGS

         fY = fX[i_xbest];
         fx = fX[i_xbest];
         poisedness_known = 0;
         poised = 0;
         Y_radius = 0;
         QZ = [];
         RZ = [];
         x = X[indfree,i_xbest];
         hessian = eye(n);
         ind_Y = [i_xbest];
         cur_degree = plin;
         model = [];
         scale = [];
         xstatus = transpose( zeros((1,m+1)) )
         xstatus[i_xbest] = 1;

         # end of current_method

      # end of restart or not

   # Initialize Delta

   Delta = Delta0;

   #  Compute the associated polynomial interpolation model.

   if ( current_method == 0 ):      # BCDFO
      
      initmodel = zeros((1,pquad),float)
      model  = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, initmodel, ind_Y, 
         0, 0, gx, scale, shift_Y, Delta, indfree, gmat, eye(n), epsilon, noisy )

      gx     = bcdfo_gradP( model, x, x, scale, shift_Y )
      normgx, pgx = bcdfo_projgrad( n, x, gx, xl[0,indfree], xu[0,indfree] )

   else:                            # BFGS

      gx = gmat[indfree,i_xbest]
      normgx, pgx = bcdfo_projgrad( n, x, gx, xl[0,indfree] ,xu[0,indfree] )
   
   if ( not restart ):

      #  Initialize the current epsilon
   
      eps_current = max( mu0 * normgx, epsilon )

   if ( verbose >= 10 ):
      bcdfo_print_vector( 'm0', model )

   ########################################
   ########  Initial printout  ############
   ########################################

   if ( verbose >= 1 ):

      #  Print the banner and header.

      if ( show_errg ):
         print( '\n')
         print( '  ***********************************************************************************************')
         print( '  *                                                                                             *')
         print( '  *             BC-DFO: bound-constrained minimization with/without derivatives                 *')
         print( '  *                                                                                             *')
         print( '  *                (c) A. Troeltzsch, Ph. L. Toint and S. Gratton, 2009-2011                    *')
         print( '  *                                                                                             *')
         print( '  ***********************************************************************************************')
         print( '')
         print( '  it    neval        fvalue           ||g||     errg      ||s||     Delta       rho    nfree type');
      else:
         print( '\n')
         print( '  *************************************************************************************')
         print( '  *                                                                                   *')
         print( '  *        BC-DFO: bound-constrained minimization with/without derivatives            *')
         print( '  *                                                                                   *')
         print( '  *           (c) A. Troeltzsch, Ph. L. Toint and S. Gratton, 2009-2011               *')
         print( '  *                                                                                   *')
         print( '  *************************************************************************************')
         print( '')
         print( '  it    neval        fvalue           ||g||     ||s||     Delta       rho    nfree type');
    
      #  Print the objective function value at the starting point.
      
      if ( ( not restart ) and ( current_method == 0 ) ):
         print '%5d  %5d  %+.14e' % (-1, 1, fx0)
#         print '{0:5d} {1:6d} {2:18f}'.format(-1,1,fx0)
         if ( verbose == 3 ):
            bcdfo_print_summary_vector( 'x0', x0 )
         elif ( verbose > 3 ):
            bcdfo_print_vector( 'x0', x0 )
          
      #  Print the information obtained after computing the first model.

      if ( show_errg ):

         #  Compute poisedness and gradient error, if unavailable so far.

         if ( not poisedness_known ):
            [ poised, Y_radius ] = bcdfo_poisedness_Y( QZ, RZ, Y, eps_L, x, lSolver, 
                      whichmodel, hardcons, xl, xu, indfree, stratLam, scale, shift_Y );
            poisedness_known     = 1;
       
         errg = poised * Y_radius / factor_CV;

         print '%5d  %5d  %+.14e  %.2e  %.2e' % (nitold, neval, fx, normgx, errg )

      else:
      
         print '%5d  %5d  %+.14e  %.2e' % (nitold, neval, fx, normgx)

      # printout to file convhist.m

#      if ( verbose ):
#         fid = fopen('convhist.m','w');
#         print(fid,'function A=history \n A=[ \n');
#         print(fid,'%6d  %+.14e %.2e \n', neval, fx, normgx);
#         fclose(fid);    

      if ( verbose == 3 ):
         bcdfo_print_summary_vector( 'x ', x )
         bcdfo_print_summary_vector( 'gx', gx )
      elif ( verbose > 3 ):
         bcdfo_print_vector( 'x ', x )
         bcdfo_print_vector( 'gx', gx ) 

   ###############################################################################
   #######################  CALL MAIN ROUTINE  ###################################
   ###############################################################################

   #pquad = pquad - n
   hessian_bfgs = []
   g_bfgs = []
   stage = 0

   [nit, i_xbest, x, fx, m, X, fX, ind_Y, Delta, eps_current, cur_degree, model, 
      gx, normgx, vstatus, xstatus, sstatus, dstatus, ndummyY, sspace_save,
      xspace_save, msg, CNTsin, gmat, hessian, hessian_bfgs, stage, neval] = bcdfo_main(
      f, nitold, nit, i_xbest, xl, xu, m, X, fX, ind_Y, QZ, RZ, Delta,
      cur_degree, neval, maxeval, maxit, model, gx, normgx, show_errg, pquad, pdiag, 
      plin, stallfact, eps_rho, Deltamax, rep_degree, epsilon, verbose, eta1, eta2, 
      gamma1, gamma2, gamma3, interpol_TR, factor_CV, Lambda_XN, Lambda_CP, factor_FPU,
      factor_FPR, Lambda_FP, criterion_S, criterion_FP, criterion_CP, mu, theta, 
      eps_TR, eps_L, lSolver, stratLam, eps_current, vstatus, xstatus, transpose(sstatus),
      dstatus, ndummyY, sspace_save, xspace_save, xfix, fxmax, poised_model, 
      kappa_ill, kappa_th, eps_bnd, poised, Y_radius, c, 'toplevel', whichmodel, 
      hardcons, noisy, scaleX, scalefacX, CNTsin, shrink_Delta, gmat, hessian, 
      hessian_bfgs, g_bfgs, useMethod, stage, Deltamin, scale, shift_Y)

   ###############################################################################
   #######################  RETURN FROM MAIN  ####################################
   ###############################################################################

#   if ( savef > 0 ):
#      savok = bcdfo_save( save_filename, epsilon, eps_current, maxeval, neval,   
#                    maxit, nit, verbose, show_errg, rep_degree, eta1, eta2,    
#                    gamma1, gamma2, gamma3, interpol_TR, factor_CV, factor_FPU,
#                    factor_FPR, Lambda_XN, Lambda_CP, Lambda_FP, criterion_S,  
#                    criterion_FP, criterion_CP, mu0, mu, theta, eps_TR,        
#                    eps_L, shift_Y, X, fX, Delta, Deltamax, fxmax, lSolver,    
#                    stratLam, whichmodel, kappa_ill, kappa_th, eps_bnd,        
#                    i_xbest, hardcons, noisy, vstatus, dstatus, sspace_save,   
#                    xspace_save, scaleX, scalefacX, shrink_Delta, gmat,        
#                    useMethod, current_method, Deltamin );
#      if ( not savok ):
#         print [ '  BC-DFO saving error: could not save on file ', save_filename ] 
    
 

   #  Add closing bracket in file convhist.m

#   if ( verbose ):
#      fid = fopen('convhist.m','a');
#      print(fid,'];');
#      fclose(fid);

   #  Re-assemble gradient at return

   if ( nfix > 0 ):
      I  = eye( n + nfix );
      x  = I[:,indfix] * xfix[indfix] + I[:,indfree] * x;
      gx = I[:,indfix] * zeros(nfix,1) + I[:,indfree] * gx;

   #  Rescale x if necessary

   if ( scaleX ):
      x = x / scalefacX;   # verify elementwise division!

   #  Return to the calling routine

   return x, fx, gx, msg, nit, neval

