#!/usr/local/bin/python

from bcdfo import *
from bcdfo_functions import *
from funlist import *


#fList = [ "linquad","Aquadratic", "quadratic", "rosenbrock", "cubic", "Rosenbrock2D"]
fList = [ "Rosenbrock2D" ]
nprob = size(fList)

for f_name in fList:
   funf = eval(f_name)
   dfunf = eval("df"+f_name)
   print "Function ",f_name
   
   #  Initialize test problem

   #   [ x0, xl, xu, nprob ] = funlist_init( funf );

   if ( f_name == 'Rosenbrock2D' ):
      n = 2
      xl = numpy.array([[-1., -1.]])
      xu = numpy.array([[9., 9.]])
      x0 = numpy.array([[2., 5.]])
      print xl, xu, x0
   else:
      n  = 4
      x0 = -0.5 * ones( (1,n), float )
      xl = -10 * ones( (1,n), float )
      xu = 10  * ones( (1,n), float )
      print xl, xu, x0

   #  Possibility to modify default values of parameters

   tolmg      = 1e-5      # stopping criterion for gradient norm
   maxeval    = 200       # stopping criterion for nbr. of f evaluations
   whichmodel = 0         # choice of local interpolation models
   hardcons   = 0         # respect bounds as hard constraints
   noisy      = 0         # function is supposed to be noisy
   Deltamin   = 1e-8      # stop at Deltamin
   

   [ x, fx, gx, msg, nit, nf ] = bcdfo( funf, dfunf, x0,
   'lbounds', xl, 'ubounds', xu,
   'epsilon', tolmg,
   'maxeval', maxeval, 'maxit', maxeval, 
   'whichmodel', whichmodel, 
   'hardcons', hardcons, 
   'noisy', noisy , 
   'shrink_Delta', 1, 
   'verbosity', 'medium',
   'restart', 'none',
   'method', 'BCDFO',
   'localSolver', 'CG',
   'Deltamin', Deltamin,
   'geometry_Y0','simplex',
   'degree_m0','linear' )

#   print "\nOutput:"
#   print "x      = ", x
#   print "f      = ", fx
#   print "nit    = ", nit
#   if ( len(gx) == 0 ):
#      print "pgnorm = ", 0
#   else:
#      print "pgnorm = ", linalg.norm(gx)
#   print "nf     = ", nf
#   print ""

