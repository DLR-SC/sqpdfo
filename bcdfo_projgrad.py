#!/usr/local/bin/python
from numpy import *

def bcdfo_projgrad(n,x,g,bl,bu):
########################################################
#  This function computes the projected gradient and its
#  infinity norm.
#
#  INPUTS:
#
#  n     : dimension
#  x     : current iterate
#  g     : gradient at x
#  bl    : lower bounds
#  bu    : upper bounds
#
#  OUTPUTS:
#
#  gnorm : infinity norm of the projected gradient
#  gn    : projected gradient vector
#
########################################################

   gnorm = 0.0
   gn = array([[]])

   for i in range(0,n):
      gi = g[0][i]
      if ( gi < 0.0 ):
         gn = append(gn,-min(abs(bu[i]-x[i]),-gi))
      else:
         gn = append(gn,min(abs(bl[i]-x[i]),gi))
      
      # infinity-norm
      gnorm = max(gnorm,abs(gn[i]))

   return gnorm, gn
