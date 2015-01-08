#!/usr/local/bin/python
from numpy import *
import helper

@helper.convertingDecorator
def bcdfo_projgrad_(n,x,g,bl,bu):
	return bcdfo_projgrad(n,x,g,bl,bu)


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
      print "g[0]\n", g[0]				
      gi = g[i][0]#g[0][i]
      if ( gi < 0.0 ):
         gn = append(gn,-min(abs(bu[i]-x[i]),-gi))
      else:

         print "x", x
         print "bl", bl
         print "bl[i]", bl[i]									
         print "x[i]", x[i]									
         print "abs(bl[i]-x[i])", abs(bl[i]-x[i])
         print "gi", gi
         print "gn", gn									
         print "min(abs(bl[i]-x[i])", append(gn,min(abs(bl[i]-x[i]),gi))									
         gn = append(gn,min(abs(bl[i]-x[i]),gi))
      
      # infinity-norm
      gnorm = max(gnorm,abs(gn[i]))

   return gnorm, gn
