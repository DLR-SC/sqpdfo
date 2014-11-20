#! /usr/bin/env python
from numpy import linalg

def bcdfo_checkZ( Z, kappa_ill ):
   """
   Check condition number of the interpolation matrix Z
   """
   
   badcond = 0

   try:
      condZ = linalg.cond(Z)      
      if (condZ > kappa_ill):
         badcond = 1
   except:
      badcond = 1

   return badcond