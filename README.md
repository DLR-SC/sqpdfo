
SQPDFO Optimizer
================

[![Build](https://github.com/dlr-sc/sqpdfo/workflows/sqpdfo/badge.svg)](https://github.com/DLR-SC/sqpdfo/actions?query=workflow%3Asqpdfo+branch%3Amaster)

This is the SQPDFO (Sequential-Quadratic-Programming Derivative-Free Optimization) 
Optimizer for generally  constrained nonlinear optimization without derivatives,
developed by A. Troeltzsch at German Aerospace Center (DLR).
See license information in LICENSE.

SQPDFO is a further development of ECDFO and it uses parts of BCDFO.

-------------------------------------------------------------------------------

To run the optimizer you will need the Python packages numpy and scipy:

>pip install numpy
>pip install scipy

-------------------------------------------------------------------------------

Then you launch SQPDFO by:

>python run_sqpdfo.py 

or 

>python run_sqpdfo_academic_examples.py

-------------------------------------------------------------------------------

To add your own optimization problems, there are two ways:

(1)
- modify the function func_f() in run_sqpdfo.py for the objective function
- modify the function func_c() in run_sqpdfo.py for the constraint functions
- modify starting point and bounds information in run_sqpdfo.py
- launch >python run_sqpdfo.py

(2)
- add the problem by following other examples in sqpdfo_init_prob.py and sqpdfo_func.py
- give it a new problem number 
- change set_prob(number_of_the_new_problem) in run_sqpdfo_academic_examples.py
- launch >python run_sqpdfo_academic_examples.py

-------------------------------------------------------------------------------

If you want to use SQPDFO to solve problems from the CUTEst library:
(tested only on Linux/Unix)

- download and install CUTEst from https://github.com/ralna/CUTEst  
- download pycutest with eg. >pip install pycutest
- set all environment variables from CUTEst and pycutest
- launch >python run_sqpdfo_cutest.py

-------------------------------------------------------------------------------

If this code is of any use to you, please cite the following papers:

Tröltzsch, Anke (2016) 
A Sequential Quadratic Programming Algorithm for Equality-Constrained Optimization 
without Derivatives. Optimization Letters, 10 (2), pages 383-399. Springer. 
DOI: 10.1007/s11590-014-0830-y ISSN 1862-4472 

Serge Gratton, Philippe L. Toint, and Anke Tröltzsch (2011) 
An active-set trust-region method for derivative-free nonlinear bound-constrained 
optimization. Optimization Methods and Software, 26 (4-5), pages 873-894. 
DOI: 10.1080/10556788.2010.549231 

-------------------------------------------------------------------------------
