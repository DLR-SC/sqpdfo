# -*- coding: utf-8 -*-

from sqpdfo.runtime import *
import sqpdfo.sqpdfo_global_variables as glob
from numpy import array, zeros, concatenate, zeros_like, inf

def f_benchmark(x, prob):
    """
    #-----------------------------------------------------------------------
    # Computation of f, ci, ce
    #-----------------------------------------------------------------------
    """

    if prob == 1:
        f = - (5 - (x[0] - 2) ** 2 - 2 * (x[1] - 1) ** 2)

    elif prob == 2:
        f = 2 * x[0] ** 2 + x[1] ** 2

    elif prob == 3:
        f = x[0] ** 2 + x[1] ** 2 + x[2] ** 2

    elif prob == 4:
        f = x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3]

    elif prob == 5:
        #  Powells function from solnp - manual
        #  x* = (-1.717, 1.5957, 1.8272, -0.7636, -0.7636)

        f = exp_(x[0] * x[1] * x[2] * x[3] * x[4])

    elif prob == 6:
        f = -(0.592 * ((exp_(1) - 1) * x[0]) / ((-0.408 * x[0] + 1) * (exp_(x[0]) - 1)) - 1)

    elif prob == 7:
        # alkyl problem found here :http://www.gamsworld.org/global/globallib/alkyl.htm
        # best known solution found by Baron f*=-1.76499964633
        # x1=-1.76499964633; x2=1.70370291772; x3=1.58470999786; x4=0.543084200389;
        # x5=3.03582206371; x6=2.0; x6=-0.882499871995;
        # x7=0.901319381076; x8=0.95; x8=-17.3201583536; x9=10.4754765591;
        # x10=1.56163636364; x11=1.53535353535;
        # x12=0.99; x12=2.06361745003; x13=0.99; x13=21.9374937958;
        # x14= 1.11111111111; x14=-0.488775780565;
        # x15=0.99; x15=10.7759591879;
        # e1=1.0; e2=-5.13314985138; e3=12.7096102404; e4=0.035;
        # e5=-0.679755732296; e6=-23.0920987324; e7=0.312989497393; e8=-7.01855236578
        f = x[0]

    elif prob == 10:  # problem 19 from Hock-Schittkowskis collection
        f = (x[0] - 10) ** 3 + (x[1] - 20) ** 3

    elif prob == 11:  # problem 21 from Hock-Schittkowskis collection
        f = 0.01 * x[0] ** 2 + x[1] ** 2 - 100

    elif prob == 12:  # problem 35 (Beale) from Hock-Schittkowskis collection
        f = 9.0 - 8 * x[0] - 6 * x[1] - 4 * x[2] + 2 * x[0] ** 2 + 2 * x[1] ** 2 + x[2] ** 2 \
            + 2 * x[0] * x[1] + 2 * x[0] * x[2]

    elif prob == 13:  # problem 76 from Hock-Schittkowskis collection
        f = x[0] ** 2 + 0.5 * x[1] ** 2 + x[2] ** 2 + 0.5 * x[3] ** 2 - x[0] * x[2] \
            + x[2] * x[3] - x[0] - 3 * x[1] + x[2] - x[3]
    elif prob == 14:  # problem 44 from Hock-Schittkowskis collection
        f = x[0] - x[1] - x[2] - x[0] * x[2] + x[0] * x[3] + x[1] * x[2] - x[1] * x[3]

    elif prob == 15:  # 2D Rosenbrock with 2 eq + 1 ineq
        f = (1 - x[0]) ** 2 + 100 * (x[1] - x[0] ** 2) ** 2

    elif prob == 16:  # 2D Rosenbrock with 1 eq + 2 ineq
        f = (1 - x[0]) ** 2 + 100 * (x[1] - x[0] ** 2) ** 2

    elif prob == 1000:
        # CUTEr problems
        cproblem = glob.get_prob_cuter()
        f = cproblem.obj(x.reshape(-1))

    else:
        raise RuntimeError("Unknown Problem number: ", prob)

    return f

def c_benchmark(x, prob):
    """
    #-----------------------------------------------------------------------
    # Computation of f, ci, ce
    #-----------------------------------------------------------------------
    """

    # Initialization
    ce = array([])
    c = array([])

    if prob == 1:
        ce = zeros(1)
        ce = x[0] + 4 * x[1] - 3
        c = ce.reshape(-1, 1)

    elif prob == 2:
        ce = zeros(1)
        ce = x[0] + x[1] - 1
        c = ce.reshape(-1, 1)

    elif prob == 3:
        ce = zeros(2)
        ce[0] = x[0] + x[1] + x[2]
        ce[1] = x[0] + 2 * x[1] + 3 * x[2] - 1
        c = ce.reshape(-1, 1)

    elif prob == 4:
        ce = zeros(3)
        ce[0] = x[0] + x[1] + x[2]
        ce[1] = x[0] + 2 * x[1] + 3 * x[2] - 1
        ce[2] = x[3] ** 3 - 1
        c = ce.reshape(-1, 1)

    elif prob == 5:
        #  Powells function from solnp - manual
        #  x* = (-1.717, 1.5957, 1.8272, -0.7636, -0.7636)

        ce = zeros(3)
        ce[0] = x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2 + x[4] ** 2 - 10
        ce[1] = x[1] * x[2] - 5 * x[3] * x[4]
        ce[2] = x[0] ** 3 + x[1] ** 3 + 1
        c = ce.reshape(-1, 1)

    elif prob == 6:
        pass

    elif prob == 7:
        # alkyl problem found here :http://www.gamsworld.org/global/globallib/alkyl.htm
        # best known solution found by Baron f*=-1.76499964633
        # x1=-1.76499964633; x2=1.70370291772; x3=1.58470999786; x4=0.543084200389;
        # x5=3.03582206371; x6=2.0; x6=-0.882499871995;
        # x7=0.901319381076; x8=0.95; x8=-17.3201583536; x9=10.4754765591;
        # x10=1.56163636364; x11=1.53535353535;
        # x12=0.99; x12=2.06361745003; x13=0.99; x13=21.9374937958;
        # x14= 1.11111111111; x14=-0.488775780565;
        # x15=0.99; x15=10.7759591879;
        # e1=1.0; e2=-5.13314985138; e3=12.7096102404; e4=0.035;
        # e5=-0.679755732296; e6=-23.0920987324; e7=0.312989497393; e8=-7.01855236578
        ce = zeros(8)
        ce[0] = 6.3 * x[4] * x[7] + x[0] - 5.04 * x[1] - 0.35 * x[2] - x[3] - 3.36 * x[5]
        ce[1] = -0.819672131147541 * x[1] + x[4] - 0.819672131147541 * x[5]
        ce[2] = 0.98 * x[3] - x[6] * (0.01 * x[4] * x[9] + x[3])
        ce[3] = -x[1] * x[8] + 10 * x[2] + x[5]
        ce[4] = x[4] * x[11] - x[1] * (1.12 + 0.13167 * x[8] - 0.0067 * x[8] * x[8])
        ce[5] = x[7] * x[12] - 0.01 * (1.098 * x[8] - 0.038 * x[8] * x[8]) - 0.325 * x[6] - 0.57425
        ce[6] = x[9] * x[13] + 22.2 * x[10] - 35.82
        ce[7] = x[10] * x[14] - 3 * x[7] + 1.33
        c = ce.reshape(-1, 1)

    elif prob == 10:  # problem 19 from Hock-Schittkowskis collection
        ci = zeros(2)
        ci[0] = (x[0] - 5) ** 2 + (x[1] - 5) ** 2 - 100
        ci[1] = -(x[1] - 5) ** 2 - (x[0] - 6) ** 2 + 82.81
        c = concatenate((ce, ci))
        c = c.reshape(-1, 1)

    elif prob == 11:  # problem 21 from Hock-Schittkowskis collection
        ci = zeros(1)
        ci[0] = 10 * x[0] - x[1] - 10  # - x[2]
        c = concatenate((ce, ci))
        c = c.reshape(-1, 1)

    elif prob == 12:  # problem 35 (Beale) from Hock-Schittkowskis collection
        ci = zeros(1)
        ci[0] = 3 - x[0] - x[1] - 2 * x[2]
        c = concatenate((ce, ci))
        c = c.reshape(-1, 1)

    elif prob == 13:  # problem 76 from Hock-Schittkowskis collection
        ci = zeros(3)
        ci[0] = 5 - x[0] - 2 * x[1] - x[2] - x[3]
        ci[1] = 4 - 3 * x[0] - x[1] - 2 * x[2] + x[3]
        ci[2] = x[1] + 4 * x[2] - 1.5
        c = concatenate((ce, ci))
        c = c.reshape(-1, 1)

    elif prob == 14:  # problem 44 from Hock-Schittkowskis collection
        ci = zeros(6)
        ci[0] = 8 - x[0] - 2 * x[1]
        ci[1] = 12 - 4 * x[0] - x[1]
        ci[2] = 12 - 3 * x[0] - 4 * x[1]
        ci[3] = 8 - 2 * x[2] - x[3]
        ci[4] = 8 - x[2] - 2 * x[3]
        ci[5] = 5 - x[2] - x[3]
        c = concatenate((ce, ci))
        c = c.reshape(-1, 1)

    elif prob == 15:  # 2D Rosenbrock with 2 eq + 1 ineq
        ce = np.zeros(2)
        ce[0] = x[0] ** 2 + x[1] ** 2 - 2
        ce[1] = - (x[0] - 1) ** 3 + x[1] - 1
        ci = np.zeros(1)
        ci[0] = - x[0] - x[1] + 2
        ce = ce.reshape(-1, 1)
        ci = ci.reshape(-1, 1)
        c = concatenate((ce, ci))

    elif prob == 16:  # 2D Rosenbrock with 1 eq + 2 ineq
        ce = np.zeros(1)
        ce[0] = x[0] ** 2 + x[1] ** 2 - 2
        ci = np.zeros(2)
        ci[0] = - x[0] - x[1] + 2
        ci[1] = - (x[0] - 1) ** 3 + x[1] - 1
        ce = ce.reshape(-1, 1)
        ci = ci.reshape(-1, 1)
        c = concatenate((ce, ci))

    elif prob == 1000:
        # CUTEr problems
        cproblem = glob.get_prob_cuter()
        (_, c) = cproblem.objcons(x.reshape(-1))

        if cproblem.m > 0:
            me = sum(cproblem.is_eq_cons)
            mi = cproblem.m - me
            li = cproblem.cl
            ui = cproblem.cu
            ce_new = []
            ci_new = []
            cnew = []

            # re-order c such that ce first and then ci
            if mi > 0:

                for i in range(0, cproblem.m):
                    if li[i] == ui[i]:  # equalities
                        ce_new.append(c[i] - li[i])
                        # print('eq')

                    else:  # inequalities
                        if li[i] == -1e20 and ui[i] == 0.0:
                            ci_new.append(-c[i])
                            # print('ineq to switch')
                        elif li[i] == -1e20 and ui[i] < 1e7:
                            ci_new.append(-c[i] + ui[i])
                            # print('ineq to switch and to change')
                        elif li[i] == 0.0 and ui[i] == 1e20:
                            ci_new.append(c[i])
                            # print('ineq good bounds')
                        elif li[i] > -1e7 and ui[i] == 1e20:
                            ci_new.append(c[i] - li[i])
                            # print('ineq to change')
                        else:
                            # Handling of two-sided inequalities !!!!
                            # print('ineq two-sided')
                            # print(li[i],ui[i])
                            if li[i] > -1e7:
                                ci_new.append(c[i] - li[i])
                            if ui[i] < 1e7:
                                ci_new.append(-c[i] + ui[i])

                cnew = concatenate((ce_new, ci_new))
                c = cnew.reshape(-1, 1)
            else:
                c = c.reshape(-1, 1)
                if sum(li) > 0 or sum(ui) > 0:
                    print('sqpdfo_func: Warning! ce must not be zero! Check li and ui!')

    else:
        raise RuntimeError("Unknown Problem number: ", prob)

    return c



def benchmark_start_values(prob):
    """
    # This function returns the dimensions of the problem:
# . n  = number of variables,
# . nb = number of variables with bounds,
# . mi = number of inequality constraints,
# . me = number of equality constraints.
    """

    # Set output variables

    x0 = array([])
    lx = array([])
    ux = array([])
    li = None
    ui = None



    # dxmin = sqrt(eps);

    if prob == 1:
        n = 2
        nb = 2
        mi = 0
        me = 1
        x0 = array([[4.6], [0.0]]).T
        lx = array([1.95, - 1e+20]).reshape(-1, 1)
        ux = array([1e+20, 0.3]).reshape(-1, 1)
    elif prob == 2:
        n = 2
        nb = 0
        mi = 0
        me = 1
        x0 = array([[- 1], [2.54378]]).T
        lx = - inf * ones_(n, 1)
        ux = inf * ones_(n, 1)
    elif prob == 3:
        nb = 2
        mi = 0
        me = 2
        x0 = array([[0.0], [0.0], [0.5]]).T
        n = length_(x0)
        lx = array([- 0.5, 0.0, - inf]).reshape(-1, 1)
        ux = array([inf, inf, inf]).reshape(-1, 1)
    elif prob == 4:
        nb = 0
        mi = 0
        me = 3
        x0 = array([[1.0], [1.0], [1.0], [0.0]]).T
        n = length_(x0)
        lx = - inf * ones_(n, 1)
        ux = inf * ones_(n, 1)
    elif prob == 5:
        nb = 0
        mi = 0
        me = 3
        x0 = array([[- 2.0], [2.0], [2.0], [1.0], [1.0]]).T
        n = 5
        lx = - inf * ones_(n, 1)
        ux = inf * ones_(n, 1)
    elif prob == 6:
        n = 1
        nb = 1
        mi = 0
        me = 0
        x0 = array([[0.6]])
        lx = array([0.5]).reshape(-1, 1)
        ux = array([0.8]).reshape(-1, 1)
    elif prob == 7:  # alkyl problem found here :http://www.gamsworld.org/global/globallib/alkyl.htm
        n = 15
        nb = 14
        me = 8
        mi = 0
        x0 = array([[-0.9, 1.745, 1.2, 1.1, 3.048, 1.974, 0.893, 0.928, 8, 3.6, 1.50, 1, 1, 1, 1]]).T
        lx = array([-inf, 0, 0, 0, 0, 0, 0.85, 0.9, 3, 1.2, 1.45, 0.99, 0.99, 0.9, 0.99]).reshape(-1, 1)
        ux = array(
            [inf, 2, 1.6, 1.2, 5, 2, 0.93, 0.95, 12, 4, 1.62, 1.01010101010101, 1.01010101010101, 1.11111111111111,
             1.01010101010101]).reshape(-1, 1)
    elif prob == 10:  # problem 19 from Hock-Schittkowskis collection
        n = 2
        nb = 4
        me = 0
        mi = 2
        x0 = array([[20.1, 5.84]])
        lx = array([[13.0, 0.0]]).reshape(-1, 1)
        ux = array([[100.0, 100.0]]).reshape(-1, 1)
    elif prob == 11:  # problem 21 from Hock-Schittkowskis collection
        n = 2
        nb = 4
        me = 0
        mi = 1
        x0 = array([[-1.0, -1.0]])
        lx = array([[2.0, -50.0]]).reshape(-1, 1)
        ux = array([[50.0, 50.0]]).reshape(-1, 1)
    elif prob == 12:  # problem 35 (Beale) from HS collection
        n = 3
        nb = 3
        me = 0
        mi = 1
        x0 = array([[0.5, 0.5, 0.5]])
        lx = array([[0.0, 0.0, 0.0]]).reshape(-1, 1)
        ux = array([[1e20, 1e20, 1e20]]).reshape(-1, 1)
    elif prob == 13:  # problem 76 from Hock-Schittkowskis collection
        n = 4
        nb = 4
        me = 0
        mi = 3
        x0 = array([[0.5, 0.5, 0.5, 0.5]])
        lx = array([[0.0, 0.0, 0.0, 0.0]]).reshape(-1, 1)
        ux = array([[1e20, 1e20, 1e20, 1e20]]).reshape(-1, 1)
    elif prob == 14:  # problem 44 from Hock-Schittkowskis collection
        n = 4
        nb = 4
        me = 0
        mi = 6
        x0 = array([[0.0, 0.0, 0.0, 0.0]])
        lx = array([[0.0, 0.0, 0.0, 0.0]]).reshape(-1, 1)
        ux = array([[1e20, 1e20, 1e20, 1e20]]).reshape(-1, 1)
    elif prob == 15:
        n = 2
        nb = 4
        me = 2
        mi = 1
        x0 = array([[-1.2, 1.0]])
        lx = array([[-5.0, -5.0]]).reshape(-1, 1)
        ux = array([[10.0, 10.0]]).reshape(-1, 1)
    elif prob == 16:
        n = 2
        nb = 4
        me = 1
        mi = 2
        x0 = array([[0.3, 0.3]])
        lx = array([[-5.0, -5.0]]).reshape(-1, 1)
        ux = array([[10.0, 10.0]]).reshape(-1, 1)
    elif prob == 1000:
        # Warning : here the CUTEst interface from this website has to be
        # installed in order to use CUTEst problems :
        # https://jfowkes.github.io/pycutest/_build/html/index.html
        # Thanks to Jaroslav Fowkes and Lindon Roberts

        cproblem = glob.get_prob_cuter()
        n = cproblem.n
        m = cproblem.m
        me = sum(cproblem.is_eq_cons)
        mi = m - me
        x0 = cproblem.x0.reshape(-1, 1)
        lx = cproblem.bl.reshape(-1, 1)
        ux = cproblem.bu.reshape(-1, 1)
        li = cproblem.cl
        ui = cproblem.cu
        nb = sum_(min_((lx[0:n] > -inf) + (inf > ux[0:n]), 1))
        # print(cproblem.eq_cons_first)
    else:
        raise RuntimeError("Unknown Problem number: ", prob)
    return x0, lx, ux, li, ui, n, nb, mi, me


def get(prob):
    """
    Returns the benchmark problem, including function, constraint function,
    bounds ...
    """
    def f_func(x):
        return f_benchmark(x, prob)

    def c_func(x):
        return c_benchmark(x, prob)

    x0, lx, ux, li, ui, n, nb, mi, me = benchmark_start_values(prob)

    return f_func, x0, lx, ux, me, mi, c_func, li, ui


def set_test_prob(prob):
    f_func, x0, lx, ux, me, mi, c_func, li, ui = get(prob)

    glob.set_filename_f(f_func)
    glob.set_filename_cons(c_func)
