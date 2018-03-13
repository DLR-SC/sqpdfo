def func_f(xvector):
    # 2D Rosenbrock function (constrained on the unitdisk if func_c.py is considered)
    # x* = (1.0, 1.0)        (constrained on the unitdisk: x* = (0.7864, 0.6177))
    # f* = 0.0               (constrained on the unitdisk: f* = 0.045674824758137236)
    f =  (1-xvector[0])**2 + 100*(xvector[1]-xvector[0]**2)**2
    msg = 0
    return f,msg