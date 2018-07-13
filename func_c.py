import numpy as np

def func_c(xvector):
    # solution on the unitdisk
    ce = np.zeros(3)
    ce[0] = xvector[0]**2 + xvector[1]**2 - 2
    ce[1] = abs(xvector[0]) - 1     # this constraint is not necessary, just added for illustration
    ce[2] = abs(xvector[1]) - 1     # this constraint is not necessary, just added for illustration
    msg = 0
    return ce, msg