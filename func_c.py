import numpy as np

def func_c(xvector):
    # solution on the unitdisk
    ce = np.zeros(1)
    ce[0] = xvector[0]**2 + xvector[1]**2 - 2;
    ce = ce.reshape(-1,1)
    msg = 0
    return ce, msg