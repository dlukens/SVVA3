
######### Integration Tool ###########

#for single integration tool, 3/8 simpson rule
#for double integration tool, normal simpson rule

import numpy as np
def f(x):
    return 1 / x
                    

def integral(f, lower_bound, upper_bound, n):
    # f = function
    # n = number of sub=intervals, must be multiple of 3
    if n % 3 >=1:
        raise ValueError("n must be a multiple of 3")
    step_size = (upper_bound - lower_bound)/n
    som = f(lower_bound) + f(upper_bound)
    for i in range(1, n):
        if i % 3 == 0 :
            som = som + 2 * f(lower_bound + i*step_size)
        else:
            som = som + 3 * f(lower_bound + i*step_size)
    R = ((3*h)/8)* som
    return (R)

def f2(x,y):
    return x*y + 2*x + 3*y

def dintegral(f, x_lower_bound, x_upper_bound, y_lower_bound, y_upper_bound, step_size):
    integral = 0
    x = x_lower_bound
    while (x < x_upper_bound):
        y = y_lower_bound
        while (y < y_upper_bound):
            mid_vol = f2(x + 0.5 * step_size, y + 0.5* step_size)* step_size**2
            trap_vol = 0.25 * step_size**2 * (f2(x,y) + f2(x+step_size, y) + f2(x, y + step_size) + f2(x + step_size, y + step_size))
            integral = integral + (2*mid_vol + trap_vol) / 3
            y = y + step_size
        x = x + step_size
    return(integral)

########################################

