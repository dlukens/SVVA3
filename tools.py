
######### Integration Tool ###########

#for single integration tool, 3/8 simpson rule
#for double integration tool, normal simpson rule
import numpy as np
import matplotlib.pyplot as plt
import interp


def f(x):
    return 1 / x
                    

def integral(f, lower_bound, upper_bound, n=60):
    # f = function
    # n = number of sub-intervals, must be multiple of 3
    if n % 3 >=1:
        raise ValueError("n must be a multiple of 3")
    step_size = (upper_bound - lower_bound)/n
    som = f(lower_bound) + f(upper_bound)
    for i in range(1, n):
        if i % 3 == 0 :
            som = som + 2 * f(lower_bound + i*step_size)
        else:
            som = som + 3 * f(lower_bound + i*step_size)
    R = ((3*step_size)/8)* som
    return (R)

def integral2(f, lower, upper, h=60):

    I = np.zeros(h)
    I2 = np.zeros(h)
    grid = np.linspace(lower, upper, h)
    
    for i in range(len(grid)-1):         
        I[i+1] = I[i] + integral(f, grid[i], grid[i+1])
    
    C = interp.interpolate1d(I, grid)
    
    for j in range(len(grid)-1):
        def f1(x):
            return C[0,j]*(x - grid[j])**3 + C[1,j]*(x - grid[j])**2 + C[2,j]*(x - grid[j]) + C[3,j]
        
        I2[j+1] = I2[j] + integral(f1, grid[j], grid[j+1])
    return I2[-1]


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



#2D and 3D plot graphers

def plot2d(gridx, gridz, data, title):
    #2D contour plot
    
    X, Z = np.meshgrid(gridx, gridz)
    Y = data
    
    plt.figure()
    cp = plt.contourf(X, Z, Y)
    plt.colorbar(cp)
    
    plt.title('{}'.format(title))
    plt.xlabel('X along wingspan [m]')
    plt.ylabel('Z along chord [m]')
    plt.show()
    
def plot3d(gridx, gridz, data, title):
    X, Z = np.meshgrid(gridx, gridz)
    Y = data
    
    #3D plot
    ax = plt.axes(projection='3d')
    
    ax.plot_surface(X, Z, Y, rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none')
    
    plt.title('{}'.format(title))
    ax.set_xlabel('X along wingspan [m]')
    ax.set_ylabel('Z along chord [m]')
    ax.set_zlabel('Distributed load [kN/m^2]')
    ax.view_init(azim=166)
    plt.show()