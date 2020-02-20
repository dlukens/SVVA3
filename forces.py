import numpy as np
import tools as t
from interpolationA03 import C0_z, C0_x, grid_z, grid_x

#Integrating once the interpolated force function
#S(z) = a(z-z_i)^3 + b(z-z_i)^2 + c(z-z_i) + d

def func_z(z):
    return C0_z[0,0,0]*(z - grid_z[0])**3 + C0_z[0,1,0]*(z - grid_z[0])**2 + C0_z[0,2,0]*(z - grid_z[0]) + C0_z[0,3,0]

    
def integrate2(C0, grid):
    n = 60
    T = np.zeros(len(grid)-1)
    
    for i in range(len(grid)-1):
        T[i] = t.integral(func_z, grid[i], grid[i+1], n)

    return T

T = integrate2(C0_z, grid_z)