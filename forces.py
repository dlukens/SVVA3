import numpy as np
import tools as t
from interpolationA03 import C0_z, C0_x, grid_z, grid_x

#Integrating once the interpolated force function
#S(z) = a(z-z_i)^3 + b(z-z_i)^2 + c(z-z_i) + d


print(C0_z[0,0,1])

def integrate2(C0, grid):
    n = 120
    I = np.zeros((len(C0[:,:,0]), len(grid)-1))
    
    for i in range(len(C0[:,:,0])):
        for j in range(len(grid)-1):
            
            def func(z):
                return C0[i,0,j]*(z - grid[j])**3 + C0[i,1,j]*(z - grid[j])**2 + C0[i,2,j]*(z - grid[j]) + C0[i,3,j]
            
            I[i, j] = t.integral(func, grid[j], grid[j+1], n)

    return I

I_z = integrate2(C0_z, grid_z)
# I_x = integrate2(C0_x, grid_x)