import numpy as np
import tools as t
from interpolationA03 import C0_z, C0_x, grid_z, grid_x
import interpolationA03 as itr

#Integrating once the interpolated force function
#S(z) = a(z-z_i)^3 + b(z-z_i)^2 + c(z-z_i) + d



def integrate2(C0, grid):
    n = 60
    I = np.zeros((len(C0[:,:,0]), len(grid)-1))
    I_sum = np.zeros(len(grid))
    
    for i in range(len(C0[:,:,0])):
        for j in range(len(C0[0,0,:])):
            
            def func(z):
                return C0[i,0,j]*(z - grid[j])**3 + C0[i,1,j]*(z - grid[j])**2 + C0[i,2,j]*(z - grid[j]) + C0[i,3,j]
            
            I[i, j] = abs(t.integral(func, grid[j], grid[j+1], n))
            
    I_sum = I.sum(axis=1)

    return I, I_sum

I_z, I_zsum = integrate2(C0_z, grid_z)
I_x, I_xsum = integrate2(C0_x, grid_x)

#interpolate first integralx

C1_z = itr.interpolate(I_z, grid_z[:-1])
C1_x = itr.interpolate(I_x, grid_x[:-1])


###Second integration

II_z, II_zsum = integrate2(C1_z, grid_z)
II_x, II_xsum = integrate2(C1_x, grid_x)


##cp location
Cp_z = -II_xsum/I_xsum[:-1]



