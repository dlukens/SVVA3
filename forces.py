import numpy as np
import tools as t
from interpolationA03 import C0_z, C0_x, grid_z, grid_x, data
import interpolationA03 as itr
import matplotlib.pyplot as plt

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



#interpolate first integral

C1_z = itr.interpolate(I_z, grid_z[:-1])
C1_x = itr.interpolate(I_x, grid_x[:-1])


#Second integration

II_z, II_zsum = integrate2(C1_z, grid_z)
II_x, II_xsum = integrate2(C1_x, grid_x)


##cp location
Cp_x = -II_xsum/I_xsum[:-1]
Cp_coeff = itr.interpolate1d(Cp_x[1::2], grid_x[:-2]) #Grid is halved to adjust for cpx datapoints


# itr.interplot(0.001, grid_x[:-2], Cp_function)


###Plotting###
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
    plt.plot(grid_x[:-2], Cp_x[1::2])
    plt.show()
    

# plot2d(grid_x, grid_z, data, 'Distributed force q')
# plot2d(grid_x[:-1], grid_z[:-1], I_x, 'First integral')
# plot2d(grid_x[:-1], grid_z[:-2], II_x, 'Second Integral')

# plot2d(grid_x, grid_z, data, 'Distributed force q')
# plot2d(grid_z[:-1], grid_x[:-1], I_z, 'First integral')
# plot2d(grid_z[:-1], grid_x[:-2], II_z, 'Second Integral')


