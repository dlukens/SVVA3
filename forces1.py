import numpy as np
import tools as t
from interpolationA03 import C0_z, C0_x, grid_z, grid_x, data
import interpolationA03 as itr
import matplotlib.pyplot as plt

#Integrating once the interpolated force function
#S(z) = a(z-z_i)^3 + b(z-z_i)^2 + c(z-z_i) + d


def polyintegrate(C0):
    C1 = np.zeros((len(C0[:,0,0]), 4, len(C0[0,0,:])))
    #S(z) = a(z-z_i)^3 + b(z-z_i)^2 + c(z-z_i) + d
    #to
    #S(z) = 1/4*a(z-z_i)^4 + 1/3*b(z-z_i)^3 + 1/2*c(z-z_i)^2 + dz
    for i in range(len(C0[:,0,0])):
    
        C1[i,0,:] = 1/4*C0[i,0,:]
        C1[i,1,:] = 1/3*C0[i,1,:]
        C1[i,2,:] = 1/2*C0[i,2,:]
        C1[i,3,:] = C0[i,3,:]
        
    return(C1)


def integrate(C0a, grid):
    I = np.zeros((len(C0a[:,:,0]), len(grid)-1))
    I_sum = np.zeros(len(grid))
    
    for i in range(len(C0a[:,:,0])):
        for j in range(len(C0a[0,0,:])):
            
            def f(z):
                return 1/4*C0a[i,0,j]*(z - grid[j])**4 + 1/3*C0a[i,1,j]*(z - grid[j])**3 + 1/2*C0a[i,2,j]*(z - grid[j])**2 + C0a[i,3,j]*z
            
            I[i, j] = abs(f(grid[j+1]) - f(grid[j]))
            
    I_sum = I.sum(axis=1)

    return I, I_sum

def integrate2(C1a, grid):
    I = np.zeros((len(C1a[:,:,0]), len(grid)-1))
    I_sum = np.zeros(len(grid))
    
    for i in range(len(C1a[:,:,0])):
        for j in range(len(C1a[0,0,:])):
            
            def f(z):
                return 1/5*C1a[i,0,j]*(z - grid[j])**5 + 1/4*C1a[i,1,j]*(z - grid[j])**4 + 1/3*C1a[i,2,j]*(z - grid[j])**3 + 1/2*C1a[i,3,j]*z**2
            
            I[i, j] = abs(f(grid[j+1]) - f(grid[j]))
            
    I_sum = I.sum(axis=1)

    return I, I_sum

Ia_z, Ia_zsum = integrate(C0_z, grid_z)
Ia_x, Ia_xsum = integrate(C0_x, grid_x)
    

#Coefficients from first integral

C1a_z = polyintegrate(C0_z)
C1a_x = polyintegrate(C0_x)

itr.interplot(0.001, grid_x, C1a_x[0,:,:])


#Second integration

IIa_z, IIa_zsum = integrate2(C1a_z, grid_z)
IIa_x, IIa_xsum = integrate2(C1a_x, grid_x)


##cp location
Cpa_x = -IIa_zsum/Ia_zsum
# Cp_coeff = itr.interpolate1d(Cp_x, grid_x[:-1]) #Grid is halved to adjust for cpx datapoints


# itr.interplot(0.001, grid_x[:-2], Cp_coeff)


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
    plt.plot(grid_x[:-1], Cpa_x)
    plt.show()
    

plot2d(grid_x, grid_z, data, 'Distributed force q')
# plot2d(grid_x[:-1], grid_z[:-1], I_x, 'First integral')
# plot2d(grid_x[:-2], grid_z[:-2], II_x, 'Second Integral')

# plot2d(grid_x, grid_z, data, 'Distributed force q')
# plot2d(grid_z[:-1], grid_x[:-1], I_z, 'First integral')
# plot2d(grid_z[:-2], grid_x[:-2], II_z, 'Second Integral')


